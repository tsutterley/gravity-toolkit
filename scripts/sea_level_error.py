#!/usr/bin/env python
u"""
sea_level_error.py
Written by Tyler Sutterley (05/2023)
Reads in sea level grid error files and converts to spherical
    harmonics for use in the least_squares_mascons.py program

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for mascon files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --runs X: Number of Monte Carlo iterations
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gen_stokes.py: Computes geoid Stokes coefficients for an input grid
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: updated inputs to spatial from_ascii function
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 07/2020: using spatial data class for input and output operations
        compute Legendre Polynomials using land mask latitudes
    Updated 04/2020: updates to reading load love numbers
    Updated 10/2019: changing Y/N flags to True/False
    Updated 12/2018: added parallel processing with multiprocessing
    Updated 06/2018: using python3 compatible octal and input
    Written 04/2018
"""
from __future__ import print_function

import sys
import os
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np

import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# program module to run with specified parameters
def sea_level_error(PROC, DREL, DSET, LMAX,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    RUNS=0,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = f'M{MMAX:d}' if MMAX and (MMAX != LMAX) else ''
    # distributing mascon mass uniformly over ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF4 file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # calculate Fully-Normalized Legendre Polynomials
    th = (90.0 - landsea.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # create index file for least_squares_mascons.py
    args = (ITERATION,dset_str,ocean_str,LMAX,order_str)
    INDEX = 'SLF_MC_ITERATION_{0}_INDEX{1}{2}_CLM_L{3:d}{4}.txt'.format(*args)
    index_file = OUTPUT_DIRECTORY.joinpath(INDEX)
    fid = index_file.open(mode='w', encoding='utf8')
    # print the path to the index file
    logging.info(str(index_file))

    # input and output file format
    file_format = 'SLF_MC_ITERATION_{0}{1}{2}{3}_L{4:d}{5}_{6:05d}.{7}'
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # for each iteration
    output_files = []
    for n in range(0, RUNS):
        # sea level file for iteration (spatial fields)
        SLF = file_format.format(ITERATION,dset_str,ocean_str,'',
            EXPANSION,'',n,suffix[DATAFORM])
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(SLF)
        # read sea level file
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            dinput = gravtk.spatial().from_ascii(INPUT_FILE,
                date=False, spacing=[dlon,dlat], nlat=nlat, nlon=nlon)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            dinput = gravtk.spatial().from_netCDF4(INPUT_FILE, date=False)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            dinput = gravtk.spatial().from_HDF5(INPUT_FILE, date=False)

        # Converting sea level field into spherical harmonics (can truncate)
        Ylms = gravtk.gen_stokes(dinput.data.T, dinput.lon, dinput.lat,
            UNITS=1, LMIN=0, LMAX=LMAX, MMAX=MMAX, LOVE=LOVE, PLM=PLM)

        # output (truncated) spherical harmonics to file
        FILE = file_format.format(ITERATION,dset_str,ocean_str,'_CLM',
            LMAX,order_str,n,suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            Ylms.to_ascii(OUTPUT_FILE, date=False)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            Ylms.to_netCDF4(OUTPUT_FILE, date=False, **attributes)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            Ylms.to_HDF5(OUTPUT_FILE, date=False, **attributes)
        # change output file permissions mode to MODE
        OUTPUT_FILE.chmod(mode=MODE)
        # add to output list
        output_files.append(OUTPUT_FILE)
        # add to output index file
        print(Ylms.compressuser(OUTPUT_FILE), file=fid)
    # close the index file
    fid.close()
    # change the permissions mode of the output index file
    index_file.chmod(mode=MODE)
    # return the list of output files
    return output_files

# PURPOSE: print a file log for the sea level harmonics calculation
def output_log_file(input_arguments, output_files):
    # format: sea_level_error_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()),os.getpid())
    LOGFILE = 'sea_level_error_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the sea level harmonics calculation
def output_error_log_file(input_arguments):
    # format: failed_sea_level_error_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()),os.getpid())
    LOGFILE = 'failed_sea_level_error_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads in sea level grid error files and converts
            to spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # sea level fingerprint parameters
    parser.add_argument('--iteration','-I',
        type=int, default=1,
        help='Sea level fingerprint iteration')
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # number of monte carlo iterations
    parser.add_argument('--runs',
        type=int, default=10000,
        help='Number of Monte Carlo iterations')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # sea_level_error_run_2002-04-01_PID-00000.log
    # sea_level_error_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run sea_level_error algorithm with parameters
        output_files = sea_level_error(
            args.center,
            args.release,
            args.product,
            args.lmax,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            RUNS=args.runs,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
