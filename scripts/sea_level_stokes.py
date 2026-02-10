#!/usr/bin/env python
u"""
sea_level_stokes.py
Written by Tyler Sutterley (05/2023)
Reads in sea level grid files and converts to spherical harmonics

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for mascon files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
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
    -F X, --format X: input data format for auxiliary files
        ascii
        netCDF4
        HDF5
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        AW13-IJ05: Geruo A IJ05-R2 GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
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
    Updated 01/2023: refactored associated legendre polynomials
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
    Updated 07/2020: compute Legendre Polynomials using land mask latitudes
    Updated 06/2020: using spatial data class for input and output operations
    Updated 04/2020: updates to reading load love numbers
    Updated 10/2019: changing Y/N flags to True/False
    Updated 12/2018: added parallel processing with multiprocessing
    Updated 11/2018: can vary the land-sea mask for ascii files
    Updated 06/2018: using python3 compatible octal and input
    Updated 04/2018: changed the names of the output log files
    Updated 03/2018: added option --mode to set the output file permissions
    Updated 02/2018: added parameter EXPANSION to specify the SLF truncation
    Updated 08/2017: running at 0.5x0.5 degrees
    Written 08/2017
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
def sea_level_stokes(PROC, DREL, DSET, LMAX,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    GIA=None,
    GIA_FILE=None,
    DATAFORM=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # GRACE/GRACE-FO months
    months = sorted(set(np.arange(START,END+1)) - set(MISSING))
    nmon = len(months)
    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = f'M{MMAX:d}' if MMAX and (MMAX != LMAX) else ''
    # input GIA stokes coefficients to get titles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    gia_str = f'_{GIA_Ylms_rate.title.upper()}' if GIA else ''
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
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

    # create index file for calc_mascon.py
    # (spherical harmonics to be removed from the GRACE data)
    args = (ITERATION,dset_str,gia_str,ocean_str,LMAX,order_str,START,END)
    INDEX = ('SLF_ITERATION_{0}_INDEX{1}{2}{3}_CLM_L{4:d}{5}_'
        '{6:03d}-{7:03d}.txt').format(*args)
    index_file = OUTPUT_DIRECTORY.joinpath(INDEX)
    fid = index_file.open(mode='w', encoding='utf8')
    # print the path to the index file
    logging.info(index_file)

    # for each month
    output_files = []
    # input and output format
    file_format = 'SLF_ITERATION_{0}{1}{2}{3}_{4}L{5:d}{6}_{7:03d}.{8}'
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    for t in range(0, nmon):
        # sea level file for month (spatial fields)
        SLF = file_format.format(ITERATION, dset_str, gia_str, ocean_str, '',
            EXPANSION, '', months[t], suffix[DATAFORM])
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(SLF)
        # read sea level file
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            dinput = gravtk.spatial().from_ascii(INPUT_FILE,
                spacing=[dlon,dlat], nlat=nlat, nlon=nlon)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            dinput = gravtk.spatial().from_netCDF4(INPUT_FILE)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            dinput = gravtk.spatial().from_HDF5(INPUT_FILE)

        # Converting sea level field into spherical harmonics (can truncate)
        Ylms = gravtk.gen_stokes(dinput.data.T, dinput.lon, dinput.lat,
            UNITS=1, LMIN=0, LMAX=LMAX, MMAX=MMAX, LOVE=LOVE, PLM=PLM)
        Ylms.time = np.copy(dinput.time)
        Ylms.month = months[t]

        # output (truncated) spherical harmonics to file
        FILE = file_format.format(ITERATION, dset_str, gia_str, ocean_str,
            'CLM_', LMAX, order_str, months[t], suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            Ylms.to_ascii(OUTPUT_FILE)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            Ylms.to_netCDF4(OUTPUT_FILE, **attributes)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            Ylms.to_HDF5(OUTPUT_FILE, **attributes)
        # change output file permissions mode to MODE
        OUTPUT_FILE.chmod(mode=MODE)
        # add to output list
        output_files.append(OUTPUT_FILE)
        # add to output index file
        print(gravtk.harmonics().compressuser(OUTPUT_FILE), file=fid)
    # close the index file
    fid.close()
    # change the permissions mode of the output index file
    index_file.chmod(mode=MODE)
    # return the list of output files
    return output_files

# PURPOSE: print a file log for the sea level harmonics calculation
def output_log_file(input_arguments, output_files):
    # format: sea_level_stokes_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'sea_level_stokes_run_{0}_PID-{1:d}.log'.format(*args)
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
    # format: failed_sea_level_stokes_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'failed_sea_level_stokes_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Reads in sea level grid files and converts
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
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
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
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
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
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # sea_level_stokes_run_2002-04-01_PID-00000.log
    # failed_sea_level_stokes_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
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
        # run sea_level_stokes algorithm with parameters
        output_files = sea_level_stokes(
            args.center,
            args.release,
            args.product,
            args.lmax,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            DATAFORM=args.format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
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
