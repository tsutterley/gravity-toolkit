#!/usr/bin/env python
u"""
regional_spherical_caps.py
Written by Tyler Sutterley (05/2023)

Computes and outputs spherical harmonics for set of regional spherical caps
    from a file listing the coordinates of each spherical cap center

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for mascon files
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
    --mascon-file X: index file of mascons spherical harmonics
    --coordinate-file X: file with spatial coordinates of mascon centers
    --header X: number of header lines in coordinate file
    --cap-radius X: spherical cap radius (degrees)
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gen_spherical_cap.py: Computes geoid Stokes coefficients for spherical cap
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use file not found error for coordinate file
    Updated 01/2021: harmonics object output from gen_spherical_cap.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: using harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2017: save coordinate file as NetCDF/HDF5 title attribute
    Updated 04/2017: changed no input file exception to IOError
    Updated 06/2016: using __future__ print function, num variable set to int
    Updated 03/2016: added output_files for log files with unique filenames
        use getopt parameters to set number of PROCESSES to run in parallel,
            whether or not to output a log file, added new help module
    Updated 08/2015: changed sys.exit to a raise exception instance
    Updated 06/2015: can have variable radii using 4th column of file
    Updated 04/2015: fix for INTERVAL 2 on plm_global
    Updated 01/2015: complete update of program
        added main definition for parameter files
        added multiprocessing for multiple parameter files
        with traceback error handling
    Written 10/2013
"""
from __future__ import print_function

import sys
import os
import time
import logging
import pathlib
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate spherical caps from a coordinate index
def regional_spherical_caps(LMAX,
    MASCON_FILE=None,
    COORDINATE_FILE=None,
    HEADER=0,
    RAD_CAP=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    OUTPUT_DIRECTORY=None,
    VERBOSE=0,
    MODE=0o775):

    # check if coordinate file exists
    COORDINATE_FILE = pathlib.Path(COORDINATE_FILE).expanduser().absolute()
    if not COORDINATE_FILE.exists():
        raise FileNotFoundError(str(COORDINATE_FILE))

    # recursively create output directories if not existent
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output file format (ascii, netCDF4, HDF5)
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # list of output files
    output_files = []

    # input coordinate file (center points of each spherical cap)
    # read coordinate file for lat/lon of spherical cap centers
    coord = np.loadtxt(COORDINATE_FILE, skiprows=HEADER)
    # column 1: cap number
    # column 2: longitude of center point
    # column 3: latitude of center point
    num = coord[:,0].astype(int)
    lon = coord[:,1]
    lat = coord[:,2]
    ncap = len(num)
    # radius of each spherical cap
    if RAD_CAP:
        # radius of all spherical caps from an argument
        RAD_CAP = np.zeros((ncap)) + RAD_CAP
    else:
        # column 4: radius of spherical cap
        RAD_CAP = coord[:,3]

    # equivalent water thickness of each spherical cap
    dinput = 1.0

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # harmonic units
    dfactor = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    # Earth Parameters
    rho_e = dfactor.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = dfactor.rad_e# Average Radius of the Earth [cm]

    # colatitude of the spherical cap centers
    th = (90.0 - lat)*np.pi/180.0
    # Legendre polynomials of the spherical caps
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # index file listing spherical caps used in region
    # open index file for spherical harmonic file names
    MASCON_FILE = pathlib.Path(MASCON_FILE).expanduser().absolute()
    fid = MASCON_FILE.open(mode='w', encoding='utf8')
    # append output files to lists
    output_files.append(MASCON_FILE)

    # for each spherical cap:
    # calculate the spherical harmonic coefficients
    # write spherical harmonics to file
    # print file path of spherical harmonics to index file
    for i,cap_number in enumerate(num):
        # plms for computing spherical cap at latitude i
        plm_i = np.squeeze(PLM[:,:,i])
        # Calculate spherical harmonic coefficients
        Ylms = gravtk.gen_spherical_cap(dinput, lon[i], lat[i], LMAX=LMAX,
            RAD_CAP=RAD_CAP[i], UNITS=1, PLM=plm_i, LOVE=LOVE)

        # calculate equivalent area after harmonic conversion (1 cmwe)
        # area = (volume*density)/(mass/area)
        ar = 4.0*np.pi*(rad_e**3.0)*rho_e*np.squeeze(Ylms.clm[0,0])/3.0
        # calculate equivalent radius (should equal RAD_CAP)
        rad = np.sqrt(ar/np.pi)/rad_e*180.0/np.pi
        # if verbose output: sanity check of radii
        args = (cap_number, RAD_CAP[i], rad)
        logging.info('{0:4d} {1:10.4f} {2:10.4f}'.format(*args))

        # output spherical harmonics file
        arg = (RAD_CAP[i], num[i], LMAX, suffix[DATAFORM])
        Ylms_file = 'SPH_CAP_RAD{0:0.1f}_{1:d}_L{2:d}.{3}'.format(*arg)
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(Ylms_file)

        # print output file of spherical harmonics to index
        # replace full path with tilde
        print(gravtk.harmonics().compressuser(OUTPUT_FILE), file=fid)

        # append output files to lists
        output_files.append(OUTPUT_FILE)

        # attributes for output files
        attributes = {}
        attributes['title'] = str(RAD_CAP[i])
        attributes['reference'] = COORDINATE_FILE.name
        # output spherical harmonic file to file format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            Ylms.to_ascii(OUTPUT_FILE, date=False, verbose=VERBOSE)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            Ylms.to_netCDF4(OUTPUT_FILE, date=False, verbose=VERBOSE,
                **attributes)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            Ylms.to_HDF5(OUTPUT_FILE, date=False, verbose=VERBOSE,
                **attributes)
        # change the permissions mode of the output file
        OUTPUT_FILE.chmod(mode=MODE)

    # close the index file
    fid.close()
    # change the permissions mode of the index file
    MASCON_FILE.chmod(mode=MODE)

    # return list of output files
    return output_files

# PURPOSE: print a file log for the regional spherical cap calculation
def output_log_file(input_arguments, output_files):
    # format: regional_spherical_cap_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE='regional_spherical_cap_run_{0}_PID-{1:d}.log'.format(*args)
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

# PURPOSE: print a error file log for the regional spherical cap calculation
def output_error_log_file(input_arguments):
    # format: regional_spherical_cap_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE='regional_spherical_cap_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Computes and outputs spherical harmonics for a
            set of spherical caps
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=pathlib.Path,
        required=True,
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--coordinate-file',
        type=pathlib.Path,
        required=True,
        help='File with spatial coordinates of mascon centers')
    # number of header lines to skip in coordinate file
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines to skip in coordinate file')
    # spherical cap radius (if using a uniform set of caps)
    parser.add_argument('--cap-radius',
        type=float,
        help='Spherical cap radius (degrees)')
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
    # Output log file for each job in forms
    # regional_spherical_cap_run_2002-04-01_PID-00000.log
    # regional_spherical_cap_failed_run_2002-04-01_PID-00000.log
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
        # run regional_spherical_caps algorithm with parameters
        output_files = regional_spherical_caps(args.lmax,
            MASCON_FILE=args.mascon_file,
            COORDINATE_FILE=args.coordinate_file,
            HEADER=args.header,
            RAD_CAP=args.cap_radius,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            OUTPUT_DIRECTORY=args.output_directory,
            VERBOSE=args.verbose,
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
