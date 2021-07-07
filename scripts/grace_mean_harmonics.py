#!/usr/bin/env python
u"""
grace_mean_harmonics.py
Written by Tyler Sutterley (06/2021)

Calculates the temporal mean of the GRACE/GRACE-FO spherical harmonics
    for a given date range from a set of parameters

COMMAND LINE OPTIONS:
    -D X, --directory X: GRACE/GRACE-FO working data directory
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    --atm-correction: Apply atmospheric jump correction coefficients
    --pole-tide: Correct for pole tide drift
    --geocenter X: Update Degree 1 coefficients with SLR or derived values
    --slr-c20 X: Replace C20 coefficients with SLR values
    --slr-21 X: Replace C21 and S21 coefficients with SLR values
    --slr-22 X: Replace C22 and S22 coefficients with SLR values
    --slr-c30 X: Replace C30 coefficients with SLR values
    --slr-c50 X: Replace C50 coefficients with SLR values
    --mean-file X: Output GRACE/GRACE-FO mean file
    --mean-format X: Output data format for GRACE/GRACE-FO mean file
    --log: Output a log file listing output files
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces C20,C21,S21,C22,S22,C30 and C50 with SLR values (if specified)
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: using the harmonics class for spherical harmonic operations
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: can replace C30 with coefficients from SLR
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 03/2018: edited output netCDF4 and HDF5 attributes
    Updated 04/2017: changed no input file exception to IOError
    Updated 05/2016: using __future__ print function
    Updated 02/2016: use getopt parameters to set options for the
        number of processes to run in parallel and to output a log file
    Updated 12/2015: output unique log file
    Updated 08/2015: changed sys.exit to a raise exception instance
        Added pole tide and GAE/GAF/GAG correction parameters
    Updated 05/2015: added MMAX parameter to study new 2015 60X30 fields
    Updated 01/2015: added error handling for multiprocessing threads
    Updated 11/2014: added HDF5 dataform option
    Updated 10/2014: updated for distributed computing of tasks
        with the multiprocessing module
    Written 05/2014
"""
from __future__ import print_function

import sys
import os
import time
import argparse
import numpy as np
import traceback
from gravity_toolkit.grace_input_months import grace_input_months
from gravity_toolkit.harmonics import harmonics
import gravity_toolkit.utilities as utilities

#-- PURPOSE: keep track of threads
def info(args):
    print(os.path.basename(sys.argv[0]))
    print(args)
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: import GRACE/GRACE-FO files for a given months range
#-- calculate the mean of the spherical harmonics and output to file
def grace_mean_harmonics(base_dir, PROC, DREL, DSET, LMAX,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    ATM=False,
    POLE_TIDE=False,
    DEG1=None,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C50=None,
    MEAN_FILE=None,
    MEANFORM=None,
    VERBOSE=False,
    MODE=0o775):

    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''

    #-- data formats for output: ascii, netCDF4, HDF5
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[MEANFORM]

    #-- reading GRACE months for input range with grace_input_months.py
    #-- replacing low-degree harmonics with SLR values if specified
    #-- include degree 1 (geocenter) harmonics if specified
    #-- correcting for Pole Tide Drift and Atmospheric Jumps if specified
    input_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX,
        SLR_21=SLR_21, SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C50=SLR_C50,
        MODEL_DEG1=False, POLE_TIDE=POLE_TIDE, ATM=ATM)
    grace_Ylms = harmonics().from_dict(input_Ylms)
    #-- descriptor string for processing parameters
    grace_str = input_Ylms['title']
    #-- calculate mean Ylms
    mean_Ylms = grace_Ylms.mean()
    mean_Ylms.time = np.mean(grace_Ylms.time)
    mean_Ylms.month = np.mean(grace_Ylms.month)

    #-- default output filename if not entering via parameter file
    if not MEAN_FILE:
        DIRECTORY = os.path.expanduser(input_Ylms['directory'])
        args = (PROC,DREL,DSET,grace_str,LMAX,order_str,START,END,suffix)
        file_format = '{0}_{1}_{2}_MEAN_CLM{3}_L{4:d}{5}_{6:03d}-{7:03d}.{8}'
        MEAN_FILE = os.path.join(DIRECTORY,file_format.format(*args))
    else:
        DIRECTORY = os.path.dirname(MEAN_FILE)
    #-- recursively create output directory if non-existent
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY, MODE)

    #-- output spherical harmonics for the static field
    if (MEANFORM == 'ascii'):
        #-- output mean field to ascii
        mean_Ylms.to_ascii(MEAN_FILE,verbose=VERBOSE)
    elif (MEANFORM == 'netCDF4'):
        #-- output mean field to netCDF4
        mean_Ylms.to_netCDF4(MEAN_FILE,verbose=VERBOSE)
    elif (MEANFORM == 'HDF5'):
        #-- output mean field to HDF5
        mean_Ylms.to_HDF5(MEAN_FILE,verbose=VERBOSE)
    #-- change the permissions mode
    os.chmod(MEAN_FILE, MODE)

    #-- return the output file
    return MEAN_FILE

#-- PURPOSE: print a file log for the GRACE/GRACE-FO mean program
def output_log_file(arguments,output_file):
    #-- format: GRACE_mean_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_run_{0}_PID-{1:d}.log'.format(*args)
    #-- create a unique log and open the log file
    DIRECTORY = os.path.expanduser(arguments.directory)
    fid = utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    #-- print argument values sorted alphabetically
    print('ARGUMENTS:', file=fid)
    for arg, value in sorted(vars(arguments).items()):
        print('{0}: {1}'.format(arg, value), file=fid)
    #-- print output files
    print('\n\nOUTPUT FILE:',file=fid)
    print('{0}'.format(output_file),file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE/GRACE-FO mean program
def output_error_log_file(arguments):
    #-- format: GRACE_mean_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_failed_run_{0}_PID-{1:d}.log'.format(*args)
    #-- create a unique log and open the log file
    DIRECTORY = os.path.expanduser(arguments.directory)
    fid = utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    #-- print argument values sorted alphabetically
    print('ARGUMENTS:', file=fid)
    for arg, value in sorted(vars(arguments).items()):
        print('{0}: {1}'.format(arg, value), file=fid)
    #-- print traceback error
    print('\n\nTRACEBACK ERROR:', file=fid)
    traceback.print_exc(file=fid)
    #-- close the log file
    fid.close()

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates the temporal mean of the GRACE/GRACE-FO
            spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, default='CSR',
        choices=['CSR','GFZ','JPL','CNES'],
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    #-- GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    #-- start and end GRACE/GRACE-FO months
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
    #-- use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    #-- correct for pole tide drift follow Wahr et al. (2015)
    parser.add_argument('--pole-tide',
        default=False, action='store_true',
        help='Correct for pole tide drift')
    #-- Update Degree 1 coefficients with SLR or derived values
    #-- Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
    #--     https://grace.jpl.nasa.gov/data/get-data/geocenter/
    #-- SLR: satellite laser ranging from CSR
    #--     ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
    #-- SLF: Sutterley and Velicogna, Remote Sensing (2019)
    #--     https://www.mdpi.com/2072-4292/11/18/2108
    #-- Swenson: GRACE-derived coefficients from Sean Swenson
    #--     https://doi.org/10.1029/2007JB005338
    #-- GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    #--     http://gravis.gfz-potsdam.de/corrections
    parser.add_argument('--geocenter',
        metavar='DEG1', type=str,
        choices=['Tellus','SLR','SLF','Swenson','GFZ'],
        help='Update Degree 1 coefficients with SLR or derived values')
    #-- replace low degree harmonics with values from Satellite Laser Ranging
    parser.add_argument('--slr-c20',
        type=str, default='GSFC', choices=['CSR','GFZ','GSFC'],
        help='Replace C20 coefficients with SLR values')
    parser.add_argument('--slr-21',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C21 and S21 coefficients with SLR values')
    parser.add_argument('--slr-22',
        type=str, default=None, choices=['CSR'],
        help='Replace C22 and S22 coefficients with SLR values')
    parser.add_argument('--slr-c30',
        type=str, default='GSFC', choices=['CSR','GFZ','GSFC','LARES'],
        help='Replace C30 coefficients with SLR values')
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    #-- mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Output GRACE/GRACE-FO mean file')
    #-- input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Output data format for GRACE/GRACE-FO mean file')
    #-- Output log file for each job in forms
    #-- GRACE_mean_run_2002-04-01_PID-00000.log
    #-- GRACE_mean_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args,_ = parser.parse_known_args()

    #-- try to run the analysis with listed parameters
    try:
        info(args) if args.verbose else None
        #-- run grace_mean_harmonics algorithm with parameters
        output_file = grace_mean_harmonics(
            args.directory,
            args.center,
            args.release,
            args.product,
            args.lmax,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            ATM=args.atm_correction,
            POLE_TIDE=args.pole_tide,
            DEG1=args.geocenter,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C50=args.slr_c50,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
        if args.log:#-- write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:#-- write successful job completion log file
            output_log_file(args,output_file)

#-- run main program
if __name__ == '__main__':
    main()
