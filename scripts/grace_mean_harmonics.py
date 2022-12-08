#!/usr/bin/env python
u"""
grace_mean_harmonics.py
Written by Tyler Sutterley (12/2022)

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
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        SLR: satellite laser ranging coefficients from CSR
        UCI: Sutterley and Velicogna coefficients, Remote Sensing (2019)
        Swenson: GRACE-derived coefficients from Sean Swenson
        GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    --slr-c20 X: Replace C20 coefficients with SLR values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GFZ: use values from GFZ
        GSFC: use values from GSFC (TN-14)
    --slr-21 X: Replace C21 and S21 coefficients with SLR values
        CSR: use values from CSR
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC
    --slr-22 X: Replace C22 and S22 coefficients with SLR values
        CSR: use values from CSR
    --slr-c30 X: Replace C30 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC (TN-14)
    --slr-c40 X: Replace C40 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --slr-c50 X: Replace C50 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --mean-file X: Output GRACE/GRACE-FO mean file
    --mean-format X: Output data format for GRACE/GRACE-FO mean file
        ascii
        netCDF4
        HDF5
        gfc
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
        Replaces low-degree harmonics with SLR values (if specified)
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 04/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
        option to specify a specific geocenter correction file
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: simplified file exports using wrappers in harmonics
        added option to output in gravity field coefficients (gfc) format
        remove choices for argparse processing centers
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
import logging
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: import GRACE/GRACE-FO files for a given months range
# calculate the mean of the spherical harmonics and output to file
def grace_mean_harmonics(base_dir, PROC, DREL, DSET, LMAX,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    ATM=False,
    POLE_TIDE=False,
    DEG1=None,
    DEG1_FILE=None,
    MODEL_DEG1=False,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C40=None,
    SLR_C50=None,
    MEAN_FILE=None,
    MEANFORM=None,
    VERBOSE=0,
    MODE=0o775):

    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''

    # data formats for output: ascii, netCDF4, HDF5, gfc
    suffix = dict(ascii='txt',netCDF4='nc',HDF5='H5',gfc='gfc')[MEANFORM]

    # reading GRACE months for input date range
    # replacing low-degree harmonics with SLR values if specified
    # include degree 1 (geocenter) harmonics if specified
    # correcting for Pole Tide Drift and Atmospheric Jumps if specified
    input_Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX, SLR_21=SLR_21,
        SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C40=SLR_C40, SLR_C50=SLR_C50,
        DEG1_FILE=DEG1_FILE, MODEL_DEG1=MODEL_DEG1, ATM=ATM,
        POLE_TIDE=POLE_TIDE)
    grace_Ylms = gravtk.harmonics().from_dict(input_Ylms)
    # descriptor string for processing parameters
    grace_str = input_Ylms['title']
    # calculate mean Ylms
    mean_Ylms = mean().from_harmonics(grace_Ylms.mean())
    mean_Ylms.time = np.mean(grace_Ylms.time)
    mean_Ylms.month = np.mean(grace_Ylms.month)
    # number of months
    nt = grace_Ylms.shape[-1]
    # calculate RMS of harmonic errors
    mean_Ylms.eclm = np.sqrt(np.sum(input_Ylms['eclm']**2,axis=2)/nt)
    mean_Ylms.eslm = np.sqrt(np.sum(input_Ylms['eslm']**2,axis=2)/nt)
    # product information
    mean_Ylms.center = PROC
    mean_Ylms.release = DREL
    mean_Ylms.product = DSET

    # default output filename if not entering via parameter file
    if not MEAN_FILE:
        DIRECTORY = os.path.expanduser(input_Ylms['directory'])
        args = (PROC,DREL,DSET,grace_str,LMAX,order_str,START,END,suffix)
        file_format = '{0}_{1}_{2}_MEAN_CLM{3}_L{4:d}{5}_{6:03d}-{7:03d}.{8}'
        MEAN_FILE = os.path.join(DIRECTORY,file_format.format(*args))
    else:
        DIRECTORY = os.path.dirname(MEAN_FILE)
    # recursively create output directory if non-existent
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY, MODE)

    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # output spherical harmonics for the static field
    if (MEANFORM == 'gfc'):
        # output mean field to gfc format
        mean_Ylms.to_gfc(MEAN_FILE, verbose=VERBOSE)
    else:
        # output mean field to specified file format
        mean_Ylms.to_file(MEAN_FILE, format=MEANFORM,
            verbose=VERBOSE, **attributes)
    # change the permissions mode
    os.chmod(MEAN_FILE, MODE)

    # return the output file
    return MEAN_FILE

# PURPOSE: additional routines for the harmonics module
class mean(gravtk.harmonics):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.center=None
        self.release='RLxx'
        self.product=None
        self.eclm=None
        self.eslm=None

    def from_harmonics(self, temp):
        """
        Convert a harmonics object to a new mean object
        """
        self = mean(lmax=temp.lmax, mmax=temp.mmax)
        # try to assign variables to self
        for key in ['clm','slm','eclm','eslm','shape','ndim','filename',
            'center','release','product']:
            try:
                val = getattr(temp, key)
                setattr(self, key, np.copy(val))
            except AttributeError:
                pass
        # assign ndim and shape attributes
        self.update_dimensions()
        return self

    def to_gfc(self, filename, **kwargs):
        """
        Write a harmonics object to gfc file
        Inputs: full path of output gfc file
        Options:
            harmonics objects contain date information
            keyword arguments for gfc output
        """
        self.filename = os.path.expanduser(filename)
        # set default verbosity
        kwargs.setdefault('verbose',False)
        logging.info(self.filename)
        # open the output file
        fid = open(self.filename, mode='w', encoding='utf8')
        # print the header informat
        self.print_header(fid)
        # output file format
        file_format = ('{0:3} {1:4d} {2:4d} {3:+18.12E} {4:+18.12E}  '
            '{5:11.5E}  {6:11.5E}')
        # write to file for each spherical harmonic degree and order
        for m in range(0, self.mmax+1):
            for l in range(m, self.lmax+1):
                args = ('gfc', l, m, self.clm[l,m], self.slm[l,m],
                    self.eclm[l,m], self.eslm[l,m])
                print(file_format.format(*args), file=fid)
        # close the output file
        fid.close()

    # PURPOSE: print gfc header to top of file
    def print_header(self, fid):
        # print header
        fid.write('{0} {1}\n'.format('begin_of_head',73*'='))
        fid.write('{0:30}{1}\n'.format('product_type','gravity_field'))
        fid.write('{0:30}{1}\n'.format('center',self.center))
        fid.write('{0:30}{1}\n'.format('release',self.release))
        fid.write('{0:30}{1}\n'.format('product',self.product))
        fid.write('{0:30}{1:+16.10E}\n'.format('earth_gravity_constant',
            3.986004415E+14))
        fid.write('{0:30}{1:+16.10E}\n'.format('radius',6.378136300E+06))
        fid.write('{0:30}{1:d}\n'.format('max_degree',self.lmax))
        fid.write('{0:30}{1}\n'.format('errors','uncalibrated'))
        fid.write('{0:30}{1}\n'.format('norm','fully_normalized'))
        args = ('key','L','M','C','S','sigma C','sigma S')
        fid.write('\n{0:7}{1:5}{2:10}{3:20}{4:15}{5:13}{6:7}\n'.format(*args))
        fid.write('{0} {1}\n'.format('end_of_head',75*'='))

# PURPOSE: print a file log for the GRACE/GRACE-FO mean program
def output_log_file(input_arguments, output_file):
    # format: GRACE_mean_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILE:')
    logging.info(output_file)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the GRACE/GRACE-FO mean program
def output_error_log_file(input_arguments):
    # format: GRACE_mean_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
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
        description="""Calculates the temporal mean of the GRACE/GRACE-FO
            spherical harmonics
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # Data processing center or satellite mission
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
    # use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    # correct for pole tide drift follow Wahr et al. (2015)
    parser.add_argument('--pole-tide',
        default=False, action='store_true',
        help='Correct for pole tide drift')
    # Update Degree 1 coefficients with SLR or derived values
    # Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
    #     https://grace.jpl.nasa.gov/data/get-data/geocenter/
    # SLR: satellite laser ranging from CSR
    #     ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
    # UCI: Sutterley and Velicogna, Remote Sensing (2019)
    #     https://www.mdpi.com/2072-4292/11/18/2108
    # Swenson: GRACE-derived coefficients from Sean Swenson
    #     https://doi.org/10.1029/2007JB005338
    # GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    #     http://gravis.gfz-potsdam.de/corrections
    parser.add_argument('--geocenter',
        metavar='DEG1', type=str,
        choices=['Tellus','SLR','SLF','UCI','Swenson','GFZ'],
        help='Update Degree 1 coefficients with SLR or derived values')
    parser.add_argument('--geocenter-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Specific geocenter file if not default')
    parser.add_argument('--interpolate-geocenter',
        default=False, action='store_true',
        help='Least-squares model missing Degree 1 coefficients')
    # replace low degree harmonics with values from Satellite Laser Ranging
    parser.add_argument('--slr-c20',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C20 coefficients with SLR values')
    parser.add_argument('--slr-21',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C21 and S21 coefficients with SLR values')
    parser.add_argument('--slr-22',
        type=str, default=None, choices=['CSR','GSFC'],
        help='Replace C22 and S22 coefficients with SLR values')
    parser.add_argument('--slr-c30',
        type=str, default=None, choices=['CSR','GFZ','GSFC','LARES'],
        help='Replace C30 coefficients with SLR values')
    parser.add_argument('--slr-c40',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C40 coefficients with SLR values')
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Output GRACE/GRACE-FO mean file')
    # input data format (ascii, netCDF4, HDF5, gfc)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Output data format for GRACE/GRACE-FO mean file')
    # Output log file for each job in forms
    # GRACE_mean_run_2002-04-01_PID-00000.log
    # GRACE_mean_failed_run_2002-04-01_PID-00000.log
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run grace_mean_harmonics algorithm with parameters
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
            DEG1_FILE=args.geocenter_file,
            MODEL_DEG1=args.interpolate_geocenter,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C40=args.slr_c40,
            SLR_C50=args.slr_c50,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as e:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_file)

# run main program
if __name__ == '__main__':
    main()
