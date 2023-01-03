#!/usr/bin/env python
u"""
convert_harmonics.py
Written by Tyler Sutterley (01/2023)
Converts a file from the spatial domain into the spherical harmonic domain

CALLING SEQUENCE:
    python convert_harmonics.py -F 2 --lmax 60 -U 1 infile outfile

COMMAND LINE OPTIONS:
    --help: list the command line options
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -U X, --units X: input units
        1: cm of water thickness (cmwe)
        2: Gigatonnes (Gt)
        3: mm of water thickness kg/m^2
    -S X, --spacing X: spatial resolution of input data (dlon,dlat)
    -I X, --interval X: input grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
    -f X, --fill-value X: set fill_value for input spatial fields
    --header X: number of header rows to skip in input ascii files
    -F X, --format X: input and output data format
        ascii
        netCDF4
        HDF5
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
        iterate over spatial objects versus indexing
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: fix to use fill values for input ascii files
        use functions for converting to and from GRACE months
    Updated 08/2021: fix spherical harmonic orders if not set
    Updated 06/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: updates to reading load love numbers
    Written 10/2019
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: converts from the spatial domain into the spherical harmonic domain
def convert_harmonics(INPUT_FILE, OUTPUT_FILE,
    LMAX=None,
    MMAX=None,
    UNITS=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DDEG=None,
    INTERVAL=None,
    FILL_VALUE=None,
    HEADER=None,
    DATAFORM=None,
    MODE=0o775):

    # verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    # Grid spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    # Grid dimensions
    if (INTERVAL == 1):# (0:360, 90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):# degree spacing/2
        nlon = np.int64((360.0/dlon))
        nlat = np.int64((180.0/dlat))

    # read spatial file in data format
    # expand dimensions
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        input_spatial = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,
            nlon=nlon,fill_value=FILL_VALUE).from_ascii(INPUT_FILE,
            header=HEADER).expand_dims()
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        input_spatial = gravtk.spatial().from_netCDF4(INPUT_FILE).expand_dims()
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        input_spatial = gravtk.spatial().from_HDF5(INPUT_FILE).expand_dims()
    # convert missing values to zero
    input_spatial.replace_invalid(0.0)
    # input data shape
    nlat,nlon,nt = input_spatial.shape

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # calculate associated Legendre polynomials
    th = (90.0 - input_spatial.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # create list of harmonics objects
    Ylms_list = []
    for i,spatial_data in enumerate(input_spatial):
        # convert spatial field to spherical harmonics
        output_Ylms = gravtk.gen_stokes(spatial_data.data.T,
            spatial_data.lon, spatial_data.lat, UNITS=UNITS,
            LMIN=0, LMAX=LMAX, MMAX=MMAX, PLM=PLM, LOVE=LOVE)
        output_Ylms.time = np.copy(spatial_data.time)
        output_Ylms.month = gravtk.time.calendar_to_grace(spatial_data.time)
        # append to list
        Ylms_list.append(output_Ylms)
    # convert Ylms list for output spherical harmonics
    Ylms = gravtk.harmonics().from_list(Ylms_list, clear=True)

    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # outputting data to file
    Ylms.to_file(OUTPUT_FILE, format=DATAFORM, **attributes)
    # change output permissions level to MODE
    os.chmod(OUTPUT_FILE, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts a file from the spatial domain into the
            spherical harmonic domain
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input spatial file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output harmonic file')
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
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3],
        help='Output units')
    # output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2],
        help='Input grid interval (1: global, 2: centered global)')
    # fill value for ascii
    parser.add_argument('--fill-value','-f',
        type=float,
        help='Set fill_value for input spatial fields')
    # ascii parameters
    parser.add_argument('--header',
        type=int,
        help='Number of header rows to skip in input ascii files')
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the output files (octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program with parameters
    try:
        info(args)
        convert_harmonics(args.infile, args.outfile,
            LMAX=args.lmax,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            FILL_VALUE=args.fill_value,
            HEADER=args.header,
            DATAFORM=args.format,
            MODE=args.mode)
    except Exception as e:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
