#!/usr/bin/env python
u"""
convert_harmonics.py
Written by Tyler Sutterley (09/2020)
Converts a file from the spatial domain into the spherical harmonic domain

CALLING SEQUENCE:
    python convert_harmonics.py -F 2 --lmax 60 -U 1 infile outfile

COMMAND LINE OPTIONS:
    --help: list the command line options
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -r X, --reference X: Reference frame for load love numbers
    -U X, --units X: input units
        1: cm of water thickness (cmwe)
        2: Gigatonnes (Gt)
        3: mm of water thickness kg/m^2
    -S X, --spacing X: spatial resolution of input data (dlon,dlat)
    -I X, --interval X: input grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
    --missing: input spatial fields have missing values
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
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    ncdf_read_stokes.py: reads spherical harmonic netcdf files
    ncdf_stokes.py: writes output spherical harmonic data to netcdf
    hdf5_read_stokes.py: reads spherical harmonic HDF5 files
    hdf5_stokes.py: writes output spherical harmonic data to HDF5
    spatial.py: spatial data class for reading, writing and processing data
    ncdf_read.py: reads input spatial data from netCDF4 files
    hdf5_read.py: reads input spatial data from HDF5 files
    ncdf_write.py: writes output spatial data to netCDF4
    hdf5_write.py: writes output spatial data to HDF5
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 09/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: updates to reading load love numbers
    Written 10/2019
"""
from __future__ import print_function

import sys
import os
import re
import argparse
import numpy as np

from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.utilities import get_data_path

#-- PURPOSE: read load love numbers for the range of spherical harmonic degrees
def load_love_numbers(LMAX, REFERENCE='CF'):
    #-- load love numbers file
    love_numbers_file = get_data_path(['data','love_numbers'])
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    if (LMAX > 696):
        #-- Creates arrays of kl, hl, and ll Love Numbers
        hl = np.zeros((LMAX+1))
        kl = np.zeros((LMAX+1))
        ll = np.zeros((LMAX+1))
        hl[:697],kl[:697],ll[:697] = read_love_numbers(love_numbers_file,
            FORMAT='tuple', REFERENCE=REFERENCE)
        #-- for degrees greater than 696
        for l in range(697,LMAX+1):
            hl[l] = 2.0*hl[l-1] - hl[l-2]#-- linearly extrapolating hl
            kl[l] = 2.0*kl[l-1] - kl[l-2]#-- linearly extrapolating kl
            ll[l] = 2.0*ll[l-1] - ll[l-2]#-- linearly extrapolating ll
    else:
        #-- read arrays of kl, hl, and ll Love Numbers
        hl,kl,ll=read_love_numbers(love_numbers_file, FORMAT='tuple',
            REFERENCE=REFERENCE)
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- PURPOSE: converts from the spatial domain into the spherical harmonic domain
def convert_harmonics(INPUT_FILE, OUTPUT_FILE, LMAX=None, MMAX=None, UNITS=None,
    REFERENCE=None, DDEG=None, INTERVAL=None, MISSING=False, FILL_VALUE=None,
    HEADER=None, DATAFORM=None, VERBOSE=False, MODE=0o775):

    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- Grid spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    #-- Grid dimensions
    if (INTERVAL == 1):#-- (0:360, 90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
    elif (INTERVAL == 2):#-- degree spacing/2
        nlon = np.int((360.0/dlon))
        nlat = np.int((180.0/dlat))

    #-- read spatial file in data format
    #-- expand dimensions
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        input_spatial = spatial(spacing=[dlon,dlat],nlat=nlat,
            nlon=nlon).from_ascii(INPUT_FILE,header=HEADER).expand_dims()
    elif (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        input_spatial = spatial().from_netCDF4(INPUT_FILE).expand_dims()
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        input_spatial = spatial().from_HDF5(INPUT_FILE).expand_dims()
    #-- convert missing values to zero
    input_spatial.replace_invalid(0.0)
    #-- input data shape
    nlat,nlon,nt = input_spatial.shape

    #-- read arrays of kl, hl, and ll Love Numbers
    LOVE = load_love_numbers(LMAX, REFERENCE=REFERENCE)

    #-- calculate associated Legendre polynomials
    th = (90.0 - input_spatial.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(th))
    #-- date count array
    counter = np.arange(nt)

    #-- allocate for output spherical harmonics
    Ylms = harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.time = input_spatial.time.copy()
    Ylms.month = np.array(12.0*(Ylms.time - 2002.0),dtype='i') + 1
    Ylms.clm = np.zeros((LMAX+1,LMAX+1,nt))
    Ylms.slm = np.zeros((LMAX+1,LMAX+1,nt))
    for t in range(nt):
        #-- convert spatial field to spherical harmonics
        output_Ylms = gen_stokes(input_spatial.data[:,:,t].T,
            input_spatial.lon, input_spatial.lat, UNITS=UNITS,
            LMIN=0, LMAX=LMAX, MMAX=MMAX, PLM=PLM, LOVE=LOVE)
        Ylms.clm[:,:,t] = output_Ylms['clm'][:,:].copy()
        Ylms.slm[:,:,t] = output_Ylms['slm'][:,:].copy()

    #-- if verbose output: print input and output file names
    if VERBOSE:
        print('{0}:'.format(os.path.basename(sys.argv[0])))
        print('{0} -->\n\t{1}'.format(INPUT_FILE,OUTPUT_FILE))
    #-- outputting data to file
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        Ylms.to_ascii(OUTPUT_FILE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4 (.nc)
        Ylms.to_netCDF4(OUTPUT_FILE)
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        Ylms.to_HDF5(OUTPUT_FILE)
    #-- change output permissions level to MODE
    os.chmod(OUTPUT_FILE,MODE)

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Converts a file from the spatial domain into the
            spherical harmonic domain
            """
    )
    #-- command line parameters
    #-- input and output file
    parser.add_argument('infile',
        type=os.path.expanduser, nargs='?',
        help='Input spatial file')
    parser.add_argument('outfile',
        type=os.path.expanduser, nargs='?',
        help='Output harmonic file')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference','-r',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3],
        help='Output units')
    #-- output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help='Output grid interval (1: global, 2: centered global)')
    #-- fill value
    parser.add_argument('--missing',
        default=False, action='store_true',
        help='Input spatial fields have missing values')
    parser.add_argument('--fill-value','-f',
        type=float,
        help='Set fill_value for input spatial fields')
    #-- ascii parameters
    parser.add_argument('--header',
        type=int,
        help='Number of header rows to skip in input ascii files')
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the output files (octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- run program with parameters
    convert_harmonics(args.infile, args.outfile, LMAX=args.lmax, MMAX=args.mmax,
        UNITS=args.units, DDEG=args.spacing, INTERVAL=args.interval,
        MISSING=args.missing, FILL_VALUE=args.fill_value, HEADER=args.header,
        DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
