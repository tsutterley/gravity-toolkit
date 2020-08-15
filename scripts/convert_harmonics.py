#!/usr/bin/env python
u"""
convert_harmonics.py
Written by Tyler Sutterley (08/2020)
Converts a file from the spatial domain into the spherical harmonic domain

CALLING SEQUENCE:
    python convert_harmonics.py -F 2 --lmax=60 -U 1 input_file output_file

COMMAND LINE OPTIONS:
    --help: list the command line options
    --lmax=X: maximum spherical harmonic degree
    --mmax=X: maximum spherical harmonic order
    --reference=X: Reference frame for load love numbers
    -U X, --units=X: input units
        1: cm of water thickness (cmwe)
        2: Gigatonnes (Gt)
        3: mm of water thickness kg/m^2
    -S X, --spacing=X: spatial resolution of input data (dlon,dlat)
    -I X, --interval=X: input grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
    --missing: input spatial fields have missing values
    --fill-value=X: set fill_value for input spatial fields
    --header=X: number of header rows to skip in input ascii files
    --delimiter=X: delimiter in input ascii files
    -F X, --format=X: input and output data format
        1: ascii format
        2: netCDF4 format
        3: HDF5 format
    -V, --verbose: verbose output of processing run
    -M X, --mode=X: Permissions mode of the files created

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
    ncdf_read.py: reads spatial data from netCDF4 files
    hdf5_read.py: reads spatial data from HDF5 files
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: updates to reading load love numbers
    Written 10/2019
"""
from __future__ import print_function

import sys
import os
import re
import getopt
import numpy as np

from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.hdf5_read import hdf5_read
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
    HEADER=None, DELIMITER=None, DATAFORM=None, VERBOSE=False, MODE=0o775):

    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- python dictionary with input spatial variables
    input_spatial = {}

    #-- Grid spacing
    if (np.ndim(DDEG) == 0): #-- dlon=dlat
        dlon = DDEG
        dlat = DDEG
    else: #-- dlon ne dlat
        dlon = DDEG[0]
        dlat = DDEG[1]

    #-- Grid dimensions
    if (INTERVAL == 1):#-- (0:360, 90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
        input_spatial['lon'] = dlon*np.arange(0,nlon)
        input_spatial['lat'] = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):#-- degree spacing/2
        nlon = np.int((360.0/dlon))
        nlat = np.int((180.0/dlat))
        input_spatial['lon'] = np.arange(dlon/2.0,360+dlon/2.0,dlon)
        input_spatial['lat'] = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)

    #-- read input spherical harmonic coefficients from file in DATAFORM
    if (DATAFORM == 1):
        file_input = np.loadtxt(INPUT_FILE,skiprows=HEADER,delimiter=DELIMITER)
        input_spatial['data'] = np.zeros((1,nlat,nlon))
        input_spatial['time'] = 1.0
        #-- number of files
        file_lines = len(file_input)
        #-- for each file line
        for j,line in enumerate(file_input):
            #-- calculating the lon/lat indice
            ilon = np.int(line[0]/dlon)
            ilat = np.int((90.0-line[1])/dlat)
            input_spatial['data'][0,ilat,ilon] = line[2]
        #-- value for missing data points
        missing = FILL_VALUE if FILL_VALUE else 0.0
    elif DATAFORM in (2,3):
        if (DATAFORM == 2):
            #-- read input netCDF4 file (.nc)
            input_spatial = ncdf_read(INPUT_FILE, DATE=True, MISSING=MISSING)
        elif (DATAFORM == 3):
            #-- read input HDF5 file (.H5)
            input_spatial = hdf5_read(INPUT_FILE, DATE=True, MISSING=MISSING)
        #-- extract time
        tdec = input_spatial['time']
        #-- value for missing data points
        if FILL_VALUE:
            missing = FILL_VALUE
        elif MISSING:
            missing = input_spatial['attributes']['_FillValue']
        else:
            missing = 0.0

    #-- convert missing values to zero
    ii,jj,kk = np.nonzero(input_spatial['data'] == missing)
    input_spatial['data'][ii,jj,kk] = 0.0
    nt,nlat,nlon = np.shape(input_spatial['data'])

    #-- read arrays of kl, hl, and ll Love Numbers
    LOVE = load_love_numbers(LMAX, REFERENCE=REFERENCE)

    #-- calculate associated Legendre polynomials
    th = (90.0 - input_spatial['lat'])*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(th))
    #-- date count array
    counter = np.arange(nt)

    #-- allocate for output spherical harmonics
    Ylms = harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.time = input_spatial['time'].copy()
    Ylms.month = np.array(12.0*(Ylms.time - 2002.0),dtype='i') + 1
    Ylms.clm = np.zeros((LMAX+1,LMAX+1,nt))
    Ylms.slm = np.zeros((LMAX+1,LMAX+1,nt))
    for t in range(nt):
        #-- convert spatial field to spherical harmonics
        output_Ylms = gen_stokes(input_spatial['data'][t,:,:].T,
            input_spatial['lon'], input_spatial['lat'], UNITS=UNITS,
            LMIN=0, LMAX=LMAX, MMAX=MMAX, PLM=PLM, LOVE=LOVE)
        Ylms.clm[:,:,t] = output_Ylms['clm'][:,:].copy()
        Ylms.slm[:,:,t] = output_Ylms['slm'][:,:].copy()

    #-- if verbose output: print input and output file names
    if VERBOSE:
        print('{0}:'.format(os.path.basename(sys.argv[0])))
        print('{0} -->\n\t{1}'.format(INPUT_FILE,OUTPUT_FILE))
    #-- outputting data to file
    if (DATAFORM == 1):
        #-- ascii (.txt)
        Ylms.to_ascii(OUTPUT_FILE)
    elif (DATAFORM == 2):
        #-- netCDF4 (.nc)
        Ylms.to_netCDF4(OUTPUT_FILE)
    elif (DATAFORM == 3):
        #-- HDF5 (.H5)
        Ylms.to_HDF5(OUTPUT_FILE)
    #-- change output permissions level to MODE
    os.chmod(OUTPUT_FILE,MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' --lmax=X\t\tMaximum spherical harmonic degree')
    print(' --mmax=X\t\tMaximum spherical harmonic order')
    print(' -U X, --units=X\tInput units')
    print('\t1: cm w.e.\n\t2: Gigatonnes (Gt)\n\t3: kg/m^2')
    print(' --reference=X\t\tReference frame for load Love numbers')
    print(' -S X, --spacing=X\tSpatial resolution of input data (dlon,dlat)')
    print(' -I X, --interval=X\tInput grid interval')
    print('\t1: (0:360, 90:-90)\n\t2: (degree spacing/2)')
    print(' --missing\t\tInput spatial fields have missing values')
    print(' --fill-value\t\tSet fill_value for input spatial fields')
    print(' --header\t\tNumber of header rows to skip in input ascii files')
    print(' --delimiter\t\tDelimiter in input ascii files')
    print(' -F X, --format=X\tInput and output data format')
    print('\t1: ascii\n\t2: netcdf\n\t3: HDF5')
    print('-V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermissions mode of the output files\n')

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program and run the analyses
    #-- with the specific parameters.
    short_options = 'hU:S:I:F:VM:'
    long_options = ['help','lmax=','mmax=','units=','reference=','spacing=',
        'interval=','missing','fill-value=','header=','delimiter=','format=',
        'verbose','mode=']
    optlist, arglist = getopt.getopt(sys.argv[1:], short_options, long_options)

    #-- command line parameters
    #-- maximum spherical harmonic degree and order
    LMAX = 60
    MMAX = None
    #-- output units
    UNITS = 1
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    REFERENCE = 'CF'
    #-- input grid parameters
    DDEG = [0.5,0.5]
    INTERVAL = 1
    #-- fill value
    MISSING = False
    FILL_VALUE = None
    #-- ascii parameters
    HEADER = 0
    DELIMITER = None
    #-- input and output data format (1: ascii, 2: netcdf, 3: HDF5)
    DATAFORM = 2
    #-- verbose output
    VERBOSE = False
    #-- permissions mode of the output files (octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-L","--lmax"):
            LMAX = np.int(arg)
        elif opt in ("-U","--units"):
            UNITS = np.int(arg)
        elif opt in ("--reference"):
            REFERENCE = arg.upper()
        elif opt in ("-S","--spacing"):
            DDEG = np.array(arg.split(','),dtype=np.float)
        elif opt in ("-I","--interval"):
            INTERVAL = np.int(arg)
        elif opt in ("-F","--format"):
            DATAFORM = np.int(arg)
        elif opt in ("--missing"):
            MISSING = True
        elif opt in ("--fill-value"):
            FILL_VALUE = np.float(arg)
        elif opt in ("--header"):
            HEADER = np.int(arg)
        elif opt in ("--delimiter"):
            DELIMITER = arg
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- verify inputs
    #-- Input and Output Files
    if not arglist:
        raise Exception('No input and output files specified')

    #-- input spherical harmonic data
    INPUT_FILE = os.path.expanduser(arglist[0])
    #-- output spatial field
    OUTPUT_FILE = os.path.expanduser(arglist[1])
    #-- run program with parameters
    convert_harmonics(INPUT_FILE, OUTPUT_FILE, LMAX=LMAX, MMAX=MMAX,
        UNITS=UNITS, REFERENCE=REFERENCE, DDEG=DDEG, INTERVAL=INTERVAL,
        MISSING=MISSING, FILL_VALUE=FILL_VALUE, HEADER=HEADER,
        DELIMITER=DELIMITER, DATAFORM=DATAFORM, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
