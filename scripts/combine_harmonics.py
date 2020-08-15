#!/usr/bin/env python
u"""
combine_harmonics.py
Written by Tyler Sutterley (08/2020)
Converts a file from the spherical harmonic domain into the spatial domain

CALLING SEQUENCE:
    python combine_harmonics.py -F 2 --lmax=60 -U 1 input_file output_file

COMMAND LINE OPTIONS:
    --help: list the command line options
    --lmax=X: maximum spherical harmonic degree
    --mmax=X: maximum spherical harmonic order
    --reference=X: Reference frame for load love numbers
    -R X, --radius=X: Gaussian smoothing radius (km)
    -D, --destripe: use a decorrelation filter (destriping filter)
    -U X, --units=X: output units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: Pa, equivalent surface pressure in Pascals
    -S X, --spacing=X: spatial resolution of output data (dlon,dlat)
    -I X, --interval=X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set bounds with --bounds)
    -B X, --bounds=X: bounding box for interval 3 (minlon,maxlon,minlat,maxlat)
    -O, --ocean: redistribute total mass over the ocean
    --mask=X: input land-sea function (netCDF4) with variable LSMASK as mask
    --mean=X: mean file to remove from the harmonic data
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
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    harmonic_summation.py: calculates a spatial field from spherical harmonics
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
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for input and output operations
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 01/2020: output time in ascii files as the 4th column
    Updated 10/2019: changing Y/N flags to True/False. file can be a time series
        can output a non-global grid by setting bounding box parameters
    Written 07/2018
"""
from __future__ import print_function

import sys
import os
import re
import getopt
import numpy as np

from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.units import units
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
        hl,kl,ll = read_love_numbers(read_love_numbers, FORMAT='tuple',
            REFERENCE=REFERENCE)
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- PURPOSE: converts from the spherical harmonic domain into the spatial domain
def combine_harmonics(INPUT_FILE, OUTPUT_FILE, LMAX=None, MMAX=None,
    REFERENCE=None, RAD=None, DESTRIPE=False, UNITS=None, DDEG=None,
    INTERVAL=None, BOUNDS=None, REDISTRIBUTE=False, LSMASK=None,
    MEAN_FILE=None, DATAFORM=None, VERBOSE=False, MODE=0o775):

    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- read input spherical harmonic coefficients from file in DATAFORM
    if (DATAFORM == 1):
        input_Ylms = harmonics().from_ascii(INPUT_FILE)
    elif (DATAFORM == 2):
        #-- read input netCDF4 file (.nc)
        input_Ylms = harmonics().from_netCDF4(INPUT_FILE)
    elif (DATAFORM == 3):
        #-- read input HDF5 file (.H5)
        input_Ylms = harmonics().from_HDF5(INPUT_FILE)
    #-- reform harmonic dimensions to be l,m,t
    #-- truncate to degree and order LMAX, MMAX
    input_Ylms = input_Ylms.truncate(lmax=LMAX, mmax=MMAX).expand_dims()

    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = load_love_numbers(LMAX, REFERENCE=REFERENCE)

    #-- distribute total mass uniformly over the ocean
    if REDISTRIBUTE:
        #-- read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = ocean_stokes(LSMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))
        #-- calculate ratio between total mass and a uniformly distributed
        #-- layer of water over the ocean
        ratio = input_Ylms.clm[0,0,:]/ocean_Ylms['clm'][0,0]
        #-- for each spherical harmonic
        for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
            for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                #-- remove the ratio*ocean Ylms from Ylms
                #-- note: x -= y is equivalent to x = x - y
                input_Ylms.clm[l,m,:] -= ratio*ocean_Ylms['clm'][l,m]
                input_Ylms.slm[l,m,:] -= ratio*ocean_Ylms['slm'][l,m]

    #-- if using a decorrelation filter (Isabella's destriping Routine)
    if DESTRIPE:
        input_Ylms = input_Ylms.destripe()

    #-- Gaussian smoothing
    if (RAD != 0):
        wt = 2.0*np.pi*gauss_weights(RAD,LMAX)
    else:
        wt = np.ones((LMAX+1))

    #-- Output spatial data
    grid = spatial()
    grid.time = np.copy(input_Ylms.time)
    grid.month = np.copy(input_Ylms.month)

    #-- Output Degree Spacing
    if (len(DDEG) == 1):
        #-- dlon == dlat
        dlon = DDEG
        dlat = DDEG
    else:
        #-- dlon != dlat
        dlon,dlat = DDEG

    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (0:360,90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
        grid.lon = dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        #-- (Degree spacing)/2
        grid.lon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    elif (INTERVAL == 3):
        #-- non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0,maxlon+dlon/2.0,dlon)
        grid.lat = np.arange(maxlat-dlat/2.0,minlat-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)

    #-- Setting units factor for output
    #-- dfactor computes the degree dependent coefficients
    if (UNITS == 1):
        #-- 1: cmwe, centimeters water equivalent
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).cmwe
    elif (UNITS == 2):
        #-- 2: mmGH, mm geoid height
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmGH
    elif (UNITS == 3):
        #-- 3: mmCU, mm elastic crustal deformation
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmCU
    elif (UNITS == 4):
        #-- 4: micGal, microGal gravity perturbations
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).microGal
    elif (UNITS == 5):
        #-- 5: Pa, equivalent surface pressure in Pascals
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).Pa
    else:
        raise ValueError(('UNITS is invalid:\n1: cmwe\n2: mmGH\n3: mmCU '
            '(elastic)\n4:microGal\n5: Pa'))

    #-- Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- output spatial grid
    nt = len(input_Ylms.time)
    grid.data = np.zeros((nlat,nlon,nt))
    #-- converting harmonics to truncated, smoothed coefficients in output units
    for t in range(nt):
        #-- spherical harmonics for time t
        Ylms = input_Ylms.index(t)
        Ylms.convolve(dfactor*wt)
        #-- convert spherical harmonics to output spatial grid
        grid.data[:,:,t] = harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, PLM=PLM).T

    #-- if verbose output: print input and output file names
    if VERBOSE:
        print('{0}:'.format(os.path.basename(sys.argv[0])))
        print('{0} -->\n\t{1}\n'.format(INPUT_FILE,OUTPUT_FILE))
    #-- outputting data to file
    output_data(grid.squeeze(), FILENAME=OUTPUT_FILE,
        DATAFORM=DATAFORM, UNITS=UNITS)
    #-- change output permissions level to MODE
    os.chmod(OUTPUT_FILE,MODE)

#-- PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, DATAFORM=None, UNITS=None):
    #-- output units and units longname
    unit_short = ['cmwe', 'mmGH', 'mmCU', 'microGal', 'Pa']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']
    if (DATAFORM == 1):
        #-- ascii (.txt)
        data.to_ascii(FILENAME)
    elif (DATAFORM == 2):
        #-- netcdf (.nc)
        data.to_netCDF4(FILENAME, units=unit_short[UNITS-1],
            longname=unit_name[UNITS-1])
    elif (DATAFORM == 3):
        #-- HDF5 (.H5)
        data.to_HDF5(FILENAME, units=unit_short[UNITS-1],
            longname=unit_name[UNITS-1])

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' --lmax=X\t\tMaximum spherical harmonic degree')
    print(' --mmax=X\t\tMaximum spherical harmonic order')
    print(' --reference=X\t\tReference frame for load Love numbers')
    print(' -R X, --radius=X\tGaussian smoothing radius (km)')
    print(' -D, --destripe\t\tUse a decorrelation filter')
    print(' -U X, --units=X\tOutput units')
    print('\t1: cm w.e.\n\t2: mm geoid\n\t3: mm uplift\n\t4: microGal\n\t5: Pa')
    print(' -S X, --spacing=X\tSpatial resolution of output data (dlon,dlat)')
    print(' -I X, --interval=X\tOutput grid interval')
    print('\t1: (0:360, 90:-90)\n\t2: (degree spacing/2)\n\t3: non-global grid')
    print(' -B X, --bounds=X\tBounding box for interval 3')
    print(' -O, --ocean\t\tRedistribute total mass over the ocean')
    print(' --mean=X\t\tMean file to remove from the harmonic data')
    print(' -F X, --format=X\tInput and output data format')
    print('\t1: ascii\n\t2: netcdf\n\t3: HDF5')
    print('-V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermissions mode of the output files\n')

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program and run the analyses
    #-- with the specific parameters.
    short_options = 'hR:DU:S:I:B:OF:VM:'
    long_options = ['help','lmax=','mmax=','reference=','radius=','destripe',
        'units=','spacing=','interval=','bounds=','ocean','mean=','format=',
        'verbose','mode=']
    optlist, arglist = getopt.getopt(sys.argv[1:], short_options, long_options)

    #-- command line parameters
    #-- maximum spherical harmonic degree and order
    LMAX = 60
    MMAX = None
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    REFERENCE = 'CF'
    #-- Gaussian smoothing radius (km)
    RAD = 0
    #-- Use a decorrelation (destriping) filter
    DESTRIPE = False
    #-- output units
    UNITS = 1
    #-- output grid parameters
    DDEG = [0.5,0.5]
    INTERVAL = 1
    BOUNDS = None
    #-- redistribute total mass over the ocean
    REDISTRIBUTE = False
    #-- land-sea mask for redistributing over the ocean
    LSMASK = None
    #-- mean file to remove
    MEAN_FILE = None
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
        elif opt in ("--lmax"):
            LMAX = np.int(arg)
        elif opt in ("--mmax"):
            MMAX = np.int(arg)
        elif opt in ("--reference"):
            REFERENCE = arg.upper()
        elif opt in ("-R","--radius"):
            RAD = np.float(arg)
        elif opt in ("-D","--destripe"):
            DESTRIPE = True
        elif opt in ("-U","--units"):
            UNITS = np.int(arg)
        elif opt in ("-S","--spacing"):
            DDEG = np.array(arg.split(','),dtype=np.float)
        elif opt in ("-I","--interval"):
            INTERVAL = np.int(arg)
        elif opt in ("-B","--bounds"):
            BOUNDS = np.array(arg.split(','),dtype=np.float)
        elif opt in ("-O","--ocean"):
            REDISTRIBUTE = True
        elif opt in ("--mask"):
            LSMASK = os.path.expanduser(arg)
        elif opt in ("--mean"):
            MEAN_FILE = os.path.expanduser(arg)
        elif opt in ("-F","--format"):
            DATAFORM = np.int(arg)
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- Input and Output Files
    if not arglist:
        raise Exception('No input and output files specified')

    #-- input spherical harmonic data
    INPUT_FILE = os.path.expanduser(arglist[0])
    #-- output spatial field
    OUTPUT_FILE = os.path.expanduser(arglist[1])
    #-- run program with parameters
    combine_harmonics(INPUT_FILE, OUTPUT_FILE, LMAX=LMAX, MMAX=MMAX,
        RAD=RAD, DESTRIPE=DESTRIPE, UNITS=UNITS, DDEG=DDEG,
        INTERVAL=INTERVAL, BOUNDS=BOUNDS, REDISTRIBUTE=REDISTRIBUTE,
        LSMASK=LSMASK, MEAN_FILE=MEAN_FILE, DATAFORM=DATAFORM,
        VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
