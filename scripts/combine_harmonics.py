#!/usr/bin/env python
u"""
combine_harmonics.py
Written by Tyler Sutterley (04/2022)
Converts a file from the spherical harmonic domain into the spatial domain

CALLING SEQUENCE:
    python combine_harmonics.py -F 2 --lmax 60 -U 1 infile outfile

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
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use a decorrelation filter (destriping filter)
    -U X, --units X: output units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: millibars equivalent surface pressure
    -S X, --spacing X: spatial resolution of output data (dlon,dlat)
    -I X, --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    -B X, --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    --redistribute-mass: redistribute total mass over the ocean
    --mask X: input land-sea function (netCDF4) with variable LSMASK as mask
    --mean X: mean file to remove from the harmonic data
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
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: update grid attributes after allocating for data
    Updated 08/2021: fix spherical harmonic orders if not set
    Updated 07/2021: dded path to default land-sea mask for mass redistribution
    Updated 06/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
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
import logging
import argparse
import traceback
import numpy as np

import gravity_toolkit.utilities as utilities
from gravity_toolkit.read_love_numbers import load_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.units import units

#-- PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: converts from the spherical harmonic domain into the spatial domain
def combine_harmonics(INPUT_FILE, OUTPUT_FILE,
    LMAX=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    RAD=None,
    DESTRIPE=False,
    UNITS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    REDISTRIBUTE=False,
    LANDMASK=None,
    MEAN_FILE=None,
    DATAFORM=None,
    MODE=0o775):

    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- read input spherical harmonic coefficients from file in DATAFORM
    if (DATAFORM == 'ascii'):
        input_Ylms = harmonics().from_ascii(INPUT_FILE)
    elif (DATAFORM == 'netCDF4'):
        #-- read input netCDF4 file (.nc)
        input_Ylms = harmonics().from_netCDF4(INPUT_FILE)
    elif (DATAFORM == 'HDF5'):
        #-- read input HDF5 file (.H5)
        input_Ylms = harmonics().from_HDF5(INPUT_FILE)
    #-- reform harmonic dimensions to be l,m,t
    #-- truncate to degree and order LMAX, MMAX
    input_Ylms = input_Ylms.truncate(lmax=LMAX, mmax=MMAX).expand_dims()
    #-- remove mean file from input Ylms
    if MEAN_FILE and (DATAFORM == 'ascii'):
        mean_Ylms = harmonics().from_ascii(MEAN_FILE,date=False)
        input_Ylms.subtract(mean_Ylms)
    elif MEAN_FILE and (DATAFORM == 'netCDF4'):
        #-- read input netCDF4 file (.nc)
        mean_Ylms = harmonics().from_netCDF4(MEAN_FILE,date=False)
        input_Ylms.subtract(mean_Ylms)
    elif MEAN_FILE and (DATAFORM == 'HDF5'):
        #-- read input HDF5 file (.H5)
        mean_Ylms = harmonics().from_HDF5(MEAN_FILE,date=False)
        input_Ylms.subtract(mean_Ylms)

    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    #-- distribute total mass uniformly over the ocean
    if REDISTRIBUTE:
        #-- read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = ocean_stokes(LANDMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))
        #-- calculate ratio between total mass and a uniformly distributed
        #-- layer of water over the ocean
        ratio = input_Ylms.clm[0,0,:]/ocean_Ylms.clm[0,0]
        #-- for each spherical harmonic
        for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
            for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                #-- remove the ratio*ocean Ylms from Ylms
                #-- note: x -= y is equivalent to x = x - y
                input_Ylms.clm[l,m,:] -= ratio*ocean_Ylms.clm[l,m]
                input_Ylms.slm[l,m,:] -= ratio*ocean_Ylms.slm[l,m]

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
    nt = len(input_Ylms.time)

    #-- Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (0:360,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
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
    #-- output spatial grid
    grid.data = np.zeros((nlat,nlon,nt))
    #-- update attributes
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    #-- Setting units factor for output
    #-- dfactor computes the degree dependent coefficients
    if (UNITS == 1):
        #-- 1: cmwe, centimeters water equivalent
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).cmwe
    elif (UNITS == 2):
        #-- 2: mmGH, millimeters geoid height
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmGH
    elif (UNITS == 3):
        #-- 3: mmCU, millimeters elastic crustal deformation
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmCU
    elif (UNITS == 4):
        #-- 4: micGal, microGal gravity perturbations
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).microGal
    elif (UNITS == 5):
        #-- 5: mbar, millibars equivalent surface pressure
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mbar
    else:
        raise ValueError('Invalid units code {0:d}'.format(UNITS))

    #-- Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- converting harmonics to truncated, smoothed coefficients in output units
    for t in range(nt):
        #-- spherical harmonics for time t
        Ylms = input_Ylms.index(t)
        Ylms.convolve(dfactor*wt)
        #-- convert spherical harmonics to output spatial grid
        grid.data[:,:,t] = harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, PLM=PLM).T

    #-- outputting data to file
    output_data(grid.squeeze(), FILENAME=OUTPUT_FILE,
        DATAFORM=DATAFORM, UNITS=UNITS)
    #-- change output permissions level to MODE
    os.chmod(OUTPUT_FILE,MODE)

#-- PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, DATAFORM=None, UNITS=None):
    #-- output units and units longname
    unit_short = ['cmwe', 'mmGH', 'mmCU', 'microGal', 'mbar']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        data.to_ascii(FILENAME)
    elif (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        data.to_netCDF4(FILENAME, units=unit_short[UNITS-1],
            longname=unit_name[UNITS-1])
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        data.to_HDF5(FILENAME, units=unit_short[UNITS-1],
            longname=unit_name[UNITS-1])

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts a file from the spherical harmonic
            domain into the spatial domain
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    #-- input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input harmonic file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output spatial file')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    #-- Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3,4,5],
        help='Output units')
    #-- output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds','-B',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    #-- redistribute total mass over the ocean
    parser.add_argument('--redistribute-mass',
        default=False, action='store_true',
        help='Redistribute total mass over the ocean')
    #-- land-sea mask for redistributing over the ocean
    lsmask = utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing over the ocean')
    #-- mean file to remove
    parser.add_argument('--mean',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Mean file to remove from the harmonic data')
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    #-- permissions mode of the output files (octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    #-- run program with parameters
    try:
        info(args)
        combine_harmonics(args.infile, args.outfile,
            LMAX=args.lmax,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            RAD=args.radius,
            DESTRIPE=args.destripe,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            REDISTRIBUTE=args.redistribute_mass,
            LANDMASK=args.mask,
            MEAN_FILE=args.mean,
            DATAFORM=args.format,
            MODE=args.mode)
    except Exception as e:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        logging.critical('process id {0:d} failed'.format(os.getpid()))
        logging.error(traceback.format_exc())

#-- run main program
if __name__ == '__main__':
    main()
