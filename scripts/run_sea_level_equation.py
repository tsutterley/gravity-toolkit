#!/usr/bin/env python
u"""
run_sea_level_equation.py (04/2022)
Solves the sea level equation with the option of including polar motion feedback
Uses a Clenshaw summation to calculate the spherical harmonic summation

CALLING SEQUENCE:
    python run_sea_level_equation.py --lmax 240 --body 0 --fluid 0 \
        --polar-feedback --reference CF --iterations 6 --format netCDF4 \
        --verbose --mode 0o775 input_file output_file

INPUTS:
    input_file: input load file
    output_file: output sea level fingerprints file

COMMAND LINE OPTIONS:
    --mask X: input land-sea function (netCDF4) with variable LSMASK as mask
    -l X, --lmax X: Maximum spherical harmonic degree
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -b X, --body X: Treatment of the body tide Love number
        0: Wahr (1981) and Wahr (1985) values from PREM
        1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    -f X, --fluid X: Treatment of the fluid Love number
        0: Han and Wahr (1989) fluid love number
        1: Munk and MacDonald (1960) secular love number
        2: Munk and MacDonald (1960) fluid love number
        3: Lambeck (1980) fluid love number
    --polar-feedback: Include polar feedback
    --reference X: Reference frame for load love numbers
    -I X, --iterations X: maximum number of iterations for the solver
    -F X, --format X: Input and output data format
        ascii
        netCDF4
        HDF5
    -D, --date: input and output files have date information
    -V, --verbose: verbose output of processing run
    -M X, --mode X: permissions mode of the output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    sea_level_equation.py: pseudo-spectral sea level equation solver
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

REFERENCES:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        http://dx.doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
        added path to default land-sea mask for sea level equation
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: can run for a spherical harmonic time series
    Updated 12/2020: added more load love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for input and output operations
    Updated 04/2020: updates to reading load love numbers
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: parameters for Munk and MacDonald (1960) fluid love number
    Updated 06/2018: using python3 compatible octal and input
    Updated 04/2018: added option to set the gravitational load Love number
    Updated 03/2018: updated header text.  --format option sets data format
        added extrapolation of load love numbers for LMAX > 696
    Updated 02/2018: spherical harmonic order same as spherical harmonic degree
    Updated 09/2017: ocean mask with GEUS Greenland coastlines.  added comments
        more calculations shifted from spatial domain to spherical harmonics
    Updated 08/2017: different landsea mask (additional changes in West Antarctica)
        Clenshaw summation of spherical harmonics.  Running at 0.5x0.5 degrees
        iterate until convergence of spherical harmonic modulus residuals
    Updated 07/2017: outputs of legendre_polynomials.py include derivatives now
    Updated 04/2017: finished writing algorithms. should be in working condition
        set the permissions mode of the output files with --mode
    Updated 05/2017: slight edit to input distribution 3.  added some comments
    Updated 04/2017: set the permissions mode of the output files with --mode
    Written 09/2016
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
from gravity_toolkit.sea_level_equation import sea_level_equation
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial

#-- PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: Computes Sea Level Fingerprints including polar motion feedback
def run_sea_level_equation(INPUT_FILE, OUTPUT_FILE,
    LANDMASK=None,
    LMAX=0,
    LOVE_NUMBERS=0,
    BODY_TIDE_LOVE=0,
    FLUID_LOVE=0,
    REFERENCE=None,
    ITERATIONS=0,
    POLAR=False,
    DATAFORM=None,
    DATE=False,
    MODE=0o775):

    #-- Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    #-- 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    #-- Open the land-sea NetCDF file for reading
    landsea = spatial().from_netCDF4(LANDMASK, date=False, varname='LSMASK')
    #-- create land function
    nth,nphi = landsea.shape
    land_function = np.zeros((nth,nphi),dtype=np.float64)
    #-- calculate colatitude in radians
    th = (90.0 - landsea.lat)*np.pi/180.0
    #-- extract land function from file
    #-- combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indx,indy] = 1.0

    #-- read load love numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    #-- read spherical harmonic coefficients of the load from input DATAFORM
    if (DATAFORM == 'ascii'):
        #-- read input ascii file (.txt)
        load_Ylms = harmonics().from_ascii(INPUT_FILE, date=DATE)
    elif (DATAFORM == 'netCDF4'):
        #-- read input netCDF4 file (.nc)
        load_Ylms = harmonics().from_netCDF4(INPUT_FILE, date=DATE)
    elif (DATAFORM == 'HDF5'):
        #-- read input HDF5 file (.H5)
        load_Ylms = harmonics().from_HDF5(INPUT_FILE, date=DATE)
    #-- truncate harmonics to degree and order LMAX
    load_Ylms.truncate(lmax=LMAX, mmax=LMAX)
    #-- expand dimensions to iterate over slices
    load_Ylms.expand_dims()
    l1,m1,nt = load_Ylms.shape

    #-- calculate the legendre functions using Holmes and Featherstone relation
    PLM,dPLM = plm_holmes(LMAX, np.cos(th))

    #-- allocate for pseudo-spectral sea level equation solver
    sea_level = spatial(nlon=nphi, nlat=nth)
    sea_level.data = np.zeros((nth,nphi,nt))
    sea_level.mask = np.zeros((nth,nphi,nt), dtype=bool)
    for i in range(nt):
        #-- print iteration if running a series
        if (nt > 1):
            logging.info('Index {0:d} of {1:d}'.format(i+1,nt))
        #-- subset harmonics to indice
        Ylms = load_Ylms.index(i, date=DATE)
        #-- run pseudo-spectral sea level equation solver
        sea_level.data[:,:,i] = sea_level_equation(Ylms.clm, Ylms.slm,
            landsea.lon, landsea.lat, land_function.T, LMAX=LMAX,
            LOVE=(hl,kl,ll), BODY_TIDE_LOVE=BODY_TIDE_LOVE,
            FLUID_LOVE=FLUID_LOVE, POLAR=POLAR, PLM=PLM,
            ITERATIONS=ITERATIONS, FILL_VALUE=0).T
        sea_level.mask[:,:,i] = (sea_level.data[:,:,i] == 0)
    #-- copy dimensions
    sea_level.lon = np.copy(landsea.lon)
    sea_level.lat = np.copy(landsea.lat)
    sea_level.time = np.copy(load_Ylms.time) if DATE else None
    #-- remove singleton dimensions if necessary
    sea_level.squeeze()

    #-- save as output DATAFORM
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        #-- only print ocean points
        sea_level.fill_value = 0
        sea_level.update_mask()
        sea_level.to_ascii(OUTPUT_FILE, date=DATE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4 (.nc)
        sea_level.to_netCDF4(OUTPUT_FILE, date=DATE, units='centimeters',
            longname='Equivalent_Water_Thickness', title='Sea_Level_Fingerprint')
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        sea_level.to_HDF5(OUTPUT_FILE, date=DATE, units='centimeters',
            longname='Equivalent_Water_Thickness', title='Sea_Level_Fingerprint')
    #-- set the permissions mode of the output file
    os.chmod(OUTPUT_FILE, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Solves the sea level equation with the option of
            including polar motion feedback
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    #-- input and output file
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Input load file')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='Output sea level fingerprints file')
    #-- land mask file
    lsmask = utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for calculating sea level fingerprints')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=240,
        help='Maximum spherical harmonic degree')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- different treatments of the body tide Love numbers of degree 2
    #-- 0: Wahr (1981) and Wahr (1985) values from PREM
    #-- 1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    parser.add_argument('--body','-b',
        type=int, default=0, choices=[0,1],
        help='Treatment of the body tide Love number')
    #-- different treatments of the fluid Love number of gravitational potential
    #-- 0: Han and Wahr (1989) fluid love number
    #-- 1: Munk and MacDonald (1960) secular love number
    #-- 2: Munk and MacDonald (1960) fluid love number
    #-- 3: Lambeck (1980) fluid love number
    parser.add_argument('--fluid','-f',
        type=int, default=0, choices=[0,1,2,3],
        help='Treatment of the fluid Love number')
    #-- maximum number of iterations for the solver
    #-- 0th iteration: distribute the water in a uniform layer (barystatic)
    parser.add_argument('--iterations','-I',
        type=int, default=6,
        help='Maximum number of iterations')
    #-- option for polar feedback
    parser.add_argument('--polar-feedback',
        default=False, action='store_true',
        help='Include effects of polar feedback')
    #-- option for setting reference frame for load love numbers
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
    #-- print information about processing run
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

    #-- try to run the analysis with listed parameters
    try:
        info(args)
        #-- run sea level fingerprints program with parameters
        run_sea_level_equation(args.infile, args.outfile,
            LANDMASK=args.mask,
            LMAX=args.lmax,
            LOVE_NUMBERS=args.love,
            BODY_TIDE_LOVE=args.body,
            FLUID_LOVE=args.fluid,
            REFERENCE=args.reference,
            ITERATIONS=args.iterations,
            POLAR=args.polar_feedback,
            DATAFORM=args.format,
            DATE=args.date,
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
