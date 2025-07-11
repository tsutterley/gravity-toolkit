#!/usr/bin/env python
u"""
run_sea_level_equation.py (06/2025)
Solves the sea level equation with the option of including polar motion feedback
Uses a Clenshaw summation to calculate the spherical harmonic summation

CALLING SEQUENCE:
    python run_sea_level_equation.py --lmax 240 --body 0 --fluid 0 \
        --polar-feedback --reference CF --iterations 6 --format netCDF4 \
        --verbose --mode 0o775 input_file output_file

INPUTS:
    input_file: input load file (harmonics or spatial field)
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
    -d X, --density X: Density of water in g/cm^3
    --polar-feedback: Include polar feedback
    --reference X: Reference frame for load love numbers
    -I X, --iterations X: maximum number of iterations for the solver
    -F X, --format X: Input and output data format
        ascii
        netCDF4
        HDF5
    -T X, --input-type X: Input data type for load files
        harmonics: spherical harmonic coefficients
        spatial: spatial fields
    -D, --date: input and output files have date information
    -U X, --units X: input units for spatial files
        1: cm of water thickness (cmwe)
        2: Gigatonnes (Gt)
        3: mm of water thickness kg/m^2
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
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
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
    Updated 06/2025: added options to run from input spatial fields
        added attributes for lineage to track input files
        added option to set the density of water in g/cm^3
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
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
import copy
import logging
import pathlib
import argparse
import traceback
import numpy as np
import collections
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: Computes Sea Level Fingerprints including polar motion feedback
def run_sea_level_equation(INPUT_FILE, OUTPUT_FILE,
    LANDMASK=None,
    LMAX=0,
    LOVE_NUMBERS=0,
    BODY_TIDE_LOVE=0,
    FLUID_LOVE=0,
    DENSITY=1.0,
    REFERENCE=None,
    ITERATIONS=0,
    POLAR=False,
    DATAFORM=None,
    INPUT_TYPE=None,
    DATE=False,
    UNITS=None,
    MODE=0o775):

    # set default paths
    INPUT_FILE = pathlib.Path(INPUT_FILE).expanduser().absolute()
    OUTPUT_FILE = pathlib.Path(OUTPUT_FILE).expanduser().absolute()
    LANDMASK = pathlib.Path(LANDMASK).expanduser().absolute()
    # output attributes for spatial files
    attributes = collections.OrderedDict()
    attributes['lineage'] = INPUT_FILE.name
    attributes['land_sea_mask'] = LANDMASK.name
    attributes['title'] = 'Sea_Level_Fingerprint'
    attributes['max_degree'] = LMAX
    attributes['iterations'] = ITERATIONS

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    # create land function
    nth,nphi = landsea.shape
    land_function = np.zeros((nth, nphi), dtype=np.float64)
    # calculate colatitude in radians
    th = (90.0 - landsea.lat)*np.pi/180.0
    # extract land function from file
    # combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indx,indy] = 1.0

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth model and love numbers
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # add attributes for body tide love numbers
    if (BODY_TIDE_LOVE == 0):
        attributes['earth_body_tide'] = 'Wahr (1981)'
    elif (BODY_TIDE_LOVE == 1):
        attributes['earth_body_tide'] = 'Farrell (1972)'
    # add attributes for fluid love numbers
    if (FLUID_LOVE == 0):
        attributes['earth_fluid_love'] = 'Han and Wahr (1989)'
    elif FLUID_LOVE in (1,2):
        attributes['earth_fluid_love'] = 'Munk and MacDonald (1960)'
    elif (FLUID_LOVE == 3):
        attributes['earth_fluid_love'] = 'Lambeck (1980)'
    # add attribute for true polar wander
    if POLAR:
        attributes['polar_motion_feedback'] = 'Kendall et al. (2005)'

    # read input spherical harmonic coefficients from file
    single_file_formats = ('ascii', 'netCDF4', 'HDF5')
    index_file_formats = ('index-ascii', 'index-netCDF4', 'index-HDF5')
    if DATAFORM in single_file_formats and (INPUT_TYPE == 'spatial'):
        # read spatial data from input file format
        dataform = copy.copy(DATAFORM)
        load_spatial = gravtk.spatial().from_file(INPUT_FILE,
            format=DATAFORM, date=DATE)
        attributes['lineage'] = load_spatial.filename.name
    elif DATAFORM in index_file_formats and (INPUT_TYPE == 'spatial'):
        # read spatial data from index file
        _,dataform = DATAFORM.split('-')
        load_spatial = gravtk.spatial().from_index(INPUT_FILE,
            format=dataform, date=DATE)
        attributes['lineage'] = [f.name for f in load_spatial.filename]
    elif DATAFORM in single_file_formats:
        dataform = copy.copy(DATAFORM)
        # read spherical harmonic coefficients from input file format
        load_Ylms = gravtk.harmonics().from_file(INPUT_FILE,
            format=DATAFORM, date=DATE)
        attributes['lineage'] = load_Ylms.filename.name
    elif DATAFORM in index_file_formats:
        # read spherical harmonic coefficients from index file
        _,dataform = DATAFORM.split('-')
        load_Ylms = gravtk.harmonics().from_index(INPUT_FILE,
            format=dataform, date=DATE)
        attributes['lineage'] = [f.name for f in load_Ylms.filename]
    else:
        raise ValueError(f'Unknown input data format {DATAFORM:s} for {INPUT_TYPE:s}')

    # convert input data to be iterable over time slices
    if (INPUT_TYPE == 'spatial'):
        # expand dimensions to iterate over slices
        load_spatial.expand_dims()
        # number of time slices
        nt = load_spatial.shape[2]
    else:
        # truncate harmonics to degree and order LMAX
        load_Ylms.truncate(lmax=LMAX, mmax=LMAX)
        # expand dimensions to iterate over slices
        load_Ylms.expand_dims()
        # number of time slices
        nt = load_Ylms.shape[2]

    # calculate the legendre functions using Holmes and Featherstone relation
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # allocate for pseudo-spectral sea level equation solver
    sea_level = gravtk.spatial(nlon=nphi, nlat=nth)
    sea_level.data = np.zeros((nth,nphi,nt))
    sea_level.mask = np.zeros((nth,nphi,nt), dtype=bool)
    for i in range(nt):
        # print iteration if running a series
        if (nt > 1):
            logging.info(f'Index {i+1:d} of {nt:d}')
        # subset harmonics/spatial fields to indice
        if (INPUT_TYPE == 'spatial'):
            spatial_data = load_spatial.index(i, date=DATE)
            # convert missing values to zero
            spatial_data.replace_invalid(0.0)
            # convert spatial field to spherical harmonics
            Ylms = gravtk.gen_stokes(spatial_data.data.T,
                spatial_data.lon, spatial_data.lat, UNITS=UNITS,
                LMIN=0, LMAX=LMAX, LOVE=LOVE)
        else:
            Ylms = load_Ylms.index(i, date=DATE)
        # run pseudo-spectral sea level equation solver
        sea_level.data[:,:,i] = gravtk.sea_level_equation(Ylms.clm, Ylms.slm,
            landsea.lon, landsea.lat, land_function.T, LMAX=LMAX,
            LOVE=LOVE, BODY_TIDE_LOVE=BODY_TIDE_LOVE,
            FLUID_LOVE=FLUID_LOVE, DENSITY=DENSITY, POLAR=POLAR,
            PLM=PLM, ITERATIONS=ITERATIONS, FILL_VALUE=0).T
        sea_level.mask[:,:,i] = (sea_level.data[:,:,i] == 0)
    # copy dimensions
    sea_level.lon = np.copy(landsea.lon)
    sea_level.lat = np.copy(landsea.lat)
    # copy date variables
    if (INPUT_TYPE == 'spatial') and DATE:
        # copy time from spatial data
        sea_level.time = np.copy(load_spatial.time)
    elif DATE:
        # copy time from load Ylms
        sea_level.time = np.copy(load_Ylms.time)
    # remove singleton dimensions if necessary
    sea_level.squeeze()
    # add attributes to output spatial field
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    sea_level.attributes['ROOT'] = attributes

    # attributes for output files
    kwargs = {}
    kwargs['units'] = 'centimeters'
    kwargs['longname'] = 'Equivalent_Water_Thickness'
    # save as output DATAFORM
    if (dataform == 'ascii'):
        # ascii (.txt)
        # only print ocean points
        sea_level.fill_value = 0
        sea_level.update_mask()
        sea_level.to_ascii(OUTPUT_FILE, date=DATE)
    elif (dataform == 'netCDF4'):
        # netCDF4 (.nc)
        sea_level.to_netCDF4(OUTPUT_FILE, date=DATE, **kwargs)
    elif (dataform == 'HDF5'):
        # HDF5 (.H5)
        sea_level.to_HDF5(OUTPUT_FILE, date=DATE, **kwargs)
    # set the permissions mode of the output file
    OUTPUT_FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Solves the sea level equation with the option of
            including polar motion feedback
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # input and output file
    parser.add_argument('infile',
        type=pathlib.Path, nargs='?',
        help='Input load file')
    parser.add_argument('outfile',
        type=pathlib.Path, nargs='?',
        help='Output sea level fingerprints file')
    # land mask file
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for calculating sea level fingerprints')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=240,
        help='Maximum spherical harmonic degree')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # different treatments of the body tide Love numbers of degree 2
    # 0: Wahr (1981) and Wahr (1985) values from PREM
    # 1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    parser.add_argument('--body','-b',
        type=int, default=0, choices=[0,1],
        help='Treatment of the body tide Love number')
    # density of water in g/cm^3
    parser.add_argument('--density','-d',
        type=float, default=1.0,
        help='Density of water in g/cm^3')
    # different treatments of the fluid Love number of gravitational potential
    # 0: Han and Wahr (1989) fluid love number
    # 1: Munk and MacDonald (1960) secular love number
    # 2: Munk and MacDonald (1960) fluid love number
    # 3: Lambeck (1980) fluid love number
    parser.add_argument('--fluid','-f',
        type=int, default=0, choices=[0,1,2,3],
        help='Treatment of the fluid Love number')
    # maximum number of iterations for the solver
    # 0th iteration: distribute the water in a uniform layer (barystatic)
    parser.add_argument('--iterations','-I',
        type=int, default=6,
        help='Maximum number of iterations')
    # option for polar feedback
    parser.add_argument('--polar-feedback',
        default=False, action='store_true',
        help='Include effects of polar feedback')
    # option for setting reference frame for load love numbers
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # input and output data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=choices,
        help='Input and output data format')
    # define the input data type for the load files 
    parser.add_argument('--input-type','-T',
        type=str, default='harmonics', choices=['harmonics','spatial'],
        help='Input data type for load fields')
    # Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
    # input units
    # 1: cm of water thickness (cmwe)
    # 2: Gigatonnes (Gt)
    # 3: mm of water thickness kg/m^2
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3],
        help='Input units of spatial fields')
    # print information about processing run
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
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run sea level fingerprints program with parameters
        run_sea_level_equation(args.infile, args.outfile,
            LANDMASK=args.mask,
            LMAX=args.lmax,
            LOVE_NUMBERS=args.love,
            BODY_TIDE_LOVE=args.body,
            FLUID_LOVE=args.fluid,
            DENSITY=args.density,
            REFERENCE=args.reference,
            ITERATIONS=args.iterations,
            POLAR=args.polar_feedback,
            DATAFORM=args.format,
            INPUT_TYPE=args.input_type,
            DATE=args.date,
            UNITS=args.units,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
