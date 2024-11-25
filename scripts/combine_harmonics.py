#!/usr/bin/env python
u"""
combine_harmonics.py
Written by Tyler Sutterley (10/2023)
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
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use a decorrelation filter (destriping filter)
    -U X, --units X: output units
        0: norm, no unit conversion
        1: cmwe, centimeters water equivalent
        2: mmGH, millimeters geoid height
        3: mmCU, millimeters elastic crustal deformation
        4: micGal, microGal gravity perturbations
        5: mbar, millibars equivalent surface pressure
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
    -D, --date: input and output files have date information
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
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
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
    Updated 10/2023: add date argument to specify if data is a time series
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2023: allow units argument to be set to 0 for no unit conversion
    Updated 03/2023: add index ascii/netCDF4/HDF5 datatypes as possible inputs
        add descriptive file-level attributes to output netCDF4/HDF5 files
        use attributes from units class for writing to netCDF4/HDF5 files
    Updated 02/2023: use get function to retrieve specific units
        use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
        iterate over harmonics objects versus indexing
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 07/2022: create mask for output gridded variables
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
import copy
import logging
import pathlib
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: converts from the spherical harmonic domain into the spatial domain
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
    DATE=False,
    MODE=0o775):

    # verify inputs
    INPUT_FILE = pathlib.Path(INPUT_FILE).expanduser().absolute()
    OUTPUT_FILE = pathlib.Path(OUTPUT_FILE).expanduser().absolute()
    # verify that output directory exists
    OUTPUT_FILE.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
    # attributes for output files
    attributes = dict(ROOT={})
    attributes['ROOT']['product_type'] = 'gravity_field'
    attributes['ROOT']['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # read input spherical harmonic coefficients from file
    if DATAFORM in ('ascii', 'netCDF4', 'HDF5'):
        dataform = copy.copy(DATAFORM)
        input_Ylms = gravtk.harmonics().from_file(INPUT_FILE,
            format=DATAFORM, date=DATE)
        attributes['ROOT']['lineage'] = input_Ylms.filename.name
    elif DATAFORM in ('index-ascii', 'index-netCDF4', 'index-HDF5'):
        # read from index file
        _,dataform = DATAFORM.split('-')
        input_Ylms = gravtk.harmonics().from_index(INPUT_FILE,
            format=dataform, date=DATE)
        attributes['ROOT']['lineage'] = [f.name for f in input_Ylms.filename]
    # reform harmonic dimensions to be l,m,t
    # truncate to degree and order LMAX, MMAX
    input_Ylms = input_Ylms.truncate(lmax=LMAX, mmax=MMAX).expand_dims()

    # remove mean file from input Ylms
    if MEAN_FILE:
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=DATAFORM, date=False)
        input_Ylms.subtract(mean_Ylms)

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth parameters
    attributes['ROOT']['earth_model'] = LOVE.model
    attributes['ROOT']['earth_love_numbers'] = LOVE.citation
    attributes['ROOT']['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['ROOT']['max_degree'] = LMAX
    attributes['ROOT']['max_order'] = MMAX

    # distribute total mass uniformly over the ocean
    if REDISTRIBUTE:
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX,
            MMAX=MMAX, LOVE=LOVE)
        # calculate ratio between total mass and a uniformly distributed
        # layer of water over the ocean
        ratio = input_Ylms.clm[0,0,:]/ocean_Ylms.clm[0,0]
        # for each spherical harmonic
        for m in range(0,MMAX+1):# MMAX+1 to include MMAX
            for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                # remove the ratio*ocean Ylms from Ylms
                # note: x -= y is equivalent to x = x - y
                input_Ylms.clm[l,m,:] -= ratio*ocean_Ylms.clm[l,m]
                input_Ylms.slm[l,m,:] -= ratio*ocean_Ylms.slm[l,m]

    # if using a decorrelation filter (Isabella's destriping Routine)
    if DESTRIPE:
        input_Ylms = input_Ylms.destripe()

    # Gaussian smoothing
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
    else:
        wt = np.ones((LMAX+1))

    # Output spatial data
    grid = gravtk.spatial()
    grid.time = np.copy(input_Ylms.time)
    grid.month = np.copy(input_Ylms.month)
    nt = len(input_Ylms.time)

    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (0:360,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        grid.lon = dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        grid.lon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0, maxlon+dlon/2.0, dlon)
        grid.lat = np.arange(maxlat-dlat/2.0, minlat-dlat/2.0, -dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    # output spatial grid
    grid.data = np.zeros((nlat, nlon, nt))
    grid.mask = np.zeros((nlat, nlon, nt), dtype=bool)

    # Setting units factor for output
    # dfactor computes the degree dependent coefficients
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    # output units and units longname
    # 0: norm, no unit conversion
    # 1: cmwe, centimeters water equivalent
    # 2: mmGH, millimeters geoid height
    # 3: mmCU, millimeters elastic crustal deformation
    # 4: micGal, microGal gravity perturbations
    # 5: mbar, millibars equivalent surface pressure
    units = gravtk.units.bycode(UNITS)
    units_name, units_longname = gravtk.units.get_attributes(units)
    dfactor = factors.get(units)
    # add attributes for earth parameters
    attributes['ROOT']['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['ROOT']['earth_density'] = f'{factors.rho_e:0.3f} g/cm^3'
    attributes['ROOT']['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'

    # Computing plms for converting to spatial domain
    theta = (90.0 - grid.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # converting harmonics to truncated, smoothed coefficients in output units
    for t,Ylms in enumerate(input_Ylms):
        # convolve spherical harmonics with degree dependent factors
        Ylms.convolve(dfactor*wt)
        # convert spherical harmonics to output spatial grid
        grid.data[:,:,t] = gravtk.harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, PLM=PLM).T

    # outputting data to file
    grid.squeeze().to_file(filename=OUTPUT_FILE, format=dataform,
        units=units_name, longname=units_longname,
        attributes=attributes, date=DATE)

    # change output permissions level to MODE
    OUTPUT_FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts a file from the spherical harmonic
            domain into the spatial domain
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # input and output file
    parser.add_argument('infile',
        type=pathlib.Path, nargs='?',
        help='Input harmonic file')
    parser.add_argument('outfile',
        type=pathlib.Path, nargs='?',
        help='Output spatial file')
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
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Verbose output of run')
    # output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[0,1,2,3,4,5],
        help='Output units')
    # output grid parameters
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
    # redistribute total mass over the ocean
    parser.add_argument('--redistribute-mass',
        default=False, action='store_true',
        help='Redistribute total mass over the ocean')
    # land-sea mask for redistributing over the ocean
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing over the ocean')
    # mean file to remove
    parser.add_argument('--mean',
        type=pathlib.Path,
        help='Mean file to remove from the harmonic data')
    # input and output data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=choices,
        help='Input and output data format')
    # Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
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
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program with parameters
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
            DATE=args.date,
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
