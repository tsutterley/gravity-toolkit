#!/usr/bin/env python
u"""
grace_spatial_maps.py
Written by Tyler Sutterley (04/2022)

Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports
    monthly spatial fields

Will correct with the specified GIA model group, destripe/smooth/process,
    and export the data in specified units

Spatial output units: cm w.e., mm geoid height, mm elastic uplift,
    microgal gravity perturbation or surface pressure (mbar)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for spatial files
    -P X, --file-prefix X: prefix string for input and output files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    --lmin X: minimum spherical harmonic degree
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        AW13-IJ05: Geruo A IJ05-R2 GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    --atm-correction: Apply atmospheric jump correction coefficients
    --pole-tide: Correct for pole tide drift
    --geocenter X: Update Degree 1 coefficients with SLR or derived values
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        SLR: satellite laser ranging coefficients from CSR
        SLF: Sutterley and Velicogna coefficients, Remote Sensing (2019)
        Swenson: GRACE-derived coefficients from Sean Swenson
        GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    --interpolate-geocenter: Least-squares model missing Degree 1 coefficients
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
    --slr-c50 X: Replace C50 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    -U X, --units X: output units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: mbar equivalent surface pressure
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
        ascii
        netCDF4
        HDF5
        gfc
    --mask X: Land-sea mask for redistributing land water flux
    --remove-file X: Monthly files to be removed from the GRACE/GRACE-FO data
    --remove-format X: Input data format for files to be removed
        ascii
        netCDF4
        HDF5
        index-ascii
        index-netCDF4
        index-HDF5
    --redistribute-removed: redistribute removed mass fields over the ocean
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces C20,C21,S21,C22,S22,C30 and C50 with SLR values (if specified)
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: converts a land-sea mask to a series of spherical harmonics
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    geocenter.py: converts between spherical harmonics and geocenter variations
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
        option to specify a specific geocenter correction file
        fix default file prefix to include center and release information
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 10/2021: using python logging for handling verbose output
        add more choices for setting input format of the removed files
    Updated 07/2021: simplified file imports using wrappers in harmonics
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 02/2021: changed remove index to files with specified formats
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options and from gfc for mean files
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for output operations
    Updated 05/2020: for public release
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import numpy as np
import argparse
import traceback

import gravity_toolkit.utilities as utilities
from gravity_toolkit.grace_input_months import grace_input_months
from gravity_toolkit.read_GIA_model import read_GIA_model
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

#-- PURPOSE: import GRACE/GRACE-FO files for a given months range
#-- Converts the GRACE/GRACE-FO harmonics applying the specified procedures
def grace_spatial_maps(base_dir, PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    LMIN=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DESTRIPE=False,
    UNITS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    POLE_TIDE=False,
    DEG1=None,
    DEG1_FILE=None,
    MODEL_DEG1=False,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C50=None,
    DATAFORM=None,
    MEAN_FILE=None,
    MEANFORM=None,
    REMOVE_FILES=None,
    REMOVE_FORMAT=None,
    REDISTRIBUTE_REMOVED=False,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    VERBOSE=0,
    MODE=0o775):

    #-- recursively create output directory if not currently existing
    if not os.access(OUTPUT_DIRECTORY, os.F_OK):
        os.makedirs(OUTPUT_DIRECTORY, mode=MODE, exist_ok=True)

    #-- list object of output files for file logs (full path)
    output_files = []

    #-- file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    #-- Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gauss_weights(RAD,LMAX)
        gw_str = '_r{0:0.0f}km'.format(RAD)
    else:
        #-- else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    #-- flag for spherical harmonic order
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''

    #-- reading GRACE months for input date range
    #-- replacing low-degree harmonics with SLR values if specified
    #-- include degree 1 (geocenter) harmonics if specified
    #-- correcting for Pole-Tide and Atmospheric Jumps if specified
    Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX,
        SLR_21=SLR_21, SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C50=SLR_C50,
        DEG1_FILE=DEG1_FILE, MODEL_DEG1=MODEL_DEG1, ATM=ATM,
        POLE_TIDE=POLE_TIDE)
    #-- convert to harmonics object and remove mean if specified
    GRACE_Ylms = harmonics().from_dict(Ylms)
    GRACE_Ylms.directory = Ylms['directory']
    #-- use a mean file for the static field to remove
    if MEAN_FILE:
        #-- read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = harmonics().from_file(MEAN_FILE,format=MEANFORM,date=False)
        #-- remove the input mean
        GRACE_Ylms.subtract(mean_Ylms)
    else:
        GRACE_Ylms.mean(apply=True)
    #-- date information of GRACE/GRACE-FO coefficients
    nfiles = len(GRACE_Ylms.time)

    #-- filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        #-- destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
    else:
        #-- using standard GRACE/GRACE-FO harmonics
        ds_str = ''

    #-- input GIA spherical harmonic datafiles
    GIA_Ylms_rate = read_GIA_model(GIA_FILE,GIA=GIA,LMAX=LMAX,MMAX=MMAX)
    gia_str = '_{0}'.format(GIA_Ylms_rate['title']) if GIA else ''
    #-- calculate the monthly mass change from GIA
    GIA_Ylms = GRACE_Ylms.zeros_like()
    GIA_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    GIA_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    #-- monthly GIA calculated by gia_rate*time elapsed
    #-- finding change in GIA each month
    for t in range(nfiles):
        GIA_Ylms.clm[:,:,t] = GIA_Ylms_rate['clm']*(GIA_Ylms.time[t]-2003.3)
        GIA_Ylms.slm[:,:,t] = GIA_Ylms_rate['slm']*(GIA_Ylms.time[t]-2003.3)

    #-- default file prefix
    if not FILE_PREFIX:
        fargs = (PROC,DREL,DSET,Ylms['title'],gia_str)
        FILE_PREFIX = '{0}_{1}_{2}{3}{4}_'.format(*fargs)

    #-- Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_REMOVED:
        #-- read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = ocean_stokes(LANDMASK,LMAX,MMAX=MMAX,LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        ocean_str = ''

    #-- input spherical harmonic datafiles to be removed from the GRACE data
    #-- Remove sets of Ylms from the GRACE data before returning
    remove_Ylms = GRACE_Ylms.zeros_like()
    remove_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    remove_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    if REMOVE_FILES:
        #-- extend list if a single format was entered for all files
        if len(REMOVE_FORMAT) < len(REMOVE_FILES):
            REMOVE_FORMAT = REMOVE_FORMAT*len(REMOVE_FILES)
        #-- for each file to be removed
        for REMOVE_FILE,REMOVEFORM in zip(REMOVE_FILES,REMOVE_FORMAT):
            if REMOVEFORM in ('ascii','netCDF4','HDF5'):
                #-- ascii (.txt)
                #-- netCDF4 (.nc)
                #-- HDF5 (.H5)
                Ylms = harmonics().from_file(REMOVE_FILE, format=REMOVEFORM)
            elif REMOVEFORM in ('index-ascii','index-netCDF4','index-HDF5'):
                #-- read from index file
                _,removeform = REMOVEFORM.split('-')
                #-- index containing files in data format
                Ylms = harmonics().from_index(REMOVE_FILE, format=removeform)
            #-- reduce to GRACE/GRACE-FO months and truncate to degree and order
            Ylms = Ylms.subset(GRACE_Ylms.month).truncate(lmax=LMAX,mmax=MMAX)
            #-- distribute removed Ylms uniformly over the ocean
            if REDISTRIBUTE_REMOVED:
                #-- calculate ratio between total removed mass and
                #-- a uniformly distributed cm of water over the ocean
                ratio = Ylms.clm[0,0,:]/ocean_Ylms.clm[0,0]
                #-- for each spherical harmonic
                for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
                    for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                        #-- remove the ratio*ocean Ylms from Ylms
                        #-- note: x -= y is equivalent to x = x - y
                        Ylms.clm[l,m,:] -= ratio*ocean_Ylms.clm[l,m]
                        Ylms.slm[l,m,:] -= ratio*ocean_Ylms.slm[l,m]
            #-- filter removed coefficients
            if DESTRIPE:
                Ylms = Ylms.destripe()
            #-- add data for month t and INDEX_FILE to the total
            #-- remove_clm and remove_slm matrices
            #-- redistributing the mass over the ocean if specified
            remove_Ylms.add(Ylms)

    #-- Output spatial data object
    grid = spatial()
    #-- Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        grid.lon = -180 + dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        #-- (Degree spacing)/2
        grid.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
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

    #-- Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- Earth Parameters
    #-- output spatial units
    unit_list = ['cmwe', 'mmGH', 'mmCU', u'\u03BCGal', 'mbar']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']
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
        #-- 5: mbar, millibars equivalent surface pressure
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mbar
    else:
        raise ValueError('Invalid units code {0:d}'.format(UNITS))

    #-- output file format
    file_format = '{0}{1}_L{2:d}{3}{4}{5}_{6:03d}.{7}'
    #-- converting harmonics to truncated, smoothed coefficients in units
    #-- combining harmonics to calculate output spatial fields
    for i,grace_month in enumerate(GRACE_Ylms.month):
        #-- GRACE/GRACE-FO harmonics for time t
        Ylms = GRACE_Ylms.index(i)
        #-- Remove GIA rate for time
        Ylms.subtract(GIA_Ylms.index(i))
        #-- Remove monthly files to be removed
        Ylms.subtract(remove_Ylms.index(i))
        #-- smooth harmonics and convert to output units
        Ylms.convolve(dfactor*wt)
        #-- convert spherical harmonics to output spatial grid
        grid.data = harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMIN=LMIN, LMAX=LMAX,
            MMAX=MMAX, PLM=PLM).T
        #-- copy time variables for month
        grid.time = np.copy(Ylms.time)
        grid.month = np.copy(Ylms.month)

        #-- output monthly files to ascii, netCDF4 or HDF5
        args=(FILE_PREFIX,unit_list[UNITS-1],LMAX,order_str,gw_str,
            ds_str,grace_month,suffix[DATAFORM])
        FILE=os.path.join(OUTPUT_DIRECTORY,file_format.format(*args))
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            grid.to_ascii(FILE, date=True, verbose=VERBOSE)
        elif (DATAFORM == 'netCDF4'):
            #-- netCDF4
            grid.to_netCDF4(FILE, date=True, verbose=VERBOSE,
                units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
                title='GRACE/GRACE-FO Spatial Data')
        elif (DATAFORM == 'HDF5'):
            #-- HDF5
            grid.to_HDF5(FILE, date=True, verbose=VERBOSE,
                units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
                title='GRACE/GRACE-FO Spatial Data')
        #-- set the permissions mode of the output files
        os.chmod(FILE, MODE)
        #-- add file to list
        output_files.append(FILE)

    #-- return the list of output files
    return output_files

#-- PURPOSE: print a file log for the GRACE analysis
def output_log_file(arguments,output_files):
    #-- format: GRACE_processing_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
    #-- create a unique log and open the log file
    DIRECTORY = os.path.expanduser(arguments.output_directory)
    fid = utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    #-- print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(arguments).items()):
        logging.info('{0}: {1}'.format(arg, value))
    #-- print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info('{0}'.format(f))
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE analysis
def output_error_log_file(arguments):
    #-- format: GRACE_processing_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
    #-- create a unique log and open the log file
    DIRECTORY = os.path.expanduser(arguments.output_directory)
    fid = utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    #-- print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(arguments).items()):
        logging.info('{0}: {1}'.format(arg, value))
    #-- print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates monthly spatial maps from GRACE/GRACE-FO
            spherical harmonic coefficients
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
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for spatial files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for input and output files')
    #-- Data processing center or satellite mission
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    #-- GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    #-- minimum spherical harmonic degree
    parser.add_argument('--lmin',
        type=int, default=1,
        help='Minimum spherical harmonic degree')
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
        help='Use decorrelation (destriping) filter')
    #-- output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3,4,5],
        help='Output units')
    #-- output grid parameters
    parser.add_argument('--spacing',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    #-- GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    #-- GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    #-- full path to GIA file
    parser.add_argument('--gia-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GIA file to read')
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
    parser.add_argument('--geocenter-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Specific geocenter file if not default')
    parser.add_argument('--interpolate-geocenter',
        default=False, action='store_true',
        help='Least-squares model missing Degree 1 coefficients')
    #-- replace low degree harmonics with values from Satellite Laser Ranging
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
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    #-- input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    #-- mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    #-- input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    #-- monthly files to be removed from the GRACE/GRACE-FO data
    parser.add_argument('--remove-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Monthly files to be removed from the GRACE/GRACE-FO data')
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--remove-format',
        type=str, nargs='+', choices=choices,
        help='Input data format for files to be removed')
    parser.add_argument('--redistribute-removed',
        default=False, action='store_true',
        help='Redistribute removed mass fields over the ocean')
    #-- land-sea mask for redistributing fluxes
    lsmask = utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing land water flux')
    #-- Output log file for each job in forms
    #-- GRACE_processing_run_2002-04-01_PID-00000.log
    #-- GRACE_processing_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
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
        #-- run grace_spatial_maps algorithm with parameters
        output_files = grace_spatial_maps(
            args.directory,
            args.center,
            args.release,
            args.product,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            LMIN=args.lmin,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DESTRIPE=args.destripe,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            POLE_TIDE=args.pole_tide,
            DEG1=args.geocenter,
            DEG1_FILE=args.geocenter_file,
            MODEL_DEG1=args.interpolate_geocenter,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            REMOVE_FILES=args.remove_file,
            REMOVE_FORMAT=args.remove_format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as e:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        logging.critical('process id {0:d} failed'.format(os.getpid()))
        logging.error(traceback.format_exc())
        if args.log:#-- write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:#-- write successful job completion log file
            output_log_file(args,output_files)

#-- run main program
if __name__ == '__main__':
    main()
