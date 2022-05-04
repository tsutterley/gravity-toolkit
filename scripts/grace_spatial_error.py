#!/usr/bin/env python
u"""
grace_spatial_error.py
Written by Tyler Sutterley (04/2022)

Calculates the GRACE/GRACE-FO errors following Wahr et al. (2006)

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
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
        harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces C20,C21,S21,C22,S22,C30 and C50 with SLR values (if specified)
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units
    tssmooth.py: smoothes a time-series for seasonal effects
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

REFERENCES:
    J Wahr, S C Swenson, I Velicogna, "Accuracy of GRACE mass estimates",
        Geophysical Research Letters, 33(6), L06401, (2006)
        http://dx.doi.org/10.1029/2005GL025305

UPDATE HISTORY:
    Updated 04/2022: use wrapper function for reading load Love numbers
     use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
        option to specify a specific geocenter correction file
        fix default file prefix to include center and release information
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: simplified file imports using wrappers in harmonics
        remove choices for argparse processing centers
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 12/2020: added more love number options and from gfc for mean files
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for output operations
    Updated 04/2020: updates to reading load love numbers
        using the units class for converting normalized spherical harmonics
    Updated 03/2018: added option for output in equivalent pressure (UNITS=5)
        simplified love number extrapolation if LMAX is greater than 696
    Updated 06/2016: using __future__ print function, updated output file name
        updated UNITS comments similar to GRACE output programs
        added parameter MMAX for LMAX != MMAX
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 06/2015: added output_files for log files
    Updated 11/2014: added HDF5 output, added option interval
    Updated 06/2014: changed message to sys.exit
    Updated 09/2013: calculate error using the RMS of the residual
        calculated from the 13-month smoothing
    Updated 05/2013: algorithm updates following python processing scheme
    Written 08/2012
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
from gravity_toolkit.read_love_numbers import load_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.tssmooth import tssmooth
from gravity_toolkit.units import units

#-- PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: import GRACE files for a given months range
#-- Estimates the GRACE/GRACE-FO errors applying the specified procedures
def grace_spatial_error(base_dir, PROC, DREL, DSET, LMAX, RAD,
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
    #-- atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''

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
    #-- full path to directory for specific GRACE/GRACE-FO product
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

    #-- default file prefix
    if not FILE_PREFIX:
        FILE_PREFIX = '{0}_{1}_{2}{3}_'.format(PROC,DREL,DSET,Ylms['title'])

    #-- filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        #-- destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
    else:
        #-- using standard GRACE/GRACE-FO harmonics
        ds_str = ''

    #-- calculating GRACE error (Wahr et al 2006)
    #-- output GRACE error file (for both LMAX==MMAX and LMAX != MMAX cases)
    args = (PROC,DREL,DSET,LMAX,order_str,ds_str,atm_str,GRACE_Ylms.month[0],
        GRACE_Ylms.month[-1],suffix[DATAFORM])
    delta_format = '{0}_{1}_{2}_DELTA_CLM_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    DELTA_FILE = os.path.join(GRACE_Ylms.directory,delta_format.format(*args))
    #-- full path of the GRACE directory
    #-- if file was previously calculated, will read file
    #-- else will calculate the GRACE error
    if not os.access(DELTA_FILE, os.F_OK):
        #-- add output delta file to list object
        output_files.append(DELTA_FILE)

        #-- Delta coefficients of GRACE time series (Error components)
        delta_Ylms = harmonics(lmax=LMAX,mmax=MMAX)
        delta_Ylms.clm = np.zeros((LMAX+1,MMAX+1))
        delta_Ylms.slm = np.zeros((LMAX+1,MMAX+1))
        #-- Smoothing Half-Width (CNES is a 10-day solution)
        #-- 365/10/2 = 18.25 (next highest is 19)
        #-- All other solutions are monthly solutions (HFWTH for annual = 6)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            HFWTH = 19
        else:
            HFWTH = 6
        #-- Equal to the noise of the smoothed time-series
        #-- for each spherical harmonic order
        for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
            #-- for each spherical harmonic degree
            for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                #-- Delta coefficients of GRACE time series
                for cs,csharm in enumerate(['clm','slm']):
                    #-- Constrained GRACE Error (Noise of smoothed time-series)
                    #-- With Annual and Semi-Annual Terms
                    val1 = getattr(GRACE_Ylms, csharm)
                    smth = tssmooth(GRACE_Ylms.time, val1[l,m,:], HFWTH=HFWTH)
                    #-- number of smoothed points
                    nsmth = len(smth['data'])
                    #-- GRACE delta Ylms
                    #-- variance of data-(smoothed+annual+semi)
                    val2 = getattr(delta_Ylms, csharm)
                    val2[l,m] = np.sqrt(np.sum(smth['noise']**2)/nsmth)

        #-- save GRACE DELTA to file
        delta_Ylms.time = np.copy(nsmth)
        delta_Ylms.month = np.copy(nsmth)
        delta_Ylms.to_file(DELTA_FILE,format=DATAFORM)
        #-- set the permissions mode of the output harmonics file
        os.chmod(DELTA_FILE, MODE)
        #-- append delta harmonics file to output files list
        output_files.append(DELTA_FILE)
    else:
        #-- read GRACE DELTA spherical harmonics datafile
        delta_Ylms = harmonics().from_file(DELTA_FILE,format=DATAFORM)
        #-- truncate grace delta clm and slm to d/o LMAX/MMAX
        delta_Ylms = delta_Ylms.truncate(lmax=LMAX, mmax=MMAX)
        nsmth = np.int64(delta_Ylms.time)

    #-- Output spatial data object
    delta = spatial()
    #-- Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        delta.lon = -180 + dlon*np.arange(0,nlon)
        delta.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        #-- (Degree spacing)/2
        delta.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
        delta.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(delta.lon)
        nlat = len(delta.lat)
    elif (INTERVAL == 3):
        #-- non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        delta.lon = np.arange(minlon+dlon/2.0,maxlon+dlon/2.0,dlon)
        delta.lat = np.arange(maxlat-dlat/2.0,minlat-dlat/2.0,-dlat)
        nlon = len(delta.lon)
        nlat = len(delta.lat)

    #-- Earth Parameters
    #-- output spatial units
    unit_list = ['cmwe', 'mmGH', 'mmCU', u'\u03BCGal', 'mbar']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']
    #-- dfactor is the degree dependent coefficients
    #-- for specific spherical harmonic output units
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
        #-- 5: mbar, millibar equivalent surface pressure
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mbar

    #-- Computing plms for converting to spatial domain
    phi = delta.lon[np.newaxis,:]*np.pi/180.0
    theta = (90.0-delta.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))
    #-- square of legendre polynomials truncated to order MMAX
    mm = np.arange(0,MMAX+1)
    PLM2 = PLM[:,mm,:]**2

    #-- Calculating cos(m*phi)^2 and sin(m*phi)^2
    m = delta_Ylms.m[:,np.newaxis]
    ccos = np.cos(np.dot(m,phi))**2
    ssin = np.sin(np.dot(m,phi))**2

    #-- truncate delta harmonics to spherical harmonic range
    Ylms = delta_Ylms.truncate(LMAX,lmin=LMIN,mmax=MMAX)
    #-- convolve delta harmonics with degree dependent factors
    #-- smooth harmonics and convert to output units
    Ylms = Ylms.convolve(dfactor*wt).power(2.0).scale(1.0/nsmth)
    #-- Calculate fourier coefficients
    d_cos = np.zeros((MMAX+1,nlat))#-- [m,th]
    d_sin = np.zeros((MMAX+1,nlat))#-- [m,th]
    #-- Calculating delta spatial values
    for k in range(0,nlat):
        #-- summation over all spherical harmonic degrees
        d_cos[:,k] = np.sum(PLM2[:,:,k]*Ylms.clm, axis=0)
        d_sin[:,k] = np.sum(PLM2[:,:,k]*Ylms.slm, axis=0)

    #-- Multiplying by c/s(phi#m) to get spatial maps (lon,lat)
    delta.data=np.sqrt(np.dot(ccos.T,d_cos) + np.dot(ssin.T,d_sin)).T

    #-- output file format
    file_format = '{0}{1}_L{2:d}{3}{4}{5}_ERR_{6:03d}-{7:03d}.{8}'
    #-- output error file to ascii, netCDF4 or HDF5
    args = (FILE_PREFIX,unit_list[UNITS-1],LMAX,order_str,gw_str,ds_str,
        GRACE_Ylms.month[0],GRACE_Ylms.month[-1],suffix[DATAFORM])
    FILE = os.path.join(OUTPUT_DIRECTORY,file_format.format(*args))
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        delta.to_ascii(FILE, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4
        delta.to_netCDF4(FILE, date=False, verbose=VERBOSE,
            units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
            title='GRACE/GRACE-FO Spatial Error')
    elif (DATAFORM == 'HDF5'):
        #-- HDF5
        delta.to_HDF5(FILE, date=False, verbose=VERBOSE,
            units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
            title='GRACE/GRACE-FO Spatial Error')
    #-- set the permissions mode of the output files
    os.chmod(FILE, MODE)
    #-- add file to list
    output_files.append(FILE)

    #-- return the list of output files
    return output_files

#-- PURPOSE: print a file log for the GRACE analysis
def output_log_file(arguments,output_files):
    #-- format: GRACE_error_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_error_run_{0}_PID-{1:d}.log'.format(*args)
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
    #-- format: GRACE_error_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_error_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Calculates the GRACE/GRACE-FO spatial errors
            following Wahr et al. (2006)
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
    #-- Output log file for each job in forms
    #-- GRACE_error_run_2002-04-01_PID-00000.log
    #-- GRACE_error_failed_run_2002-04-01_PID-00000.log
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
        #-- run grace_spatial_error algorithm with parameters
        output_files = grace_spatial_error(
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
