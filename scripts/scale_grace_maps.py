#!/usr/bin/env python
u"""
scale_grace_maps.py
Written by Tyler Sutterley (04/2021)

Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports
    monthly scaled spatial fields, estimated scaling errors,
    and estimated scaled delta errors

Will correct with the specified GIA model group, destripe/smooth/process,
    and export the data in centimeters water equivalent (cm w.e.)

SYSTEM ARGUMENTS README:
    program is run as:
    python scale_grace_maps.py inp1 inp2 inp3
        where inp1, inp2 and inp3 are different inputs

        firstinput=sys.argv[1] (in this case inp1)
        secondinput=sys.argv[2] (in this case inp2)
        thirdinput=sys.argv[3] (in this case inp3)

    As python is base 0, sys.argv[0] is equal to scale_grace_maps.py
        (which is useful in some applications, but not for this program)

    For this program, the system arguments are parameter files
    The program reads the parameter file, which is separated by column as:
        Column 1: parameter name (such as LMAX)
        Column 2: parameter (e.g. 60)
        Column 3: comments (which are discarded)
    The parameters are stored in a python dictionary (variables indexed by keys)
        the keys are the parameter name (for LMAX: parameters['LMAX'] == 60)

INPUTS:
    parameter files containing specific variables for each analysis

COMMAND LINE OPTIONS:
    --help: list the command line options
    -P X, --np X: Run in parallel with X number of processes
    -D X, --directory X: Working data directory
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -r X, --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -l, --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
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
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified date range
        Includes degree 1 values (if specified)
        Replaces C20,C21,S21,C22,S22,C30 and C50 with SLR values (if specified)
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    geocenter.py: converts between spherical harmonics and geocenter variations
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    tssmooth.py: smoothes a time-series using a 13-month Loess-type algorithm
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
    time.py: utilities for calculating time operations
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    C-W Hsu and I Velicogna, "Detection of Sea Level Fingerprints derived from
        GRACE gravity data", Geophysical Research Letters, 44(17), (2017).
        https://doi.org/10.1002/2017GL074070

    F W Landerer and S C Swenson, "Accuracy of scaled GRACE terrestrial water
        storage estimates", Water Resources Research, 48(4), W04531, (2012).
        https://doi.org/10.1029/2011WR011453

    J Wahr, S C Swenson, and I Velicogna, "Accuracy of GRACE mass estimates",
        Geophysical Research Letters, 33(6), L06401, (2006).
        https://doi.org/10.1029/2005GL025305

UPDATE HISTORY:
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 02/2021: changed remove index to files with specified formats
    Updated 02/2021: for public release
"""
from __future__ import print_function, division

import sys
import os
import time
import argparse
import numpy as np
import multiprocessing
import traceback

from gravity_toolkit.grace_input_months import grace_input_months
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.units import units
from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.tssmooth import tssmooth
from gravity_toolkit.utilities import get_data_path

#-- PURPOSE: keep track of multiprocessing threads
def info(title):
    print(os.path.basename(sys.argv[0]))
    print(title)
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: read load love numbers for the range of spherical harmonic degrees
def load_love_numbers(LMAX, LOVE_NUMBERS=0, REFERENCE='CF'):
    """
    Reads PREM load Love numbers for the range of spherical harmonic degrees
    and applies isomorphic parameters

    Arguments
    ---------
    LMAX: maximum spherical harmonic degree

    Keyword arguments
    -----------------
    LOVE_NUMBERS: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth

    Returns
    -------
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement
    """
    #-- load love numbers file
    if (LOVE_NUMBERS == 0):
        #-- PREM outputs from Han and Wahr (1995)
        #-- https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
        love_numbers_file = get_data_path(['data','love_numbers'])
        header = 2
        columns = ['l','hl','kl','ll']
    elif (LOVE_NUMBERS == 1):
        #-- PREM outputs from Gegout (2005)
        #-- http://gemini.gsfc.nasa.gov/aplo/
        love_numbers_file = get_data_path(['data','Load_Love2_CE.dat'])
        header = 3
        columns = ['l','hl','ll','kl']
    elif (LOVE_NUMBERS == 2):
        #-- PREM outputs from Wang et al. (2012)
        #-- https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREM-LLNs-truncated.dat'])
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, LMAX=LMAX, HEADER=header,
        COLUMNS=columns, REFERENCE=REFERENCE, FORMAT='tuple')
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- PURPOSE: import GRACE/GRACE-FO files for a given months range
#-- Calculates monthly scaled spatial maps from GRACE/GRACE-FO
#-- spherical harmonic coefficients
def scale_grace_maps(base_dir, parameters, LOVE_NUMBERS=0, REFERENCE=None,
    VERBOSE=False, MODE=0o775):
    #-- convert parameters to variables
    #-- Data processing center
    PROC = parameters['PROC']
    #-- Data Release
    DREL = parameters['DREL']
    #-- GRACE dataset
    DSET = parameters['DSET']
    #-- Date Range and missing months
    start_mon = np.int(parameters['START'])
    end_mon = np.int(parameters['END'])
    missing = np.array(parameters['MISSING'].split(','),dtype=np.int)
    #-- spherical harmonic degree range
    LMIN = np.int(parameters['LMIN'])
    LMAX = np.int(parameters['LMAX'])
    #-- maximum spherical harmonic order
    if (parameters['MMAX'].title() == 'None'):
        MMAX = np.copy(LMAX)
    else:
        MMAX = np.int(parameters['MMAX'])
    #-- degree 1 coefficients
    #-- None: No degree 1
    #-- Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
    #--     https://grace.jpl.nasa.gov/data/get-data/geocenter/
    #-- SLR: satellite laser ranging from CSR
    #--     ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
    #-- SLF: Sutterley and Velicogna, Remote Sensing (2019)
    #--     https://www.mdpi.com/2072-4292/11/18/2108
    DEG1 = parameters['DEG1']
    MODEL_DEG1 = parameters['MODEL_DEG1'] in ('Y','y')
    #-- replace low-degree coefficients with values from SLR
    SLR_C20 = parameters['SLR_C20']
    SLR_21 = parameters['SLR_21']
    SLR_22 = parameters['SLR_22']
    SLR_C30 = parameters['SLR_C30']
    SLR_C50 = parameters['SLR_C50']
    #-- ECMWF jump corrections
    ATM = parameters['ATM'] in ('Y','y')
    #-- Pole-Tide from Wahr et al. (2015)
    POLE_TIDE = parameters['POLE_TIDE'] in ('Y','y')
    #-- Glacial Isostatic Adjustment file to read
    GIA = parameters['GIA'] if (parameters['GIA'].title() != 'None') else None
    GIA_FILE = os.path.expanduser(parameters['GIA_FILE'])
    #-- remove a set of spherical harmonics from the GRACE data
    REDISTRIBUTE_REMOVED = parameters['REDISTRIBUTE_REMOVED'] in ('Y','y')
    #-- gaussian smoothing radius
    RAD = np.float(parameters['RAD'])
    #-- filter coefficients for stripe effects
    DESTRIPE = parameters['DESTRIPE'] in ('Y','y')
    #-- output degree spacing
    #-- can enter dlon and dlat as [dlon,dlat] or a single value
    DDEG = np.squeeze(np.array(parameters['DDEG'].split(','),dtype='f'))
    #-- output degree interval (0:360, 90:-90) or (degree spacing/2)
    INTERVAL = np.int(parameters['INTERVAL'])
    #-- input/output data format (ascii, netCDF4, HDF5)
    DATAFORM = parameters['DATAFORM']
    #-- output directory and base filename
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    FILENAME = parameters['FILENAME']

    #-- recursively create output Directory if not currently existing
    if (not os.access(DIRECTORY, os.F_OK)):
        os.makedirs(DIRECTORY, mode=MODE, exist_ok=True)

    #-- list object of output files for file logs (full path)
    output_files = []

    #-- file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    #-- output file format
    file_format = '{0}{1}{2}_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'

    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    #-- atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- output spatial units
    unit_str = 'cmwe'
    unit_name = 'Equivalent Water Thickness'
    #-- invalid value
    fill_value = -9999.0

    #-- Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gauss_weights(RAD,LMAX)
        gw_str = '_r{0:0.0f}km'.format(RAD)
    else:
        #-- else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    #-- Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_REMOVED:
        #-- read Land-Sea Mask and convert to spherical harmonics
        LSMASK = os.path.expanduser(parameters['LANDMASK'])
        ocean_Ylms = ocean_stokes(LSMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))

    #-- Grid spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    #-- Grid dimensions
    if (INTERVAL == 1):#-- (0:360, 90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
    elif (INTERVAL == 2):#-- degree spacing/2
        nlon = np.int((360.0/dlon))
        nlat = np.int((180.0/dlat))

    #-- read data for input scale files (ascii, netCDF4, HDF5)
    if (DATAFORM == 'ascii'):
        kfactor = spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon).from_ascii(
            parameters['SCALE_FILE'],date=False)
        k_error = spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon).from_ascii(
            parameters['ERROR_FILE'],date=False)
        k_power = spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon).from_ascii(
            parameters['POWER_FILE'],date=False)
    elif (DATAFORM == 'netCDF4'):
        kfactor = spatial().from_netCDF4(parameters['SCALE_FILE'],date=False)
        k_error = spatial().from_netCDF4(parameters['ERROR_FILE'],date=False)
        k_power = spatial().from_netCDF4(parameters['POWER_FILE'],date=False)
    elif (DATAFORM == 'HDF5'):
        kfactor = spatial().from_HDF5(parameters['SCALE_FILE'],date=False)
        k_error = spatial().from_HDF5(parameters['ERROR_FILE'],date=False)
        k_power = spatial().from_HDF5(parameters['POWER_FILE'],date=False)
    #-- input data shape
    nlat,nlon = kfactor.shape

    #-- input GRACE/GRACE-FO spherical harmonic datafiles
    #-- reading GRACE months for input range with grace_input_months.py
    #-- replacing low-degree harmonics with SLR values if specified
    #-- include degree 1 (geocenter) harmonics if specified
    #-- correcting for Pole-Tide and Atmospheric Jumps if specified
    Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        start_mon, end_mon, missing, SLR_C20, DEG1, MMAX=MMAX,
        SLR_21=SLR_21, SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C50=SLR_C50,
        MODEL_DEG1=MODEL_DEG1, ATM=ATM, POLE_TIDE=POLE_TIDE)
    #-- full path to directory for specific GRACE/GRACE-FO product
    grace_dir = Ylms['directory']
    #-- create harmonics object from GRACE/GRACE-FO data
    GRACE_Ylms = harmonics().from_dict(Ylms)
    GRACE_Ylms.directory = Ylms['directory']
    GRACE_Ylms.title = Ylms['title']
    #-- mean parameters
    MEAN = parameters['MEAN'] in ('Y','y')
    #-- use a mean file for the static field to remove
    if (parameters['MEAN_FILE'].title() == 'None'):
        mean_Ylms = GRACE_Ylms.mean(apply=MEAN)
    else:
        #-- read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        if (parameters['MEANFORM'] == 'ascii'):
            mean_Ylms=harmonics().from_ascii(parameters['MEAN_FILE'],date=False)
        elif (parameters['MEANFORM'] == 'netCDF4'):
            mean_Ylms=harmonics().from_netCDF4(parameters['MEAN_FILE'],date=False)
        elif (parameters['MEANFORM'] == 'HDF5'):
            mean_Ylms=harmonics().from_HDF5(parameters['MEAN_FILE'],date=False)
        elif (parameters['MEANFORM'] == 'gfc'):
            mean_Ylms=harmonics().from_gfc(parameters['MEAN_FILE'])
        #-- remove the input mean
        if MEAN:
            GRACE_Ylms.subtract(mean_Ylms)
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

    #-- input spherical harmonic datafiles to be removed from the GRACE data
    #-- Remove sets of Ylms from the GRACE data before returning
    remove_Ylms = GRACE_Ylms.zeros_like()
    remove_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    remove_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    if (parameters['REMOVE_FILE'].title() != 'None'):
        #-- files to be removed and their respective formats
        REMOVE_FILES = parameters['REMOVE_FILE'].split(',')
        FORMATS = parameters['REMOVEFORM'].split(',')
        #-- extend list if a single format was entered for all files
        if len(FORMATS) < len(REMOVE_FILES):
            FORMATS = FORMATS*len(REMOVE_FILES)
        #-- for each file to be removed
        for REMOVE_FILE,REMOVEFORM in zip(REMOVE_FILES,FORMATS):
            if (REMOVEFORM == 'ascii'):
                #-- ascii (.txt)
                Ylms = harmonics().from_ascii(REMOVE_FILE)
            elif (REMOVEFORM == 'netCDF4'):
                #-- netCDF4 (.nc)
                Ylms = harmonics().from_netCDF4(REMOVE_FILE)
            elif (REMOVEFORM == 'HDF5'):
                #-- HDF5 (.h5)
                Ylms = harmonics().from_HDF5(REMOVE_FILE)
            elif (REMOVEFORM == 'index'):
                #-- index containing files in data format
                Ylms = harmonics().from_index(REMOVE_FILE, DATAFORM)
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

    #-- calculating GRACE/GRACE-FO error (Wahr et al. 2006)
    #-- output GRACE error file (for both LMAX==MMAX and LMAX != MMAX cases)
    args = (PROC,DREL,DSET,LMAX,order_str,ds_str,atm_str,GRACE_Ylms.month[0],
        GRACE_Ylms.month[-1], suffix[DATAFORM])
    delta_format = '{0}_{1}_{2}_DELTA_CLM_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    DELTA_FILE = delta_format.format(*args)
    #-- check full path of the GRACE directory for delta file
    #-- if file was previously calculated: will read file
    #-- else: will calculate the GRACE/GRACE-FO error
    if (not os.access(os.path.join(grace_dir, DELTA_FILE),os.F_OK)):
        #-- add output delta file to list object
        output_files.append(os.path.join(grace_dir,DELTA_FILE))

        #-- Delta coefficients of GRACE time series (Error components)
        delta_Ylms = harmonics(lmax=LMAX,mmax=MMAX)
        delta_Ylms.clm = np.zeros((LMAX+1,MMAX+1))
        delta_Ylms.slm = np.zeros((LMAX+1,MMAX+1))
        #-- Smoothing Half-Width (CNES is a 10-day solution)
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
                    #-- calculate GRACE Error (Noise of smoothed time-series)
                    #-- With Annual and Semi-Annual Terms
                    val1 = getattr(GRACE_Ylms, csharm)
                    smth = tssmooth(GRACE_Ylms.time, val1[l,m,:], HFWTH=HFWTH)
                    #-- number of smoothed points
                    nsmth = len(smth['data'])
                    tsmth = np.mean(smth['time'])
                    #-- GRACE delta Ylms
                    #-- variance of data-(smoothed+annual+semi)
                    val2 = getattr(delta_Ylms, csharm)
                    val2[l,m] = np.sqrt(np.sum(smth['noise']**2)/nsmth)

        #-- save GRACE/GRACE-FO delta harmonics to file
        delta_Ylms.time = np.copy(tsmth)
        delta_Ylms.month = np.int(nsmth)
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            delta_Ylms.to_ascii(os.path.join(grace_dir,DELTA_FILE))
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            delta_Ylms.to_netCDF4(os.path.join(grace_dir,DELTA_FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            delta_Ylms.to_HDF5(os.path.join(grace_dir,DELTA_FILE))
    else:
        #-- read GRACE/GRACE-FO delta harmonics from file
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            delta_Ylms=harmonics().from_ascii(os.path.join(grace_dir,DELTA_FILE))
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            delta_Ylms=harmonics().from_netCDF4(os.path.join(grace_dir,DELTA_FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            delta_Ylms=harmonics().from_HDF5(os.path.join(grace_dir,DELTA_FILE))
        #-- copy time and number of smoothed fields
        tsmth = np.squeeze(delta_Ylms.time)
        nsmth = np.int(delta_Ylms.month)

    #-- Output spatial data object
    grid = spatial()
    grid.lon = np.copy(kfactor.lon)
    grid.lat = np.copy(kfactor.lat)
    grid.time = np.zeros((nfiles))
    grid.month = np.zeros((nfiles),dtype=np.int)
    grid.data = np.zeros((nlat,nlon,nfiles))
    grid.mask = np.zeros((nlat,nlon,nfiles),dtype=bool)

    #-- Computing plms for converting to spatial domain
    phi = grid.lon[np.newaxis,:]*np.pi/180.0
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))
    #-- square of legendre polynomials truncated to order MMAX
    mm = np.arange(0,MMAX+1)
    PLM2 = PLM[:,mm,:]**2

    #-- dfactor is the degree dependent coefficients
    #-- for converting to centimeters water equivalent (cmwe)
    dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).cmwe

    #-- converting harmonics to truncated, smoothed coefficients in units
    #-- combining harmonics to calculate output spatial fields
    for i,gm in enumerate(GRACE_Ylms.month):
        #-- GRACE/GRACE-FO harmonics for time t
        Ylms = GRACE_Ylms.index(i)
        #-- Remove GIA rate for time
        Ylms.subtract(GIA_Ylms.index(i))
        #-- Remove monthly files to be removed
        Ylms.subtract(remove_Ylms.index(i))
        #-- smooth harmonics and convert to output units
        Ylms.convolve(dfactor*wt)
        #-- convert spherical harmonics to output spatial grid
        grid.data[:,:,i] = harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, MMAX=MMAX, PLM=PLM).T
        #-- copy time variables for month
        grid.time[i] = np.copy(Ylms.time)
        grid.month[i] = np.copy(Ylms.month)
    #-- update spacing and dimensions
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    #-- scale output data with kfactor
    grid = grid.scale(kfactor.data)
    grid.replace_invalid(fill_value, mask=kfactor.mask)

    #-- output monthly files to ascii, netCDF4 or HDF5
    args = (FILENAME,'',unit_str,LMAX,order_str,gw_str,ds_str,
        grid.month[0],grid.month[-1],suffix[DATAFORM])
    FILE=os.path.join(DIRECTORY,file_format.format(*args))
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        grid.to_ascii(FILE, date=True, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4
        grid.to_netCDF4(FILE, date=True, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Spatial Data')
    elif (DATAFORM == 'HDF5'):
        #-- HDF5
        grid.to_HDF5(FILE, date=True, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Spatial Data')
    #-- set the permissions mode of the output files
    os.chmod(FILE, MODE)
    #-- add file to list
    output_files.append(FILE)

    #-- calculate power of scaled GRACE/GRACE-FO data
    scaled_power = grid.sum(power=2.0).power(0.5)
    #-- calculate residual leakage errors
    #-- scaled by ratio of GRACE and synthetic power
    ratio = scaled_power.scale(k_power.power(-1).data)
    error = k_error.scale(ratio.data)

    #-- output monthly error files to ascii, netCDF4 or HDF5
    args = (FILENAME,'ERROR_',unit_str,LMAX,order_str,gw_str,ds_str,
        grid.month[0],grid.month[-1],suffix[DATAFORM])
    FILE = os.path.join(DIRECTORY,file_format.format(*args))
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        error.to_ascii(FILE, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4
        error.to_netCDF4(FILE, date=False, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Scaling Error')
    elif (DATAFORM == 'HDF5'):
        #-- HDF5
        error.to_HDF5(FILE, date=False, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Scaling Error')
    #-- set the permissions mode of the output files
    os.chmod(FILE, MODE)
    #-- add file to list
    output_files.append(FILE)

    #-- Output spatial data object
    delta = spatial()
    delta.lon = np.copy(kfactor.lon)
    delta.lat = np.copy(kfactor.lat)
    delta.time = np.copy(tsmth)
    delta.month = np.copy(nsmth)
    delta.data = np.zeros((nlat,nlon))
    delta.mask = np.zeros((nlat,nlon),dtype=bool)
    #-- calculate scaled spatial error
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
    #-- Multiplying by c/s(phi#m) to get spatial error map
    delta.data[:] = np.sqrt(np.dot(ccos.T,d_cos) + np.dot(ssin.T,d_sin)).T
    #-- update spacing and dimensions
    delta.update_spacing()
    delta.update_extents()
    delta.update_dimensions()

    #-- scale output harmonic errors with kfactor
    delta = delta.scale(kfactor.data)
    delta.replace_invalid(fill_value, mask=kfactor.mask)

    #-- output monthly files to ascii, netCDF4 or HDF5
    args = (FILENAME,'DELTA_',unit_str,LMAX,order_str,gw_str,ds_str,
        grid.month[0],grid.month[-1],suffix[DATAFORM])
    FILE=os.path.join(DIRECTORY,file_format.format(*args))
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        delta.to_ascii(FILE, date=True, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4
        delta.to_netCDF4(FILE, date=True, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Spatial Error')
    elif (DATAFORM == 'HDF5'):
        #-- HDF5
        delta.to_HDF5(FILE, date=True, verbose=VERBOSE,
            units=unit_str, longname=unit_name,
            title='GRACE/GRACE-FO Spatial Error')
    #-- set the permissions mode of the output files
    os.chmod(FILE, MODE)
    #-- add file to list
    output_files.append(FILE)

    #-- return the list of output files
    return output_files

#-- PURPOSE: print a file log for the GRACE/GRACE-FO analysis
#-- lists: the parameter file, the parameters and the output files
def output_log_file(parameters,output_files):
    #-- format: scale_grace_maps_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'scale_grace_maps_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- check if run from entering parameters or from parameter files
    if parameters['PARAMETER_FILE'] is not None:
        #-- print parameter file on top
        print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
            os.path.abspath(parameters['PARAMETER_FILE'])), file=fid)
    else:
        print('PARAMETERS:', file=fid)
    #-- print parameter values sorted alphabetically
    for p in sorted(list(set(parameters.keys())-set(['PARAMETER_FILE']))):
        print('{0}: {1}'.format(p, parameters[p]), file=fid)
    #-- print output files
    print('\n\nOUTPUT FILES:', file=fid)
    for f in output_files:
        print('{0}'.format(f), file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE/GRACE-FO analysis
#-- lists: the parameter file, the parameters and the error
def output_error_log_file(parameters):
    #-- format: scale_grace_maps_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'scale_grace_maps_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- check if run from entering parameters or from parameter files
    if parameters['PARAMETER_FILE'] is not None:
        #-- print parameter file on top
        print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
            os.path.abspath(parameters['PARAMETER_FILE'])), file=fid)
    else:
        print('PARAMETERS:', file=fid)
    #-- print parameter values sorted alphabetically
    for p in sorted(list(set(parameters.keys())-set(['PARAMETER_FILE']))):
        print('{0}: {1}'.format(p, parameters[p]), file=fid)
    #-- print traceback error
    print('\n\nTRACEBACK ERROR:', file=fid)
    traceback.print_exc(file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: open a unique log file adding a numerical instance if existing
def create_unique_logfile(filename):
    #-- split filename into fileBasename and fileExtension
    fileBasename, fileExtension = os.path.splitext(filename)
    #-- create counter to add to the end of the filename if existing
    counter = 1
    while counter:
        try:
            #-- open file descriptor only if the file doesn't exist
            fd = os.open(filename, os.O_CREAT | os.O_EXCL | os.O_RDWR)
        except OSError:
            pass
        else:
            return os.fdopen(fd, 'w+')
        #-- new filename adds counter the between fileBasename and fileExtension
        filename = '{0}_{1:d}{2}'.format(fileBasename, counter, fileExtension)
        counter += 1

#-- PURPOSE: define the analysis for multiprocessing
def define_analysis(f,base_dir,LOVE_NUMBERS=0,REFERENCE=None,LOG=False,
    VERBOSE=False,MODE=0o775):
    #-- keep track of multiprocessing threads
    info(os.path.basename(f))

    #-- variable with parameter definitions
    parameters = {}
    parameters['PARAMETER_FILE'] = f
    #-- Opening parameter file and assigning file ID number (fid)
    fid = open(os.path.expanduser(f), 'r')
    #-- for each line in the file will extract the parameter (name and value)
    for fileline in fid:
        #-- Splitting the input line between parameter name and value
        part = fileline.split()
        #-- filling the parameter definition variable
        parameters[part[0]] = part[1]
    #-- close the parameter file
    fid.close()

    #-- try to run the analysis with listed parameters
    try:
        #-- run scale spatial maps algorithm with parameters
        output_files = scale_grace_maps(base_dir, parameters,
            LOVE_NUMBERS=LOVE_NUMBERS, REFERENCE=REFERENCE,
            VERBOSE=VERBOSE, MODE=MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
        if LOG:#-- write failed job completion log file
            output_error_log_file(parameters)
    else:
        if LOG:#-- write successful job completion log file
            output_log_file(parameters,output_files)

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates scaled spatial maps from
            GRACE/GRACE-FO spherical harmonic coefficients
            """
    )
    #-- command line parameters
    parser.add_argument('parameters',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Parameter files containing specific variables for each analysis')
    #-- number of processes to run in parallel
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to run in parallel')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference','-r',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- Output log file for each job in forms
    #-- scale_grace_maps_run_2002-04-01_PID-00000.log
    #-- scale_grace_maps_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file for each job')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- use parameter files from system arguments listed after the program.
    if (args.np == 0):
        #-- run directly as series if PROCESSES = 0
        #-- for each entered parameter file
        for f in args.parameters:
            define_analysis(f,args.directory,LOVE_NUMBERS=args.love,
                REFERENCE=args.reference, LOG=args.log, VERBOSE=args.verbose,
                MODE=args.mode)
    else:
        #-- run in parallel with multiprocessing Pool
        pool = multiprocessing.Pool(processes=args.np)
        #-- for each entered parameter file
        for f in args.parameters:
            kwds=dict(LOVE_NUMBERS=args.love,REFERENCE=args.reference,
                LOG=args.log,VERBOSE=args.verbose, MODE=args.mode)
            pool.apply_async(define_analysis,args=(f,args.directory),kwds=kwds)
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()

#-- run main program
if __name__ == '__main__':
    main()
