#!/usr/bin/env python
u"""
monte_carlo_degree_one.py
Written by Tyler Sutterley (05/2022)

Calculates degree 1 errors using GRACE coefficients of degree 2 and greater,
    and ocean bottom pressure variations from OMCT/MPIOM in a Monte Carlo scheme

Relation between Geocenter Motion in mm and Normalized Geoid Coefficients
    X = sqrt(3)*a*C11
    Y = sqrt(3)*a*S11
    Z = sqrt(3)*a*C10
where a is the average radius of the Earth (6.371e9 mm)

Load Love Number of Degree 1 for the Center of Figure (CF) Reference Frame:
    kl[1] = -(hl[1]+2*ll[1])/3
Where hl and ll are the displacement Love numbers when the origin is the center
    of mass of the deformed solid Earth
Meaning that kl is defined so that the (l == 1) terms describe the offset
    between the (center of mass of the surface mass + deformed solid Earth)
    and the (center of figure of the deformed solid Earth surface)

Surface density change at any given location is a combination of its
    land and ocean components (after the atmosphere has been removed)

S[theta,phi] = L[theta,phi]*S[theta,phi] + O[theta,phi]*S[theta,phi]
    where S is the surface mass density at colatitude theta and longitude phi
    L is the land function (L==1 if land)
    O is the ocean function (O==1 if ocean)

S can be described as a summation of global spherical harmonics
    with both land and ocean components

If assumed that all degrees are known except the geocenter:
    condition with 2 equations and 3 unknowns (global C10, C11, S11)

The effects of gravitational self-attraction can be considered
    (e.g. from terrestrial hydrology and ice sheet melt)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    --runs X: number of runs in Monte Carlo simulation
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -k X, --kl X: Degree 1 Gravitational Load Love number
    --atm-correction: Apply atmospheric jump correction coefficients
    --pole-tide: Correct for pole tide drift
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
    -F X, --format X: Input/output data format for harmonics files
        ascii
        netCDF4
        HDF5
    --ocean-file X: Index file for ocean model harmonics
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
    --error-file X: Additional error files to use in monte carlo analysis
    --iterative: Iterate degree one solutions
    --fingerprint: Redistribute land-water flux using sea level fingerprints
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --mask X: Land-sea mask for calculating ocean mass and land water flux
    -p, --plot: Create output plots for components and iterations
    --log: Output log of files created for each job
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
        http://www.h5py.org/
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces C20,C21,S21,C22,S22,C30 and C50 with SLR values (if specified)
    time.py: utilities for calculating time operations
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    sea_level_equation.py: pseudo-spectral sea level equation solver
    tssmooth.py: smoothes a time-series using a 13-month Loess-type algorithm
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    harmonics.py: class for processing GRACE/GRACE-FO spherical harmonic data
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: class for reading, writing and processing spatial data
    utilities.py: download and management utilities for files

REFERENCES:
    T C Sutterley, and I Velicogna, "Improved estimates of geocenter
        variability from time-variable gravity and ocean model outputs",
        Remote Sensing, 11(18), 2108, (2019).
        https://doi.org/10.3390/rs11182108

    S Swenson, D Chambers and J Wahr, "Estimating geocenter variations
        from a combination of GRACE and ocean model output,"
        Journal of Geophysical Research: Solid Earth, 113(B08410), (2008).
        https://doi.org/10.1029/2007JB005338

    J Wahr, S C Swenson, and I Velicogna, "Accuracy of GRACE mass estimates",
        Geophysical Research Letters, 33(6), L06401, (2006).
        https://doi.org/10.1029/2005GL025305

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within documentation
        use GIA reference and citation output from GIA read program
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: add uncertainty variable columns in output yaml header
        can use variable loglevels for verbose output
    Updated 11/2021: add GSFC low-degree harmonics
        use gravity_toolkit geocenter class for operations
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: add atmospheric jump corrections to iteration
    Updated 08/2021: reorganize GRACE/GRACE-FO file import
        added option to use additional error files in monte carlo
        expand outputs to include mean harmonics and monte carlo errors
    Updated 07/2021: fix YAML headers for S11 description
        simplified file imports using wrappers in harmonics
        remove choices for argparse processing centers
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options and from gfc for mean files
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 07/2020: using spatial class for reading and operating spatial data
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: added creation date as a global attribute of netCDF4 file
        can replace C30 with coefficients from satellite laser ranging (SLR)
    Updated 06/2019: added parameter LANDMASK for setting the land-sea mask
    Updated 12/2018: added parallel processing with multiprocessing pools
        output all monte carlo iterations to a single netCDF4 file
    Written 11/2018
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import netCDF4
import argparse
import traceback
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import gravity_toolkit.utilities as utilities
from gravity_toolkit.grace_input_months import grace_input_months, \
    read_ecmwf_corrections
from gravity_toolkit.geocenter import geocenter
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
from gravity_toolkit.units import units
from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_love_numbers import load_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.sea_level_equation import sea_level_equation
from gravity_toolkit.tssmooth import tssmooth
from gravity_toolkit.time import grace_to_calendar

#-- PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: model the seasonal component of an initial degree 1 model
#-- using preliminary estimates of annual and semi-annual variations from LWM
#-- as calculated in Chen et al. (1999), doi:10.1029/1998JB900019
#-- NOTE: this is to get an accurate assessment of the land water mass for the
#-- eustatic component (not for the ocean component from GRACE)
def model_seasonal_geocenter(grace_date):
    #-- Annual amplitudes of (Soil Moisture + Snow) geocenter components (mm)
    AAx = 1.28
    AAy = 0.52
    AAz = 3.30
    #-- Annual phase of (Soil Moisture + Snow) geocenter components (degrees)
    APx = 44.0
    APy = 182.0
    APz = 43.0
    #-- Semi-Annual amplitudes of (Soil Moisture + Snow) geocenter components
    SAAx = 0.15
    SAAy = 0.56
    SAAz = 0.50
    #-- Semi-Annual phase of (Soil Moisture + Snow) geocenter components
    SAPx = 331.0
    SAPy = 312.0
    SAPz = 75.0
    #-- calculate each geocenter component from the amplitude and phase
    #-- converting the phase from degrees to radians
    X = AAx*np.sin(2.0*np.pi*grace_date + APx*np.pi/180.0) + \
        SAAx*np.sin(4.0*np.pi*grace_date + SAPx*np.pi/180.0)
    Y = AAy*np.sin(2.0*np.pi*grace_date + APy*np.pi/180.0) + \
        SAAy*np.sin(4.0*np.pi*grace_date + SAPy*np.pi/180.0)
    Z = AAz*np.sin(2.0*np.pi*grace_date + APz*np.pi/180.0) + \
        SAAz*np.sin(4.0*np.pi*grace_date + SAPz*np.pi/180.0)
    DEG1 = geocenter(X=X-X.mean(), Y=Y-Y.mean(), Z=Z-Z.mean())
    return DEG1.from_cartesian()

#-- PURPOSE: calculate the satellite error for a geocenter time-series
def monte_carlo_degree_one(base_dir, PROC, DREL, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
    RUNS=0,
    LOVE_NUMBERS=0,
    LOVE_K1=None,
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    POLE_TIDE=False,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C50=None,
    DATAFORM=None,
    MEAN_FILE=None,
    MEANFORM=None,
    ERROR_FILES=[],
    FINGERPRINT=False,
    EXPANSION=None,
    LANDMASK=None,
    PLOT=False,
    MODE=0o775):

    #-- GRACE/GRACE-FO dataset
    DSET = 'GSM'
    #-- do not import degree 1 coefficients
    DEG1 = ''

    #-- delta coefficients flag for monte carlo run
    delta_str = '_monte_carlo'
    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if MMAX and (MMAX != LMAX) else ''
    #-- atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    #-- ocean model string
    model_str = 'MPIOM' if (DREL == 'RL06') else 'OMCT'
    #-- output flag for using sea level fingerprints
    slf_str = '_SLF' if FINGERPRINT else ''
    #-- output flag for low-degree harmonic replacements
    if SLR_21 in ('CSR','GFZ','GSFC'):
        C21_str = '_w{0}_21'.format(SLR_21)
    else:
        C21_str = ''
    if SLR_22 in ('CSR','GSFC'):
        C22_str = '_w{0}_22'.format(SLR_22)
    else:
        C22_str = ''
    if SLR_C30 in ('GSFC',):
        #-- C30 replacement now default for all solutions
        C30_str = ''
    elif SLR_C30 in ('CSR','GFZ','LARES'):
        C30_str = '_w{0}_C30'.format(SLR_C30)
    else:
        C30_str = ''
    if SLR_C50 in ('CSR','GSFC','LARES'):
        C50_str = '_w{0}_C50'.format(SLR_C50)
    else:
        C50_str = ''
    #-- combine satellite laser ranging flags
    slr_str = ''.join([C21_str,C22_str,C30_str,C50_str])
    #-- suffix for input ascii, netcdf and HDF5 files
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- output directory
    DIRECTORY = os.path.join(base_dir,'geocenter')
    #-- list object of output files for file logs (full path)
    output_files = []

    #-- read load love numbers
    hl,kl,ll = load_love_numbers(EXPANSION, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE='CF')
    #-- set gravitational load love number to a specific value
    if LOVE_K1:
        kl[1] = np.copy(LOVE_K1)
    #-- maximum spherical harmonic order
    if not MMAX:
        MMAX = np.copy(LMAX)

    #-- Earth Parameters
    factors = units(lmax=LMAX).harmonic(hl,kl,ll)
    rho_e = factors.rho_e#-- Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e#-- Average Radius of the Earth [cm]
    l = factors.l
    #-- Factor for converting to Mass SH
    dfactor = factors.cmwe

    #-- Read Smoothed Ocean and Land Functions
    #-- smoothed functions are from the read_ocean_function.py program
    #-- Open the land-sea NetCDF file for reading
    landsea = spatial().from_netCDF4(LANDMASK, date=False, varname='LSMASK')
    #-- degree spacing and grid dimensions
    #-- will create GRACE spatial fields with same dimensions
    dlon,dlat = landsea.spacing
    nlat,nlon = landsea.shape
    #-- spatial parameters in radians
    dphi = dlon*np.pi/180.0
    dth = dlat*np.pi/180.0
    #-- longitude and colatitude in radians
    phi = landsea.lon[np.newaxis,:]*np.pi/180.0
    th = (90.0 - np.squeeze(landsea.lat))*np.pi/180.0
    #-- create land function
    land_function = np.zeros((nlon,nlat),dtype=np.float64)
    #-- extract land function from file
    #-- combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data.T >= 1) & (landsea.data.T <= 3))
    land_function[indx,indy] = 1.0
    #-- calculate ocean function from land function
    ocean_function = 1.0 - land_function

    #-- Calculating Legendre Polynomials using Holmes and Featherstone relation
    PLM,dPLM = plm_holmes(LMAX,np.cos(th))

    #-- calculate spherical harmonics of ocean function to degree 1
    #-- mass is equivalent to 1 cm ocean height change
    #-- eustatic ratio = -land total/ocean total
    ocean_Ylms = gen_stokes(ocean_function, landsea.lon, landsea.lat,
        UNITS=1, LMIN=0, LMAX=1, LOVE=(hl,kl,ll), PLM=PLM[:2,:2,:])

    #-- Gaussian Smoothing (Jekeli, 1981)
    if (RAD != 0):
        wt = 2.0*np.pi*gauss_weights(RAD,LMAX)
    else:
        #-- else = 1
        wt = np.ones((LMAX+1))

    #-- reading GRACE months for input date range
    #-- replacing low-degree harmonics with SLR values if specified
    #-- correcting for Pole-Tide drift if specified
    #-- atmospheric jumps will be corrected externally if specified
    Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX,
        SLR_21=SLR_21, SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C50=SLR_C50,
        POLE_TIDE=POLE_TIDE, ATM=False, MODEL_DEG1=False)
    #-- create harmonics object from GRACE/GRACE-FO data
    GSM_Ylms = harmonics().from_dict(Ylms)
    #-- use a mean file for the static field to remove
    if MEAN_FILE:
        #-- read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = harmonics().from_file(MEAN_FILE,format=MEANFORM,date=False)
        #-- remove the input mean
        GSM_Ylms.subtract(mean_Ylms)
    else:
        GSM_Ylms.mean(apply=True)
    #-- filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        #-- destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GSM_Ylms = GSM_Ylms.destripe()
    else:
        #-- using standard GRACE/GRACE-FO harmonics
        ds_str = ''
    #-- full path to directory for specific GRACE/GRACE-FO product
    GSM_Ylms.directory = Ylms['directory']
    #-- GRACE dates
    tdec = np.copy(GSM_Ylms.time)
    months = np.copy(GSM_Ylms.month)
    #-- number of months considered
    n_files = len(GSM_Ylms.month)

    #-- input GIA spherical harmonic datafiles
    GIA_Ylms_rate = read_GIA_model(GIA_FILE,GIA=GIA,LMAX=LMAX,MMAX=MMAX)
    gia_str = '_{0}'.format(GIA_Ylms_rate['title']) if GIA else ''
    #-- calculate the monthly mass change from GIA
    GIA_Ylms = GSM_Ylms.zeros_like()
    GIA_Ylms.time[:] = np.copy(GSM_Ylms.time)
    GIA_Ylms.month[:] = np.copy(GSM_Ylms.month)
    #-- monthly GIA calculated by gia_rate*time elapsed
    #-- finding change in GIA each month
    for t in range(n_files):
        GIA_Ylms.clm[:,:,t] = GIA_Ylms_rate['clm']*(GIA_Ylms.time[t]-2003.3)
        GIA_Ylms.slm[:,:,t] = GIA_Ylms_rate['slm']*(GIA_Ylms.time[t]-2003.3)
    #-- save geocenter coefficients of monthly GIA variability
    gia = geocenter().from_harmonics(GIA_Ylms)

    #-- read atmospheric jump corrections from Fagiolini et al. (2015)
    ATM_Ylms = GSM_Ylms.zeros_like()
    ATM_Ylms.time[:] = np.copy(GSM_Ylms.time)
    ATM_Ylms.month[:] = np.copy(GSM_Ylms.month)
    if ATM:
        atm_corr = read_ecmwf_corrections(base_dir,LMAX,ATM_Ylms.month)
        ATM_Ylms.clm[:,:,:] = np.copy(atm_corr['clm'])
        ATM_Ylms.slm[:,:,:] = np.copy(atm_corr['slm'])
        #-- removing the mean of the atmospheric jump correction coefficients
        ATM_Ylms.mean(apply=True)
    #-- truncate to degree and order LMAX/MMAX
    ATM_Ylms = ATM_Ylms.truncate(lmax=LMAX, mmax=MMAX)
    #-- save geocenter coefficients of the atmospheric jump corrections
    atm = geocenter().from_harmonics(ATM_Ylms)

    #-- input spherical harmonic datafiles to be used in monte carlo
    error_Ylms = []
    #-- for each file to be removed
    for ERROR_FILE in ERROR_FILES:
        #-- file in ascii, netCDF4 or HDF5 formats
        Ylms = harmonics().from_file(ERROR_FILE, format=DATAFORM)
        #-- truncate to degree and order and append to list
        error_Ylms.append(Ylms.truncate(lmax=LMAX, mmax=MMAX))

    #-- calculating GRACE/GRACE-FO error (Wahr et al. 2006)
    #-- output GRACE error file (for both LMAX==MMAX and LMAX != MMAX cases)
    args = (PROC,DREL,DSET,LMAX,order_str,ds_str,atm_str,GSM_Ylms.month[0],
        GSM_Ylms.month[-1], suffix[DATAFORM])
    delta_format = '{0}_{1}_{2}_DELTA_CLM_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    DELTA_FILE = os.path.join(GSM_Ylms.directory,delta_format.format(*args))
    #-- check full path of the GRACE directory for delta file
    #-- if file was previously calculated: will read file
    #-- else: will calculate the GRACE/GRACE-FO error
    if not os.access(DELTA_FILE, os.F_OK):
        #-- add output delta file to list object
        output_files.append(DELTA_FILE)

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
                    val1 = getattr(GSM_Ylms, csharm)
                    smth = tssmooth(tdec, val1[l,m,:], HFWTH=HFWTH)
                    #-- number of smoothed points
                    nsmth = len(smth['data'])
                    tsmth = np.mean(smth['time'])
                    #-- GRACE delta Ylms
                    #-- variance of data-(smoothed+annual+semi)
                    val2 = getattr(delta_Ylms, csharm)
                    val2[l,m] = np.sqrt(np.sum(smth['noise']**2)/nsmth)

        #-- save GRACE/GRACE-FO delta harmonics to file
        delta_Ylms.time = np.copy(tsmth)
        delta_Ylms.month = np.int64(nsmth)
        delta_Ylms.to_file(DELTA_FILE,format=DATAFORM)
    else:
        #-- read GRACE/GRACE-FO delta harmonics from file
        delta_Ylms = harmonics().from_file(DELTA_FILE,format=DATAFORM)
        #-- truncate GRACE/GRACE-FO delta clm and slm to d/o LMAX/MMAX
        delta_Ylms = delta_Ylms.truncate(lmax=LMAX, mmax=MMAX)
        tsmth = np.squeeze(delta_Ylms.time)
        nsmth = np.int64(delta_Ylms.month)

    #-- Calculating cos/sin of phi arrays
    #-- output [m,phi]
    m = GSM_Ylms.m[:, np.newaxis]
    #-- Integration factors (solid angle)
    int_fact = np.sin(th)*dphi*dth
    #-- Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))

    #-- Legendre polynomials for degree 1
    P10 = np.squeeze(PLM[1,0,:])
    P11 = np.squeeze(PLM[1,1,:])
    #-- PLM for spherical harmonic degrees 2+
    #-- converted into mass and smoothed if specified
    plmout = np.zeros((LMAX+1,MMAX+1,nlat))
    for l in range(1,LMAX+1):
        m = np.arange(0,np.min([l,MMAX])+1)
        #-- convert to smoothed coefficients of mass
        #-- Convolving plms with degree dependent factor and smoothing
        plmout[l,m,:] = PLM[l,m,:]*dfactor[l]*wt[l]

    #-- Initializing 3x3 I-Parameter matrix
    IMAT = np.zeros((3,3))
    #-- Calculating I-Parameter matrix by integrating over latitudes
    #-- I-Parameter matrix accounts for the fact that the GRACE data only
    #-- includes spherical harmonic degrees greater than or equal to 2
    for i in range(0,nlat):
        #-- C10: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,0] += np.sum(int_fact[i]*P10[i]*ccos[0,:]*ocean_function[:,i]*P10[i]*ccos[0,:])/(4.0*np.pi)
        IMAT[1,0] += np.sum(int_fact[i]*P10[i]*ccos[0,:]*ocean_function[:,i]*P11[i]*ccos[1,:])/(4.0*np.pi)
        IMAT[2,0] += np.sum(int_fact[i]*P10[i]*ccos[0,:]*ocean_function[:,i]*P11[i]*ssin[1,:])/(4.0*np.pi)
        #-- C11: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,1] += np.sum(int_fact[i]*P11[i]*ccos[1,:]*ocean_function[:,i]*P10[i]*ccos[0,:])/(4.0*np.pi)
        IMAT[1,1] += np.sum(int_fact[i]*P11[i]*ccos[1,:]*ocean_function[:,i]*P11[i]*ccos[1,:])/(4.0*np.pi)
        IMAT[2,1] += np.sum(int_fact[i]*P11[i]*ccos[1,:]*ocean_function[:,i]*P11[i]*ssin[1,:])/(4.0*np.pi)
        #-- S11: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,2] += np.sum(int_fact[i]*P11[i]*ssin[1,:]*ocean_function[:,i]*P10[i]*ccos[0,:])/(4.0*np.pi)
        IMAT[1,2] += np.sum(int_fact[i]*P11[i]*ssin[1,:]*ocean_function[:,i]*P11[i]*ccos[1,:])/(4.0*np.pi)
        IMAT[2,2] += np.sum(int_fact[i]*P11[i]*ssin[1,:]*ocean_function[:,i]*P11[i]*ssin[1,:])/(4.0*np.pi)

    #-- get seasonal variations of an initial geocenter correction
    #-- for use in the land water mass calculation
    seasonal_geocenter = model_seasonal_geocenter(tdec)

    #-- degree 1 iterations for each monte carlo run
    iteration = geocenter()
    iteration.C10 = np.zeros((n_files,RUNS))
    iteration.C11 = np.zeros((n_files,RUNS))
    iteration.S11 = np.zeros((n_files,RUNS))
    #-- for each monte carlo iteration
    for n_iter in range(0, RUNS):
        #-- calculate non-iterated terms for each file (G-matrix parameters)
        for t in range(n_files):
            #-- calculate uncertainty for time t and each degree/order
            Ylms = harmonics(lmax=LMAX, mmax=MMAX)
            Ylms.clm = (1.0-2.0*np.random.rand(LMAX+1,MMAX+1))*delta_Ylms.clm
            Ylms.slm = (1.0-2.0*np.random.rand(LMAX+1,MMAX+1))*delta_Ylms.slm
            #-- add additional uncertainty terms
            for eYlms in error_Ylms:
                Ylms.clm += (1.0-2.0*np.random.rand(LMAX+1,MMAX+1))*eYlms.clm
                Ylms.slm += (1.0-2.0*np.random.rand(LMAX+1,MMAX+1))*eYlms.slm

            #-- Removing monthly GIA signal and atmospheric correction
            GRACE_Ylms = GSM_Ylms.index(t)
            GRACE_Ylms.subtract(GIA_Ylms.index(t))
            GRACE_Ylms.subtract(ATM_Ylms.index(t))

            #-- G matrix calculates the GRACE ocean mass variations
            G = geocenter()
            G.C10 = 0.0
            G.C11 = 0.0
            G.S11 = 0.0
            #-- calculate non-iterated terms (G-matrix parameters)
            #-- calculate geocenter component of ocean mass with GRACE
            #-- allocate for product of grace and legendre polynomials
            pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
            psin = np.zeros((MMAX+1, nlat))#-[m,lat]
            #-- Summing product of plms and c/slms over all SH degrees >= 2
            for i in range(0, nlat):
                l = np.arange(2,LMAX+1)
                pcos[:,i] = np.sum(plmout[l,:,i]*(GRACE_Ylms.clm[l,:]+Ylms.clm[l,:]), axis=0)
                psin[:,i] = np.sum(plmout[l,:,i]*(GRACE_Ylms.slm[l,:]+Ylms.slm[l,:]), axis=0)
            #-- Multiplying by c/s(phi#m) to get surface density in cmH2Oeq (lon,lat)
            #-- ccos/ssin are mXphi, pcos/psin are mXtheta: resultant matrices are phiXtheta
            #-- The summation over spherical harmonic order is in this multiplication
            rmass = np.dot(np.transpose(ccos),pcos) + np.dot(np.transpose(ssin),psin)
            #-- calculate G matrix parameters through a summation of each latitude
            for i in range(0,nlat):
                #-- summation of integration factors, Legendre polynomials,
                #-- (convolution of order and harmonics) and the ocean mass at t
                G.C10 += np.sum(int_fact[i]*P10[i]*ccos[0,:]*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)
                G.C11 += np.sum(int_fact[i]*P11[i]*ccos[1,:]*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)
                G.S11 += np.sum(int_fact[i]*P11[i]*ssin[1,:]*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)

            #-- seasonal component of geocenter variation for land water
            GSM_Ylms.clm[1,0,t] = seasonal_geocenter.C10[t]
            GSM_Ylms.clm[1,1,t] = seasonal_geocenter.C11[t]
            GSM_Ylms.slm[1,1,t] = seasonal_geocenter.S11[t]
            #-- Removing monthly GIA signal and atmospheric correction
            GRACE_Ylms = GSM_Ylms.index(t)
            GRACE_Ylms.subtract(GIA_Ylms.index(t))
            GRACE_Ylms.subtract(ATM_Ylms.index(t))

            #-- allocate for product of grace and legendre polynomials
            pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
            psin = np.zeros((MMAX+1, nlat))#-[m,lat]
            #-- Summing product of plms and c/slms over all SH degrees
            for i in range(0, nlat):
                #-- for land water: use an initial seasonal geocenter estimate
                #-- from Chen et al. (1999)
                l = np.arange(1,LMAX+1)
                pcos[:,i] = np.sum(plmout[l,:,i]*(GRACE_Ylms.clm[l,:]+Ylms.clm[l,:]), axis=0)
                psin[:,i] = np.sum(plmout[l,:,i]*(GRACE_Ylms.slm[l,:]+Ylms.slm[l,:]), axis=0)

            #-- Multiplying by c/s(phi#m) to get surface density in cm w.e. (lonxlat)
            #-- this will be a spatial field similar to outputs from stokes_combine.py
            #-- ccos/ssin are mXphi, pcos/psin are mXtheta: resultant matrices are phiXtheta
            #-- The summation over spherical harmonic order is in this multiplication
            lmass = np.dot(np.transpose(ccos),pcos) + np.dot(np.transpose(ssin),psin)

            #-- use sea level fingerprints or eustatic from GRACE land components
            if FINGERPRINT:
                #-- calculate total sea level fingerprint for eustatic component
                #-- steps to calculate sea level from GRACE land-water change:
                #-- 1) calculate total land mass at time t (GRACE*land function)
                #-- NOTE: this is an unscaled GRACE estimate that uses the
                #-- buffered land function when solving the sea-level equation.
                #-- possible improvement using scaled estimate with real coastlines
                land_Ylms = gen_stokes(land_function*lmass, landsea.lon,
                    landsea.lat, UNITS=1, LMIN=0, LMAX=EXPANSION, LOVE=(hl,kl,ll))
                #-- 2) calculate sea level fingerprints of land mass at time t
                #-- use maximum of 3 iterations for computational efficiency
                sea_level = sea_level_equation(land_Ylms.clm, land_Ylms.slm,
                    landsea.lon, landsea.lat, land_function, LMAX=EXPANSION,
                    LOVE=(hl,kl,ll), BODY_TIDE_LOVE=0, FLUID_LOVE=0, ITERATIONS=3,
                    POLAR=True, FILL_VALUE=0)
                #-- 3) convert sea level fingerprints into spherical harmonics
                slf_Ylms = gen_stokes(sea_level, landsea.lon, landsea.lat,
                    UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:], LOVE=(hl,kl,ll))
                #-- 4) convert the slf degree 1 harmonics to mass with dfactor
                eustatic = geocenter().from_harmonics(slf_Ylms).scale(dfactor[1])
            else:
                #-- steps to calculate eustatic component from GRACE land-water change:
                #-- 1) calculate total mass of 1 cm of ocean height (calculated above)
                #-- 2) calculate total land mass at time t (GRACE*land function)
                #-- NOTE: possible improvement using the sea-level equation to solve
                #-- for the spatial pattern of sea level from the land water mass
                land_Ylms = gen_stokes(lmass*land_function, landsea.lon, landsea.lat,
                    UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:], LOVE=(hl,kl,ll))
                #-- 3) calculate ratio between the total land mass and the total mass
                #-- of 1 cm of ocean height (negative as positive land = sea level drop)
                #-- this converts the total land change to ocean height change
                eustatic_ratio = -land_Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
                #-- 4) scale degree one coefficients of ocean function with ratio
                #-- and convert the eustatic degree 1 harmonics to mass with dfactor
                scale_factor = eustatic_ratio*dfactor[1]
                eustatic = geocenter().from_harmonics(ocean_Ylms).scale(scale_factor)

            #-- eustatic coefficients of degree 1
            CMAT = np.array([eustatic.C10,eustatic.C11,eustatic.S11])
            #-- G Matrix for time t
            GMAT = np.array([G.C10, G.C11, G.S11])
            #-- calculate inversion for degree 1 solutions
            #-- this is mathematically equivalent to an iterative procedure
            #-- whereby the initial degree one coefficients are used to update
            #-- the G Matrix until (C10, C11, S11) converge
            DMAT = np.dot(np.linalg.inv(IMAT), (CMAT-GMAT))
            #-- could also use pseudo-inverse in least-squares
            #DMAT = np.linalg.lstsq(IMAT,(CMAT-GMAT),rcond=-1)[0]
            #-- save geocenter for iteration and time t after restoring GIA+ATM
            iteration.C10[t,n_iter] = DMAT[0]+gia.C10[t]+atm.C10[t]
            iteration.C11[t,n_iter] = DMAT[1]+gia.C11[t]+atm.C11[t]
            iteration.S11[t,n_iter] = DMAT[2]+gia.S11[t]+atm.S11[t]
        #-- remove mean of each solution for iteration
        iteration.C10[:,n_iter] -= iteration.C10[:,n_iter].mean()
        iteration.C11[:,n_iter] -= iteration.C11[:,n_iter].mean()
        iteration.S11[:,n_iter] -= iteration.S11[:,n_iter].mean()

    #-- calculate mean degree one time series through all iterations
    MEAN = geocenter()
    MEAN.C10 = np.mean(iteration.C10,axis=1)
    MEAN.C11 = np.mean(iteration.C11,axis=1)
    MEAN.S11 = np.mean(iteration.S11,axis=1)

    #-- calculate RMS off of mean time series
    RMS = geocenter()
    RMS.C10 = np.zeros((n_files))
    RMS.C11 = np.zeros((n_files))
    RMS.S11 = np.zeros((n_files))
    for t in range(n_files):
        RMS.C10[t] = np.sqrt(np.sum((iteration.C10[t,:]-MEAN.C10[t])**2)/RUNS)
        RMS.C11[t] = np.sqrt(np.sum((iteration.C11[t,:]-MEAN.C11[t])**2)/RUNS)
        RMS.S11[t] = np.sqrt(np.sum((iteration.S11[t,:]-MEAN.S11[t])**2)/RUNS)

    #-- Convert inverted solutions into fully normalized spherical harmonics
    #-- for each of the geocenter solutions (C10, C11, S11)
    DEG1 = MEAN.scale(1.0/dfactor[1])
    #-- convert estimated monte carlo errors into fully normalized harmonics
    ERROR = RMS.scale(1.0/dfactor[1])

    #-- output degree 1 coefficients
    file_format = '{0}_{1}_{2}{3}{4}{5}{6}{7}.{8}'
    output_format = ('{0:11.4f}{1:14.6e}{2:14.6e}{3:14.6e}'
        '{4:14.6e}{5:14.6e}{6:14.6e} {7:03d}\n')
    #-- public file format in fully normalized spherical harmonics
    #-- local version with all descriptor flags
    a1=(PROC,DREL,model_str,slf_str,'',gia_str,delta_str,ds_str,'txt')
    FILE1=os.path.join(DIRECTORY,file_format.format(*a1))
    fid1 = open(FILE1,'w')
    #-- print headers for cases with and without dealiasing
    print_header(fid1)
    print_harmonic(fid1,kl[1])
    print_global(fid1,PROC,DREL,model_str.replace('_',' '),GIA_Ylms_rate,
        SLR_C20,SLR_21,months)
    print_variables(fid1,'single precision','fully normalized')
    #-- for each GRACE/GRACE-FO month
    for t,mon in enumerate(months):
        #-- output geocenter coefficients to file
        fid1.write(output_format.format(tdec[t],
            DEG1.C10[t],DEG1.C11[t],DEG1.S11[t],
            ERROR.C10[t],ERROR.C11[t],ERROR.S11[t],mon))
    #-- close the output file
    fid1.close()
    #-- set the permissions mode of the output file
    os.chmod(FILE1, MODE)
    output_files.append(FILE1)

    #-- output all degree 1 coefficients as a netCDF4 file
    a2=(PROC,DREL,model_str,slf_str,'',gia_str,delta_str,ds_str,'nc')
    FILE2 = os.path.join(DIRECTORY,file_format.format(*a2))
    fileID = netCDF4.Dataset(FILE2,'w',format="NETCDF4")
    #-- Defining the NetCDF4 dimensions
    fileID.createDimension('run', RUNS)
    fileID.createDimension('time', n_files)
    #-- defining the NetCDF4 variables
    nc = {}
    nc['time'] = fileID.createVariable('time',tdec.dtype,('time',))
    nc['month'] = fileID.createVariable('month',months.dtype,('time',))
    nc['C10'] = fileID.createVariable('C10',iteration.C10.dtype,
        ('time','run',), zlib=True)
    nc['C11'] = fileID.createVariable('C11',iteration.C11.dtype,
        ('time','run',), zlib=True)
    nc['S11'] = fileID.createVariable('S11',iteration.S11.dtype,
        ('time','run',), zlib=True)
    #-- filling NetCDF4 variables
    nc['time'][:] = tdec[:].copy()
    nc['month'][:] = months[:].copy()
    nc['C10'][:] = iteration.C10[:,:]/dfactor[1]
    nc['C11'][:] = iteration.C11[:,:]/dfactor[1]
    nc['S11'][:] = iteration.S11[:,:]/dfactor[1]
    #-- defining the NetCDF4 attributes
    nc['time'].units = 'years'
    nc['time'].long_name = 'Date_in_Decimal_Years'
    nc['month'].long_name = 'GRACE_month'
    nc['month'].units = 'months since 2001-12-01'
    nc['month'].calendar = 'standard'
    nc['C10'].units = 'fully_normalized'
    nc['C10'].long_name = 'cosine_spherical_harmonic_of_degree_1,_order_0'
    nc['C11'].units = 'fully_normalized'
    nc['C11'].long_name = 'cosine_spherical_harmonic_of_degree_1,_order_1'
    nc['S11'].units = 'fully_normalized'
    nc['S11'].long_name = 'sine_spherical_harmonic_of_degree_1,_order_1'
    #-- define global attributes
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
    #-- close the output file
    fileID.close()
    #-- set the permissions mode of the output file
    os.chmod(FILE2, MODE)
    output_files.append(FILE2)

    #-- create plot showing monte carlo iterations
    if PLOT:
        #-- 3 row plot (C10, C11 and S11)
        ax = {}
        fig,(ax[0],ax[1],ax[2])=plt.subplots(nrows=3,sharex=True,figsize=(6,9))
        #-- show solutions for each iteration
        plot_colors = iter(cm.rainbow(np.linspace(0,1,RUNS)))
        for j in range(n_iter):
            color_j = next(plot_colors)
            #-- C10, C11 and S11
            ax[0].plot(months,10.0*iteration.C10[:,j],color=color_j)
            ax[1].plot(months,10.0*iteration.C11[:,j],color=color_j)
            ax[2].plot(months,10.0*iteration.S11[:,j],color=color_j)
        #-- mean C10, C11 and S11
        ax[0].plot(months,10.0*MEAN.C10,color='k',lw=1.5)
        ax[1].plot(months,10.0*MEAN.C11,color='k',lw=1.5)
        ax[2].plot(months,10.0*MEAN.S11,color='k',lw=1.5)
        #-- labels and set limits
        ax[0].set_ylabel('mm', fontsize=14)
        ax[1].set_ylabel('mm', fontsize=14)
        ax[2].set_ylabel('mm', fontsize=14)
        ax[2].set_xlabel('Grace Month', fontsize=14)
        ax[2].set_xlim(np.floor(months[0]/10.)*10.,np.ceil(months[-1]/10.)*10.)
        ax[2].xaxis.set_minor_locator(MultipleLocator(5))
        ax[2].xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
        #-- add axis labels and adjust font sizes for axis ticks
        for i,lbl in enumerate(['C10','C11','S11']):
            #-- axis label
            ax[i].add_artist(AnchoredText(lbl, pad=0.0, frameon=False, loc=2,
                prop=dict(size=16,weight='bold')))
            #-- axes tick adjustments
            for tick in ax[i].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax[i].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            #-- adjust ticks
            ax[i].get_xaxis().set_tick_params(which='both', direction='in')
            ax[i].get_yaxis().set_tick_params(which='both', direction='in')
        #-- adjust locations of subplots and save to file
        fig.subplots_adjust(left=0.12,right=0.94,bottom=0.06,top=0.98,hspace=0.1)
        args = (PROC,DREL,model_str,ds_str)
        FILE = 'Geocenter_Monte_Carlo_{0}_{1}_{2}{3}.pdf'.format(*args)
        plt.savefig(os.path.join(DIRECTORY,FILE), format='pdf')
        plt.clf()
        #-- set the permissions mode of the output files
        os.chmod(os.path.join(DIRECTORY,FILE), MODE)
        output_files.append(os.path.join(DIRECTORY,FILE))

    #-- return the list of output files
    return output_files

#-- PURPOSE: print YAML header to top of file
def print_header(fid):
    #-- print header
    fid.write('{0}:\n'.format('header'))
    #-- data dimensions
    fid.write('  {0}:\n'.format('dimensions'))
    fid.write('    {0:22}: {1:d}\n'.format('degree',1))
    fid.write('    {0:22}: {1:d}\n'.format('order',1))
    fid.write('\n')

#-- PURPOSE: print spherical harmonic attributes to YAML header
def print_harmonic(fid,kl):
    #-- non-standard attributes
    fid.write('  {0}:\n'.format('non-standard_attributes'))
    #-- load love number
    fid.write('    {0:22}:\n'.format('love_number'))
    long_name = 'Gravitational Load Love Number of Degree 1 (k1)'
    fid.write('      {0:20}: {1}\n'.format('long_name',long_name))
    fid.write('      {0:20}: {1:0.3f}\n'.format('value',kl))
    #-- data format
    data_format = '(f11.4,3e14.6,i4)'
    fid.write('    {0:22}: {1}\n'.format('formatting_string',data_format))
    fid.write('\n')

#-- PURPOSE: print global attributes to YAML header
def print_global(fid,PROC,DREL,MODEL,GIA,SLR,S21,month):
    fid.write('  {0}:\n'.format('global_attributes'))
    MISSION = dict(RL05='GRACE',RL06='GRACE/GRACE-FO')
    title = '{0} Geocenter Coefficients {1} {2}'.format(MISSION[DREL],PROC,DREL)
    fid.write('    {0:22}: {1}\n'.format('title',title))
    summary = []
    summary.append(('Geocenter coefficients derived from {0} mission '
        'measurements and {1} ocean model outputs.').format(MISSION[DREL],MODEL))
    summary.append(('  These coefficients represent the largest-scale '
        'variability of hydrologic, cryospheric, and solid Earth '
        'processes.  In addition, the coefficients represent the '
        'atmospheric and oceanic processes not captured in the {0} {1} '
        'de-aliasing product.').format(MISSION[DREL],DREL))
    #-- get GIA parameters
    summary.append(('  Glacial Isostatic Adjustment (GIA) estimates from '
        '{0} have been restored.').format(GIA['citation']))
    if (DREL == 'RL05'):
        summary.append(('  ECMWF corrections from Fagiolini et al. (2015) have '
        'been restored.'))
    fid.write('    {0:22}: {1}\n'.format('summary',''.join(summary)))
    project = []
    project.append('NASA Gravity Recovery And Climate Experiment (GRACE)')
    project.append('GRACE Follow-On (GRACE-FO)') if (DREL == 'RL06') else None
    fid.write('    {0:22}: {1}\n'.format('project',', '.join(project)))
    keywords = []
    keywords.append('GRACE')
    keywords.append('GRACE-FO') if (DREL == 'RL06') else None
    # keywords.append('Level-2')
    keywords.append('Spherical Harmonic Model')
    keywords.append('Gravitational Field')
    keywords.append('Geopotential')
    keywords.append('Time Variable Gravity')
    keywords.append('Mass Transport')
    keywords.append('Satellite Geodesy')
    fid.write('    {0:22}: {1}\n'.format('keywords',', '.join(keywords)))
    vocabulary = 'NASA Global Change Master Directory (GCMD) Science Keywords'
    fid.write('    {0:22}: {1}\n'.format('keywords_vocabulary',vocabulary))
    hist = '{0} Level-3 Data created at UC Irvine'.format(MISSION[DREL])
    fid.write('    {0:22}: {1}\n'.format('history',hist))
    src = 'An inversion using {0} measurements and {1} ocean model outputs.'
    args = (MISSION[DREL],MODEL,DREL)
    fid.write('    {0:22}: {1}\n'.format('source',src.format(*args)))
    # fid.write('    {0:22}: {1}\n'.format('platform','GRACE-A, GRACE-B'))
    # vocabulary = 'NASA Global Change Master Directory platform keywords'
    # fid.write('    {0:22}: {1}\n'.format('platform_vocabulary',vocabulary))
    # fid.write('    {0:22}: {1}\n'.format('instrument','ACC,KBR,GPS,SCA'))
    # vocabulary = 'NASA Global Change Master Directory instrument keywords'
    # fid.write('    {0:22}: {1}\n'.format('instrument_vocabulary',vocabulary))
    fid.write('    {0:22}: {1:d}\n'.format('processing_level',3))
    ack = []
    ack.append(('Work was supported by an appointment to the NASA Postdoctoral '
        'Program at NASA Goddard Space Flight Center, administered by '
        'Universities Space Research Association under contract with NASA'))
    ack.append('GRACE is a joint mission of NASA (USA) and DLR (Germany)')
    if (DREL == 'RL06'):
        ack.append('GRACE-FO is a joint mission of NASA (USA) and GFZ (Germany)')
    fid.write('    {0:22}: {1}\n'.format('acknowledgement','.  '.join(ack)))
    PRODUCT_VERSION = 'Release-{0}'.format(DREL[2:])
    fid.write('    {0:22}: {1}\n'.format('product_version',PRODUCT_VERSION))
    fid.write('    {0:22}:\n'.format('references'))
    reference = []
    #-- geocenter citations
    reference.append(('T. C. Sutterley, and I. Velicogna, "Improved estimates '
        'of geocenter variability from time-variable gravity and ocean model '
        'outputs", Remote Sensing, 11(18), 2108, (2019). '
        'https://doi.org/10.3390/rs11182108'))
    reference.append(('S. C. Swenson, D. P. Chambers, and J. Wahr, "Estimating '
        'geocenter variations from a combination of GRACE and ocean model '
        'output", Journal of Geophysical Research - Solid Earth, 113(B08410), '
        '(2008). https://doi.org/10.1029/2007JB005338'))
    #-- GIA citation
    reference.append(GIA['reference'])
    #-- ECMWF jump corrections citation
    if (DREL == 'RL05'):
        reference.append(('E. Fagiolini, F. Flechtner, M. Horwath, H. Dobslaw, '
            '''"Correction of inconsistencies in ECMWF's operational '''
            '''analysis data during de-aliasing of GRACE gravity models", '''
            'Geophysical Journal International, 202(3), 2150, (2015). '
            'https://doi.org/10.1093/gji/ggv276'))
    #-- SLR citation for a given solution
    if (SLR == 'CSR'):
        reference.append(('M. Cheng, B. D. Tapley, and J. C. Ries, '
            '''"Deceleration in the Earth's oblateness", Journal of '''
            'Geophysical Research: Solid Earth, 118(2), 740-747, (2013). '
            'https://doi.org/10.1002/jgrb.50058'))
    elif (SLR == 'GSFC'):
        reference.append(('B. D. Loomis, K. E. Rachlin, and S. B. Luthcke, '
            '"Improved Earth Oblateness Rate Reveals Increased Ice Sheet Losses '
            'and Mass-Driven Sea Level Rise", Geophysical Research Letters, '
            '46(12), 6910-6917, (2019). https://doi.org/10.1029/2019GL082929'))
        reference.append(('B. D. Loomis, K. E. Rachlin, D. N. Wiese, '
            'F. W. Landerer, and S. B. Luthcke, "Replacing GRACE/GRACE-FO C30 '
            'with satellite laser ranging: Impacts on Antarctic Ice Sheet mass '
            'change", Geophysical Research Letters, 47(3), (2020). '
            'https://doi.org/10.1029/2019GL085488'))
    elif (SLR == 'GFZ'):
        reference.append(('R. Koenig, P. Schreiner, and C. Dahle, "Monthly '
            'estimates of C(2,0) generated by GFZ from SLR satellites based '
            'on GFZ GRACE/GRACE-FO RL06 background models." V. 1.0. GFZ Data '
            'Services, (2019). http://doi.org/10.5880/GFZ.GRAVIS_06_C20_SLR'))
    if (S21 == 'CSR'):
        reference.append(('M. Cheng, J. C. Ries, and B. D. Tapley, '
            '''"Variations of the Earth's figure axis from satellite laser '''
            'ranging and GRACE", Journal of Geophysical Research: Solid Earth, '
            '116, B01409, (2011). https://doi.org/10.1029/2010JB000850'))
    elif (S21 == 'GFZ'):
        reference.append(('C. Dahle and M. Murboeck, "Post-processed '
            'GRACE/GRACE-FO Geopotential GSM Coefficients GFZ RL06 '
            '(Level-2B Product)." V. 0002. GFZ Data Services, (2019). '
            'http://doi.org/10.5880/GFZ.GRAVIS_06_L2B'))
    #-- print list of references
    for ref in reference:
        fid.write('      - {0}\n'.format(ref))
    creators = 'Tyler C. Sutterley and Isabella Velicogna'
    fid.write('    {0:22}: {1}\n'.format('creator_name', creators))
    emails = 'tsutterl@uw.edu and isabella@uci.edu'
    fid.write('    {0:22}: {1}\n'.format('creator_email', emails))
    url = 'https://www.ess.uci.edu/~velicogna/index.html'
    fid.write('    {0:22}: {1}\n'.format('creator_url', url))
    fid.write('    {0:22}: {1}\n'.format('creator_type', 'group'))
    inst = 'University of Washington; University of California, Irvine'
    fid.write('    {0:22}: {1}\n'.format('creator_institution',inst))
    #-- date range and date created
    calendar_year,calendar_month = grace_to_calendar(month)
    start_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[0],calendar_month[0])
    fid.write('    {0:22}: {1}\n'.format('time_coverage_start', start_time))
    end_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[-1],calendar_month[-1])
    fid.write('    {0:22}: {1}\n'.format('time_coverage_end', end_time))
    today = time.strftime('%Y-%m-%d',time.localtime())
    fid.write('    {0:22}: {1}\n'.format('date_created', today))
    fid.write('\n')

#-- PURPOSE: print variable descriptions to YAML header
def print_variables(fid,data_precision,data_units):
    #-- variables
    fid.write('  {0}:\n'.format('variables'))
    #-- time
    fid.write('    {0:22}:\n'.format('mid-epoch_time'))
    long_name = 'mid-date of each measurement epoch'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
    fid.write('      {0:20}: {1}\n'.format('units', 'decimal-years'))
    fid.write('      {0:20}: {1}\n'.format('comment', '1st column'))
    #-- C10
    fid.write('    {0:22}:\n'.format('C10'))
    long_name = 'C10 coefficient; cosine coefficient for degree 1 and order 0'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '2nd column'))
    #-- C11
    fid.write('    {0:22}:\n'.format('C11'))
    long_name = 'C11 coefficient; cosine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '3rd column'))
    #-- S11
    fid.write('    {0:22}:\n'.format('S11'))
    long_name = 'S11 coefficient; sine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '4th column'))
    #-- eC10
    fid.write('    {0:22}:\n'.format('eC10'))
    long_name = 'eC10 uncertainty; cosine coefficient for degree 1 and order 0'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '5th column'))
    #-- eC11
    fid.write('    {0:22}:\n'.format('eC11'))
    long_name = 'eC11 uncertainty; cosine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '6th column'))
    #-- eS11
    fid.write('    {0:22}:\n'.format('eS11'))
    long_name = 'eS11 uncertainty; sine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '7th column'))
    #-- GRACE month
    fid.write('    {0:22}:\n'.format('month'))
    long_name = 'GRACE month of each measurement epoch'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    description = 'months starting 2002-01-01'
    fid.write('      {0:20}: {1}\n'.format('description', description))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
    fid.write('      {0:20}: {1}\n'.format('units', 'month'))
    fid.write('      {0:20}: {1}\n'.format('comment', '8th column'))
    #-- end of header
    fid.write('\n\n# End of YAML header\n')

#-- PURPOSE: print a file log for the GRACE degree one analysis
def output_log_file(arguments,output_files):
    #-- format: monte_carlo_degree_one_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'monte_carlo_degree_one_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.join(arguments.directory,'geocenter')
    #-- create a unique log and open the log file
    fid = utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    #-- print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(arguments).items()):
        logging.info('{0}: {1}'.format(arg, value))
    #-- print number of monte carlo iterations used in calculation
    logging.info('\n\nNUMBER OF ITERATIONS: {0:d}'.format(arguments.runs))
    #-- print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info('{0}'.format(f))
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE degree one analysis
def output_error_log_file(arguments):
    #-- format: monte_carlo_degree_one_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'monte_carlo_degree_one_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.join(arguments.directory,'geocenter')
    #-- create a unique log and open the log file
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
        description="""Calculates degree 1 errors using GRACE/GRACE-FO
            coefficients of degree 2 and greater, and ocean bottom pressure
            variations from OMCT/MPIOM in a Monte Carlo scheme
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
    #-- GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
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
    #-- number of monte carlo iterations
    parser.add_argument('--runs',
        type=int, default=10000,
        help='Number of Monte Carlo iterations')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    parser.add_argument('--kl','-k',
        type=float, default=0.021,
        help='Degree 1 gravitational Load Love number')
    #-- Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    #-- Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
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
        help='Input/Output data format for delta harmonics file')
    #-- mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    #-- input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    #-- additional error files to be used in the monte carlo run
    parser.add_argument('--error-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', default=[],
        help='Additional error files to use in Monte Carlo analysis')
    #-- run with sea level fingerprints
    parser.add_argument('--fingerprint',
        default=False, action='store_true',
        help='Redistribute land-water flux using sea level fingerprints')
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    #-- land-sea mask for calculating ocean mass and land water flux
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Land-sea mask for calculating ocean mass and land water flux')
    #-- create output plots
    parser.add_argument('--plot','-p',
        default=False, action='store_true',
        help='Create output plots for Monte Carlo iterations')
    #-- Output log file for each job in forms
    #-- monte_carlo_degree_one_run_2002-04-01_PID-00000.log
    #-- monte_carlo_degree_one_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    #-- print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
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
        #-- run monte_carlo_degree_one algorithm with parameters
        output_files = monte_carlo_degree_one(
            args.directory,
            args.center,
            args.release,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            RUNS=args.runs,
            LOVE_NUMBERS=args.love,
            LOVE_K1=args.kl,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            POLE_TIDE=args.pole_tide,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            ERROR_FILES=args.error_file,
            FINGERPRINT=args.fingerprint,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            PLOT=args.plot,
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
