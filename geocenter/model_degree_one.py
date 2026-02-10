#!/usr/bin/env python
u"""
model_degree_one.py
Written by Tyler Sutterley (01/2025)

Calculates degree 1 variations using synthetic coefficients of degree 2 and
    greater for testing the reliability of the algorithm

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

INPUTS:
    Input index file for model harmonics

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: Output data directory
    -P X, --file-prefix X: prefix string for output files
    -D, --date: Model harmonics are a time series
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    -k X, --kl X: Degree 1 Gravitational Load Love number
    -F X, --format X: Input/output data format
        ascii
        netCDF4
        HDF5
    --iterative: Iterate degree one solutions
    -s X, --solver X: Least squares solver for degree one solutions
        inv: matrix inversion
        lstsq: least squares solution
        gelsy: complete orthogonal factorization
        gelss: singular value decomposition (SVD)
        gelsd: singular value decomposition (SVD) with divide and conquer method
    --fingerprint: Redistribute land-water flux using sea level fingerprints
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --mask X: Land-sea mask for calculating ocean mass and land water flux
    --log: Output log of files created for each job
    -p, --plot: Create output plots for components and iterations
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://scipy.org
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
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    sea_level_equation.py: pseudo-spectral sea level equation solver
    geocenter.py: converts between spherical harmonics and geocenter variations
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

UPDATE HISTORY:
    Updated 01/2025: fixed deprecated tick label resizing
    Updated 09/2023: simplify I-matrix and G-matrix calculations
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2023: add options for least-squares solver
        place matplotlib import within try/except statement
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: set default land-sea mask file in arguments
    Updated 07/2022: set plot tick formatter to not use offsets
    Updated 05/2022: use argparse descriptions within documentation
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: use gravity_toolkit geocenter class for operations
    Updated 10/2021: using python logging for handling verbose output
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 07/2020: using spatial class for reading and operating spatial data
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2019: added parameter LANDMASK for setting the land-sea mask
    Written 10/2018
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import pathlib
import argparse
import warnings
import traceback
import numpy as np
import scipy.linalg
import gravity_toolkit as gravtk

# attempt imports
try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.offsetbox
    from matplotlib.ticker import MultipleLocator
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.warn("matplotlib not available", ImportWarning)

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: model the seasonal component of an initial degree 1 model
# using preliminary estimates of annual and semi-annual variations from LWM
# as calculated in Chen et al. (1999), doi:10.1029/1998JB900019
# NOTE: this is to get an accurate assessment of the land water mass for the
# eustatic component (not for the ocean component)
def model_seasonal_geocenter(grace_date):
    # Annual amplitudes of (Soil Moisture + Snow) geocenter components (mm)
    AAx = 1.28
    AAy = 0.52
    AAz = 3.30
    # Annual phase of (Soil Moisture + Snow) geocenter components (degrees)
    APx = 44.0
    APy = 182.0
    APz = 43.0
    # Semi-Annual amplitudes of (Soil Moisture + Snow) geocenter components
    SAAx = 0.15
    SAAy = 0.56
    SAAz = 0.50
    # Semi-Annual phase of (Soil Moisture + Snow) geocenter components
    SAPx = 331.0
    SAPy = 312.0
    SAPz = 75.0
    # calculate each geocenter component from the amplitude and phase
    # converting the phase from degrees to radians
    X = AAx*np.sin(2.0*np.pi*grace_date + APx*np.pi/180.0) + \
        SAAx*np.sin(4.0*np.pi*grace_date + SAPx*np.pi/180.0)
    Y = AAy*np.sin(2.0*np.pi*grace_date + APy*np.pi/180.0) + \
        SAAy*np.sin(4.0*np.pi*grace_date + SAPy*np.pi/180.0)
    Z = AAz*np.sin(2.0*np.pi*grace_date + APz*np.pi/180.0) + \
        SAAz*np.sin(4.0*np.pi*grace_date + SAPz*np.pi/180.0)
    DEG1 = gravtk.geocenter(X=X-X.mean(), Y=Y-Y.mean(), Z=Z-Z.mean())
    return DEG1.from_cartesian()

# PURPOSE: calculate a geocenter time-series
def model_degree_one(input_file, LMAX, RAD,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    LOVE_K1=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    DATE=False,
    DATAFORM=None,
    ITERATIVE=False,
    SOLVER=None,
    FINGERPRINT=False,
    EXPANSION=None,
    LANDMASK=None,
    PLOT=False,
    MODE=0o775):

    # create output directory if currently non-existent
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # list object of output files for file logs (full path)
    output_files = []
    # output flag for using sea level fingerprints
    slf_str = 'SLF_' if FINGERPRINT else ''

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE='CF', FORMAT='class')
    # set gravitational load love number to a specific value
    if LOVE_K1:
        LOVE.kl[1] = np.copy(LOVE_K1)
    # maximum spherical harmonic order
    if not MMAX:
        MMAX = np.copy(LMAX)

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    rho_e = factors.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e# Average Radius of the Earth [cm]
    l = factors.l
    # Factor for converting to Mass SH
    dfactor = factors.get('cmwe')

    # Read Smoothed Ocean and Land Functions
    # smoothed functions are from the read_ocean_function.py program
    # Open the land-sea NetCDF file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False, varname='LSMASK')
    # degree spacing and grid dimensions
    # will create spatial fields with same dimensions
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # spatial parameters in radians
    dphi = dlon*np.pi/180.0
    dth = dlat*np.pi/180.0
    # longitude and colatitude in radians
    phi = landsea.lon[np.newaxis,:]*np.pi/180.0
    th = (90.0 - np.squeeze(landsea.lat))*np.pi/180.0
    # create land function
    land_function = np.zeros((nlon, nlat),dtype=np.float64)
    # extract land function from file
    # combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data.T >= 1) & (landsea.data.T <= 3))
    land_function[indx,indy] = 1.0
    # calculate ocean function from land function
    ocean_function = 1.0 - land_function

    # Calculating Legendre Polynomials using Holmes and Featherstone relation
    PLM, dPLM = gravtk.plm_holmes(EXPANSION, np.cos(th))

    # calculate spherical harmonics of ocean function to degree 1
    # mass is equivalent to 1 cm ocean height change
    # eustatic ratio = -land total/ocean total
    ocean_Ylms = gravtk.gen_stokes(ocean_function, landsea.lon, landsea.lat,
        UNITS=1, LMIN=0, LMAX=1, LOVE=LOVE, PLM=PLM[:2,:2,:])

    # Gaussian Smoothing (Jekeli, 1981)
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
    else:
        # else = 1
        wt = np.ones((LMAX+1))

    # input spherical harmonic datafiles
    # input spherical harmonic datafile index
    # correspond file names with GRACE month
    # this allows additional months to be in the index
    data_Ylms = gravtk.harmonics().from_index(input_file, format=DATAFORM)
    # number of files within the index
    n_files = data_Ylms.shape[-1]
    # truncate to degree and order
    data_Ylms.truncate(lmax=LMAX,mmax=MMAX)
    # save original geocenter for file in mm w.e.
    DATA = gravtk.geocenter().from_harmonics(data_Ylms).scale(10.0*(rho_e*rad_e))
    DATA.mean(apply=True)
    # extract date variables
    if DATE:
        tdec = np.copy(data_Ylms.time)
        mon = np.copy(data_Ylms.month)
    # filter data coefficients
    if DESTRIPE:
        data_Ylms = data_Ylms.destripe()
        ds_str = '_FL'
    else:
        ds_str = ''

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    m = data_Ylms.m[:, np.newaxis]
    # Integration factors (solid angle)
    int_fact = np.sin(th)*dphi*dth
    # Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))

    # Legendre polynomials for degree 1
    P10 = np.squeeze(PLM[1,0,:])
    P11 = np.squeeze(PLM[1,1,:])
    # PLM for spherical harmonic degrees 2+ up to LMAX
    # converted into mass and smoothed if specified
    plmout = np.zeros((LMAX+1, MMAX+1, nlat))
    for l in range(1,LMAX+1):
        m = np.arange(0,np.min([l,MMAX])+1)
        # convert to smoothed coefficients of mass
        # Convolving plms with degree dependent factor and smoothing
        plmout[l,m,:] = PLM[l,m,:]*dfactor[l]*wt[l]

    # Initializing 3x3 I-Parameter matrix
    IMAT = np.zeros((3,3))
    # Calculating I-Parameter matrix by integrating over latitudes
    # I-Parameter matrix accounts for the fact that the harmonic data only
    # includes spherical harmonic degrees greater than or equal to 2
    for i in range(0,nlat):
        # C10, C11, S11
        PC10 = P10[i]*ccos[0,:]
        PC11 = P11[i]*ccos[1,:]
        PS11 = P11[i]*ssin[1,:]
        # C10: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,0] += np.sum(int_fact[i]*PC10*ocean_function[:,i]*PC10)/(4.0*np.pi)
        IMAT[1,0] += np.sum(int_fact[i]*PC10*ocean_function[:,i]*PC11)/(4.0*np.pi)
        IMAT[2,0] += np.sum(int_fact[i]*PC10*ocean_function[:,i]*PS11)/(4.0*np.pi)
        # C11: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,1] += np.sum(int_fact[i]*PC11*ocean_function[:,i]*PC10)/(4.0*np.pi)
        IMAT[1,1] += np.sum(int_fact[i]*PC11*ocean_function[:,i]*PC11)/(4.0*np.pi)
        IMAT[2,1] += np.sum(int_fact[i]*PC11*ocean_function[:,i]*PS11)/(4.0*np.pi)
        # S11: C10, C11, S11 (see equations 12 and 13 of Swenson et al., 2008)
        IMAT[0,2] += np.sum(int_fact[i]*PS11*ocean_function[:,i]*PC10)/(4.0*np.pi)
        IMAT[1,2] += np.sum(int_fact[i]*PS11*ocean_function[:,i]*PC11)/(4.0*np.pi)
        IMAT[2,2] += np.sum(int_fact[i]*PS11*ocean_function[:,i]*PS11)/(4.0*np.pi)

    # get seasonal variations of an initial geocenter correction
    # for use in the land water mass calculation
    seasonal_geocenter = model_seasonal_geocenter(tdec)

    # iterate solutions: if not single iteration
    n_iter = 0
    eps = np.inf
    eps_max = 1e-6
    if ITERATIVE:
        iter_str = 'iter_'
        max_iter = 15
    else:
        iter_str = ''
        max_iter = 1

    # Calculating data matrices
    # Allocate for G matrix parameters
    # G matrix calculates the ocean water mass variations
    G = gravtk.geocenter()
    G.C10 = np.zeros((n_files))
    G.C11 = np.zeros((n_files))
    G.S11 = np.zeros((n_files))
    # DMAT is the degree one matrix ((C10,C11,S11) x Time) in terms of mass
    DMAT = np.zeros((3,n_files))
    # degree 1 iterations
    iteration = gravtk.geocenter()
    iteration.C10 = np.zeros((n_files, max_iter))
    iteration.C11 = np.zeros((n_files, max_iter))
    iteration.S11 = np.zeros((n_files, max_iter))
    # calculate non-iterated terms for each file (G-matrix parameters)
    for t in range(n_files):
        # calculate geocenter component of ocean mass
        # allocate for product of grace and legendre polynomials
        pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
        psin = np.zeros((MMAX+1, nlat))#-[m,lat]
        # Summing product of plms and c/slms over all SH degrees >= 2
        for i in range(0, nlat):
            l = np.arange(2,LMAX+1)
            pcos[:,i] = np.sum(plmout[l,:,i]*data_Ylms.clm[l,:,t], axis=0)
            psin[:,i] = np.sum(plmout[l,:,i]*data_Ylms.slm[l,:,t], axis=0)
        # Multiplying by c/s(phi#m) to get surface density in cmwe (lon,lat)
        # ccos/ssin are mXphi, pcos/psin are mXtheta: resultant matrices are phiXtheta
        # The summation over spherical harmonic order is in this multiplication
        rmass = np.dot(np.transpose(ccos),pcos) + np.dot(np.transpose(ssin),psin)
        # calculate G matrix parameters through a summation of each latitude
        for i in range(0,nlat):
            # C10, C11, S11
            PC10 = P10[i]*ccos[0,:]
            PC11 = P11[i]*ccos[1,:]
            PS11 = P11[i]*ssin[1,:]
            # summation of integration factors, Legendre polynomials,
            # (convolution of order and harmonics) and the ocean mass at t
            G.C10[t] += np.sum(int_fact[i]*PC10*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)
            G.C11[t] += np.sum(int_fact[i]*PC11*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)
            G.S11[t] += np.sum(int_fact[i]*PS11*ocean_function[:,i]*rmass[:,i])/(4.0*np.pi)

    # calculate degree one solution for each iteration (or single if not)
    while (eps > eps_max) and (n_iter < max_iter):
        # for each file
        for t in range(n_files):
            # calculate eustatic component (can iterate)
            if (n_iter == 0):
                # for first iteration (will be only iteration if not ITERATIVE):
                # seasonal component of geocenter variation for land water
                data_Ylms.clm[1,0,t] = seasonal_geocenter.C10[t]
                data_Ylms.clm[1,1,t] = seasonal_geocenter.C11[t]
                data_Ylms.slm[1,1,t] = seasonal_geocenter.S11[t]
            else:
                # for all others: use previous iteration of inversion
                # for each of the geocenter solutions (C10, C11, S11)
                data_Ylms.clm[1,0,t] = iteration.C10[t,n_iter-1]
                data_Ylms.clm[1,1,t] = iteration.C11[t,n_iter-1]
                data_Ylms.slm[1,1,t] = iteration.S11[t,n_iter-1]

            # allocate for product of grace and legendre polynomials
            pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
            psin = np.zeros((MMAX+1, nlat))#-[m,lat]
            # Summing product of plms and c/slms over all SH degrees
            for i in range(0, nlat):
                l = np.arange(1,LMAX+1)
                pcos[:,i] = np.sum(plmout[l,:,i]*data_Ylms.clm[l,:,t], axis=0)
                psin[:,i] = np.sum(plmout[l,:,i]*data_Ylms.slm[l,:,t], axis=0)

            # Multiplying by c/s(phi#m) to get surface density in cm w.e. (lonxlat)
            # this will be a spatial field similar to outputs from stokes_combine.py
            # ccos/ssin are mXphi, pcos/psin are mXtheta: resultant matrices are phiXtheta
            # The summation over spherical harmonic order is in this multiplication
            lmass = np.dot(np.transpose(ccos),pcos) + np.dot(np.transpose(ssin),psin)

            # use sea level fingerprints or eustatic from land components
            if FINGERPRINT:
                # calculate total sea level fingerprint for eustatic component
                # steps to calculate sea level from land-water change:
                # 1) calculate total land mass at time t (data*land function)
                # NOTE: this is an unscaled estimate that uses the
                # buffered land function when solving the sea-level equation.
                # possible improvement using scaled estimate with real coastlines
                land_Ylms = gravtk.gen_stokes(land_function*lmass,
                    landsea.lon, landsea.lat, UNITS=1, LMIN=0, LMAX=EXPANSION,
                    PLM=PLM, LOVE=LOVE)
                # 2) calculate sea level fingerprints of land mass at time t
                # use maximum of 3 iterations for computational efficiency
                sea_level = gravtk.sea_level_equation(land_Ylms.clm, land_Ylms.slm,
                    landsea.lon, landsea.lat, land_function, LMAX=EXPANSION,
                    LOVE=LOVE, BODY_TIDE_LOVE=0, FLUID_LOVE=0, ITERATIONS=3,
                    POLAR=True, PLM=PLM, FILL_VALUE=0)
                # 3) convert sea level fingerprints into spherical harmonics
                slf_Ylms = gravtk.gen_stokes(sea_level, landsea.lon, landsea.lat,
                    UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:], LOVE=LOVE)
                # 4) convert the slf degree 1 harmonics to mass with dfactor
                eustatic = gravtk.geocenter().from_harmonics(slf_Ylms).scale(dfactor[1])
            else:
                # steps to calculate eustatic component from land-water change:
                # 1) calculate total mass of 1 cm of ocean height (calculated above)
                # 2) calculate total land mass at time t (data*land function)
                # NOTE: possible improvement using the sea-level equation to solve
                # for the spatial pattern of sea level from the land water mass
                land_Ylms = gravtk.gen_stokes(land_function*lmass,
                    landsea.lon, landsea.lat, UNITS=1, LMIN=0, LMAX=1,
                    PLM=PLM[:2,:2,:], LOVE=LOVE)
                # 3) calculate ratio between the total land mass and the total mass
                # of 1 cm of ocean height (negative as positive land = sea level drop)
                # this converts the total land change to ocean height change
                eustatic_ratio = -land_Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
                # 4) scale degree one coefficients of ocean function with ratio
                # and convert the eustatic degree 1 harmonics to mass with dfactor
                scale_factor = eustatic_ratio*dfactor[1]
                eustatic = gravtk.geocenter().from_harmonics(ocean_Ylms).scale(scale_factor)

            # eustatic coefficients of degree 1
            CMAT = np.array([eustatic.C10,eustatic.C11,eustatic.S11])
            # G Matrix for time t
            GMAT = np.array([G.C10[t], G.C11[t], G.S11[t]])
            # calculate degree 1 solution for iteration
            # this is mathematically equivalent to an iterative procedure
            # whereby the initial degree one coefficients are used to update
            # the G Matrix until (C10, C11, S11) converge
            # calculates min(eustatic from land - measured ocean)
            if (SOLVER == 'inv'):
                DMAT[:,t] = np.dot(np.linalg.inv(IMAT), (CMAT-GMAT))
            elif (SOLVER == 'lstsq'):
                DMAT[:,t] = np.linalg.lstsq(IMAT, (CMAT-GMAT), rcond=-1)[0]
            elif SOLVER in ('gelsd', 'gelsy', 'gelss'):
                DMAT[:,t], res, rnk, s = scipy.linalg.lstsq(IMAT, (CMAT-GMAT),
                    lapack_driver=SOLVER)
            # save geocenter for iteration and time t
            iteration.C10[t,n_iter] = DMAT[0,t]/dfactor[1]
            iteration.C11[t,n_iter] = DMAT[1,t]/dfactor[1]
            iteration.S11[t,n_iter] = DMAT[2,t]/dfactor[1]
        # remove mean of each solution for iteration
        iteration.C10[:,n_iter] -= iteration.C10[:,n_iter].mean()
        iteration.C11[:,n_iter] -= iteration.C11[:,n_iter].mean()
        iteration.S11[:,n_iter] -= iteration.S11[:,n_iter].mean()
        # calculate difference between original geocenter coefficients and the
        # calculated coefficients for each of the geocenter solutions
        sigma_C10 = np.sum((data_Ylms.clm[1,0,:] - iteration.C10[:,n_iter])**2)
        sigma_C11 = np.sum((data_Ylms.clm[1,1,:] - iteration.C11[:,n_iter])**2)
        sigma_S11 = np.sum((data_Ylms.slm[1,1,:] - iteration.S11[:,n_iter])**2)
        power = data_Ylms.clm[1,0,t]**2 + data_Ylms.clm[1,1,t]**2 + data_Ylms.slm[1,1,t]**2
        eps = np.sqrt(sigma_C10 + sigma_C11 + sigma_S11)/np.sqrt(np.sum(power))
        # add 1 to n_iter counter
        n_iter += 1

    # Convert inverted solutions into fully normalized spherical harmonics
    # for each of the geocenter solutions (C10, C11, S11)
    # for the iterative case this will be the final iteration
    DEG1 = gravtk.geocenter()
    DEG1.C10,DEG1.C11,DEG1.S11 = DMAT/dfactor[1]
    DEG1.mean(apply=True)
    # calculate harmonics in mm w.e.
    mmwe = DEG1.scale(10.0*(rho_e*rad_e))

    # output degree 1 coefficients
    # 1: mm water equivalent of recovered
    # 2: mm water equivalent of actual
    args = (FILE_PREFIX,slf_str,iter_str,ds_str)
    FILE1 = OUTPUT_DIRECTORY.joinpath('{0}{1}{2}mmwe{3}.txt'.format(*args))
    fid1 = FILE1.open(mode='w', encoding='utf8')
    output_files.append(FILE1)
    FILE2 = OUTPUT_DIRECTORY.joinpath('{0}Actual_mmwe.txt'.format(FILE_PREFIX))
    fid2 = FILE2.open(mode='w', encoding='utf8')
    output_files.append(FILE2)
    # print Swenson file headers
    print("  Degree 1 coefficients, mm equivalent water thickness", file=fid1)
    print("  Degree 1 coefficients, mm equivalent water thickness", file=fid2)
    print("  format='(4f9.2,i9)'", file=fid1)
    print("  format='(4f9.2,i9)'", file=fid2)
    args = ['Time','C10','C11','S11','Month']
    print(''.join('{:>9}'.format(s) for s in args), file=fid1)
    print(''.join('{:>9}'.format(s) for s in args), file=fid2)
    # for each file
    for t in range(n_files):
        # output geocenter coefficients to file
        print('{0:9.2f}{1:9.2f}{2:9.2f}{3:9.2f}      {4:03d}'.format(tdec[t],
            mmwe.C10[t], mmwe.C11[t], mmwe.S11[t], mon[t]), file=fid1)
        print('{0:9.2f}{1:9.2f}{2:9.2f}{3:9.2f}      {4:03d}'.format(tdec[t],
            DATA.C10[t], DATA.C11[t], DATA.S11[t], mon[t]), file=fid2)
    # close the output files
    fid1.close()
    fid2.close()
    # set the permissions mode of the output files
    FILE1.chmod(mode=MODE)
    FILE2.chmod(mode=MODE)

    # create plot of recovered versus actual geocenter
    if PLOT:
        # 3 row plot (C10, C11 and S11)
        ax = {}
        fig, (ax[0], ax[1], ax[2]) = plt.subplots(num=1, nrows=3, sharex=True,
            sharey=True, figsize=(6,9))
        # plot original geocenter
        ax[0].plot(tdec, DATA.C10, 'b', lw=2)
        ax[1].plot(tdec, DATA.C11, 'b', lw=2)
        ax[2].plot(tdec, DATA.S11, 'b', lw=2)
        # plot resultant geocenter
        ax[0].plot(tdec, mmwe.C10, color='r', lw=2)
        ax[1].plot(tdec, mmwe.C11, color='r', lw=2)
        ax[2].plot(tdec, mmwe.S11, color='r', lw=2)
        # labels and set limits to Swenson range
        ax[0].set_ylabel('[mm]', fontsize=14)
        ax[1].set_ylabel('[mm]', fontsize=14)
        ax[2].set_ylabel('[mm]', fontsize=14)
        ax[2].set_xlabel('Time [Yr]', fontsize=14)
        #ax[2].set_xlim(2003,2007)
        #ax[2].set_ylim(-6,6)
        #ax[2].xaxis.set_ticks(np.arange(2003,2008,1))
        #ax[2].xaxis.set_minor_locator(MultipleLocator(0.25))
        #ax[2].yaxis.set_ticks(np.arange(-6,8,2))
        ax[2].xaxis.get_major_formatter().set_useOffset(False)
        # add axis labels and adjust font sizes for axis ticks
        for i,lbl in enumerate(['C10','C11','S11']):
            # axis label
            artist = matplotlib.offsetbox.AnchoredText(lbl, pad=0.0,
                frameon=False, loc=2, prop=dict(size=16,weight='bold'))
            ax[i].add_artist(artist)
            # axes tick adjustments
            ax[i].tick_params(axis='both', which='both',
                labelsize=14, direction='in')
        # adjust locations of subplots and save to file
        fig.subplots_adjust(left=0.1,right=0.96,bottom=0.06,top=0.98,hspace=0.1)
        args = (FILE_PREFIX,slf_str,iter_str,ds_str)
        FILE = '{0}{1}{2}Comparison{3}.pdf'.format(*args)
        PLOT1 = OUTPUT_DIRECTORY.joinpath(FILE)
        plt.savefig(PLOT1, format='pdf')
        output_files.append(PLOT1)
        plt.clf()
        # set the permissions mode of the output files
        PLOT1.chmod(mode=MODE)

    # if ITERATIVE: create plot showing iteration solutions
    if PLOT and ITERATIVE:
        # 3 row plot (C10, C11 and S11)
        ax = {}
        fig, (ax[0],ax[1],ax[2]) = plt.subplots(num=1, nrows=3, sharex=True,
            figsize=(6,9))
        # show solutions for each iteration
        plot_colors = iter(cm.rainbow(np.linspace(0,1,n_iter)))
        for j in range(n_iter):
            color_j = next(plot_colors)
            # C10, C11 and S11
            ax[0].plot(mon,10.*(iteration.C10[:,j]*dfactor[1]),color=color_j)
            ax[1].plot(mon,10.*(iteration.C11[:,j]*dfactor[1]),color=color_j)
            ax[2].plot(mon,10.*(iteration.S11[:,j]*dfactor[1]),color=color_j)
        # labels and set limits
        ax[0].set_ylabel('mm', fontsize=14)
        ax[1].set_ylabel('mm', fontsize=14)
        ax[2].set_ylabel('mm', fontsize=14)
        ax[2].set_xlabel('Grace Month', fontsize=14)
        ax[2].set_xlim(np.floor(mon[0]/10.)*10.,np.ceil(mon[-1]/10.)*10.)
        ax[2].xaxis.set_minor_locator(MultipleLocator(5))
        ax[2].xaxis.get_major_formatter().set_useOffset(False)
        # add axis labels and adjust font sizes for axis ticks
        fig_labels = ['C10','C11','S11']
        for i in range(3):
            # axis label
            artist = matplotlib.offsetbox.AnchoredText(fig_labels[i], pad=0.,
                frameon=False, loc=2, prop=dict(size=16,weight='bold'))
            ax[i].add_artist(artist)
            # axes tick adjustments
            ax[i].tick_params(axis='both', which='both',
                labelsize=14, direction='in')
        # adjust locations of subplots and save to file
        fig.subplots_adjust(left=0.12,right=0.94,bottom=0.06,top=0.98,hspace=0.1)
        args = (FILE_PREFIX,slf_str,ds_str)
        FILE = '{0}{1}Geocenter_Iterative{2}.pdf'.format(*args)
        PLOT2 = OUTPUT_DIRECTORY.joinpath(FILE)
        plt.savefig(PLOT2, format='pdf')
        output_files.append(PLOT2)
        plt.clf()
        # set the permissions mode of the output files
        PLOT2.chmod(mode=MODE)

    # return the list of output files and the number of iterations
    return (output_files, n_iter)

# PURPOSE: print a file log for the model degree one analysis
def output_log_file(input_arguments, output_files, n_iter):
    # format: model_degree_one_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'model_degree_one_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    # create a unique log and open the log file
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print number of iterations used in calculation
    if arguments.iterative:
        logging.info('\n\nNUMBER OF ITERATIONS: {0:d}'.format(n_iter))
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the model degree one analysis
def output_error_log_file(input_arguments):
    # format: model_degree_one_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'model_degree_one_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    # create a unique log and open the log file
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates degree 1 variations using synthetic
            coefficients of degree 2 and greater for testing the reliability
            of the algorithm
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='Input index file with spherical harmonic data files')
    # output working data directory
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for output files')
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Model harmonics are a time series')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    parser.add_argument('--kl','-k',
        type=float, default=0.021,
        help='Degree 1 gravitational Load Love number')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for ocean models')
    # run with iterative scheme
    parser.add_argument('--iterative',
        default=False, action='store_true',
        help='Iterate degree one solutions')
    # least squares solver
    choices = ('inv','lstsq','gelsd', 'gelsy', 'gelss')
    parser.add_argument('--solver','-s',
        type=str, default='lstsq', choices=choices,
        help='Least squares solver for degree one solutions')
    # run with sea level fingerprints
    parser.add_argument('--fingerprint',
        default=False, action='store_true',
        help='Redistribute land-water flux using sea level fingerprints')
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # land-sea mask for calculating ocean mass and land water flux
    land_mask_file = gravtk.utilities.get_data_path(['data','land_fcn_300km.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path,
        default=land_mask_file,
        help='Land-sea mask for calculating ocean mass and land water flux')
    # create output plots
    parser.add_argument('--plot','-p',
        default=False, action='store_true',
        help='Create output plots for components and iterations')
    # Output log file for each job in forms
    # model_degree_one_run_2002-04-01_PID-00000.log
    # model_degree_one_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
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
        # run model_degree_one algorithm with parameters
        output_files,n_iter = model_degree_one(
            args.infile,
            args.lmax,
            args.radius,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            LOVE_K1=args.kl,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            DATE=args.date,
            DATAFORM=args.format,
            ITERATIVE=args.iterative,
            SOLVER=args.solver,
            FINGERPRINT=args.fingerprint,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            PLOT=args.plot,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files,n_iter)

# run main program
if __name__ == '__main__':
    main()
