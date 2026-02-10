#!/usr/bin/env python
u"""
kernel_degree_one.py
Written by Tyler Sutterley (06/2024)

Calculates the sensitivity of geocenter calculations for each degree and order
Can be used to estimate the geocenter uncertainties for sets of harmonics

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
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    -k X, --kl X: Degree 1 Gravitational Load Love number
    -F X, --format X: Output data format
        ascii
        netCDF4
        HDF5
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
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: converts a land-sea mask to a series of spherical harmonics
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
    Updated 06/2024: use wrapper to importlib for optional dependencies
    Updated 09/2023: simplify I-matrix and G-matrix calculations
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2023: add options for least-squares solver
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: set default land-sea mask file in arguments
    Updated 05/2022: use argparse descriptions within documentation
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: additional file level attributes in output file
        can use variable loglevels for verbose output
    Written 11/2021
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
import scipy.linalg
import gravity_toolkit as gravtk

# attempt imports
netCDF4 = gravtk.utilities.import_dependency('netCDF4')

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate a geocenter time-series
def kernel_degree_one(base_dir, LMAX, RAD,
    MMAX=None,
    LOVE_NUMBERS=0,
    LOVE_K1=None,
    SOLVER=None,
    FINGERPRINT=False,
    EXPANSION=None,
    LANDMASK=None,
    MODE=0o775):

    # output directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('geocenter')
    # create output directory if non-existent
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # list object of output files for file logs (full path)
    output_files = []

    # read load love numbers
    LOVE = gravtk.load_love_numbers(EXPANSION,
        LOVE_NUMBERS=LOVE_NUMBERS, REFERENCE='CF',
        FORMAT='class')
    # set gravitational load love number to a specific value
    if LOVE_K1:
        LOVE.kl[1] = np.copy(LOVE_K1)
    # minimum spherical harmonic degree
    LMIN = 1
    # maximum spherical harmonic order
    if not MMAX:
        MMAX = np.copy(LMAX)

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    rho_e = factors.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e# Average Radius of the Earth [cm]
    l = factors.l
    # Factor for converting to Mass SH
    dfactor = factors.get('cmwe')

    # Read Smoothed Ocean and Land Functions
    # Open the land-sea NetCDF file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False, varname='LSMASK')
    # degree spacing and grid dimensions
    # will create GRACE spatial fields with same dimensions
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
    # calculate up to degree and order of spherical harmonic expansion for SLF
    PLM, dPLM = gravtk.plm_holmes(EXPANSION, np.cos(th))

    # calculate spherical harmonics of ocean function to degree 1
    # mass is equivalent to 1 cm ocean height change
    # eustatic ratio = -land total/ocean total
    ocean_Ylms = gravtk.gen_stokes(ocean_function, landsea.lon, landsea.lat,
        UNITS=1, LMIN=0, LMAX=1, LOVE=LOVE, PLM=PLM[:2,:2,:])

    # Gaussian Smoothing (Jekeli, 1981)
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    m = np.arange(0,MMAX+1)[:, np.newaxis]
    # Integration factors (solid angle)
    int_fact = np.sin(th)*dphi*dth
    # Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))

    # Legendre polynomials for degree 1
    P10 = np.squeeze(PLM[1,0,:])
    P11 = np.squeeze(PLM[1,1,:])

    # Initializing 3x3 I-Parameter matrix
    IMAT = np.zeros((3,3))
    # Calculating I-Parameter matrix by integrating over latitudes
    # I-Parameter matrix accounts for the fact that the GRACE data only
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

    # output flag for using sea level fingerprints
    slf_str = '_SLF' if FINGERPRINT else ''
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''

    # dictionary with output variables
    output = {}
    # restructured degree and order
    output['l'] = np.zeros((n_harm), dtype=np.int32)
    output['m'] = np.zeros((n_harm), dtype=np.int32)
    output['cs'] = np.zeros((n_harm), dtype=np.int32)
    # 1-based index of geocenter
    output['geocenter'] = 1 + np.arange(3)
    # Legendre polynomials for each degree and order
    plm = np.zeros((n_harm,nlat))
    # cosine and sine factors
    mphi = np.zeros((n_harm,nlon))
    # factor for converting to coefficients of mass
    fit_factor = np.zeros((n_harm))
    # ii is a counter variable for building the column array
    ii = 0
    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # copy dimensions
                output['l'][ii] = l
                output['m'][ii] = m
                output['cs'][ii] = cs
                # legendre polynomials
                plm[ii,:] = np.copy(PLM[l,m,:])
                # degree dependent factor to convert to mass
                fit_factor[ii] = wt[l]*(2.0*l + 1.0)/(1.0 + LOVE.kl[l])
                # cosine and sine factors
                if (csharm == 'clm'):
                    mphi[ii,:] = ccos[m,:]
                elif (csharm == 'slm'):
                    mphi[ii,:] = ssin[m,:]
                # add 1 to counter
                ii += 1

    # degree one covariance matrix for each geocenter component
    output['covariance'] = np.zeros((n_harm, n_harm, 3))
    for i in range(n_harm):
        # setting kern_i equal to 1 for d/o
        kern_i = np.zeros((n_harm,nlat))
        kern_i[i,:] = 1.0*fit_factor[i]
        # calculate land mass maps
        lmass = np.dot(mphi.T,plm*kern_i)
        # Calculating data matrices
        # GRACE Eustatic degree 1 from land variations
        eustatic = gravtk.geocenter()
        # use sea level fingerprints or eustatic from GRACE land components
        if FINGERPRINT:
            # calculate total sea level fingerprint for eustatic component
            # steps to calculate sea level from GRACE land-water change:
            # 1) calculate total land mass
            # NOTE: this is an unscaled GRACE estimate that uses the
            # buffered land function when solving the sea-level equation.
            # possible improvement using scaled estimate with real coastlines
            land_Ylms = gravtk.gen_stokes(lmass*land_function,
                landsea.lon, landsea.lat, UNITS=1, LMIN=0,
                LMAX=EXPANSION, PLM=PLM, LOVE=LOVE)
            # 2) calculate sea level fingerprints of land mass
            # use maximum of 3 iterations for computational efficiency
            sea_level = gravtk.sea_level_equation(land_Ylms.clm, land_Ylms.slm,
                landsea.lon, landsea.lat, land_function, LMAX=EXPANSION,
                LOVE=LOVE_K1, BODY_TIDE_LOVE=0, FLUID_LOVE=0, ITERATIONS=3,
                POLAR=True, PLM=PLM, ASTYPE=np.float64, SCALE=1e-32, FILL_VALUE=0)
            # 3) convert sea level fingerprints into spherical harmonics
            slf_Ylms = gravtk.gen_stokes(sea_level, landsea.lon, landsea.lat,
                UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:], LOVE=LOVE)
            # 4) convert the slf degree 1 harmonics to mass with dfactor
            eustatic.C10 = slf_Ylms.clm[1,0]*dfactor[1]
            eustatic.C11 = slf_Ylms.clm[1,1]*dfactor[1]
            eustatic.S11 = slf_Ylms.slm[1,1]*dfactor[1]
        else:
            # steps to calculate eustatic component from GRACE land-water change:
            # 1) calculate total mass of 1 cm of ocean height (calculated above)
            # 2) calculate total land mass
            # NOTE: possible improvement using the sea-level equation to solve
            # for the spatial pattern of sea level from the land water mass
            land_Ylms = gravtk.gen_stokes(lmass*land_function, landsea.lon,
                landsea.lat, UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:],
                LOVE=LOVE)
            # 3) calculate ratio between the total land mass and the total mass
            # of 1 cm of ocean height (negative as positive land = sea level drop)
            # this converts the total land change to ocean height change
            eustatic_ratio = -land_Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
            # 4) scale degree one coefficients of ocean function with ratio
            # and convert the eustatic degree 1 harmonics to mass with dfactor
            scale_factor = eustatic_ratio*dfactor[1]
            eustatic.C10 = ocean_Ylms.clm[1,0]*scale_factor
            eustatic.C11 = ocean_Ylms.clm[1,1]*scale_factor
            eustatic.S11 = ocean_Ylms.slm[1,1]*scale_factor

        # eustatic coefficients of degree 1
        CMAT = np.array([eustatic.C10, eustatic.C11, eustatic.S11])

        # calculate ocean mass for each degree and order
        for j in range(n_harm):
            # setting kern_j equal to 1 for d/o
            kern_j = np.zeros((n_harm,nlat))
            # skipping C10, C11 and S11 for ocean mass
            if (j >= 3):
                kern_j[j,:] = 1.0*fit_factor[j]
            # calculate ocean mass maps
            rmass = np.dot(mphi.T,plm*kern_j)
            # Allocate for G matrix parameters
            # G matrix calculates the GRACE ocean mass variations
            G = gravtk.geocenter()
            G.C10 = 0.0
            G.C11 = 0.0
            G.S11 = 0.0
            # calculate G matrix parameters through a summation of each latitude
            for n in range(0,nlat):
                # C10, C11, S11
                PC10 = P10[i]*ccos[0,:]
                PC11 = P11[i]*ccos[1,:]
                PS11 = P11[i]*ssin[1,:]
                # summation of integration factors, Legendre polynomials,
                # (convolution of order and harmonics) and the ocean mass
                G.C10 += np.sum(int_fact[n]*PC10*ocean_function[:,n]*rmass[:,n])/(4.0*np.pi)
                G.C11 += np.sum(int_fact[n]*PC11*ocean_function[:,n]*rmass[:,n])/(4.0*np.pi)
                G.S11 += np.sum(int_fact[n]*PS11*ocean_function[:,n]*rmass[:,n])/(4.0*np.pi)

            # G Matrix for time t
            GMAT = np.array([G.C10, G.C11, G.S11])
            # calculate degree 1 solution for harmonic
            # this is mathematically equivalent to an iterative procedure
            # whereby the initial degree one coefficients are used to update
            # the G Matrix until (C10, C11, S11) converge
            if (SOLVER == 'inv'):
                DMAT = np.dot(np.linalg.inv(IMAT), (CMAT-GMAT))
            elif (SOLVER == 'lstsq'):
                DMAT = np.linalg.lstsq(IMAT, (CMAT-GMAT), rcond=-1)[0]
            elif SOLVER in ('gelsd', 'gelsy', 'gelss'):
                DMAT, res, rnk, s = scipy.linalg.lstsq(IMAT, (CMAT-GMAT),
                    lapack_driver=SOLVER)
            # normalize covariances and save to matrix
            output['covariance'][i,j,:] = rho_e*rad_e*DMAT/(4.0*np.pi*dfactor[1])

    # save to file with all descriptor flags
    args = ('COV',LMAX,order_str,gw_str,slf_str,'nc')
    FILE = DIRECTORY.joinpath('{0}_L{1:d}{2}{3}{4}.{5}'.format(*args))
    ncdf_covariance(output, FILENAME=FILE)
    # change the permissions mode
    FILE.chmod(mode=MODE)

    # return the list of output files
    output_files.append(FILE)
    return output_files

# PURPOSE: Write spherical harmonic covariance coefficients to file
def ncdf_covariance(output, **kwargs):
    """
    Writes spherical harmonic covariance coefficients to netCDF4 files

    Arguments
    ---------
    output: dictionary with output variables

    Keyword arguments
    -----------------
    FILENAME: netCDF4 filename
    UNITS: spherical harmonic units
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing netCDF4 file
    DATE: harmonics have date information
    """
    # set default keyword arguments
    kwargs.setdefault('FILENAME',None)
    kwargs.setdefault('UNITS','Geodesy_Normalization')
    kwargs.setdefault('TITLE',None)
    kwargs.setdefault('DATE',True)
    kwargs.setdefault('CLOBBER',True)

    # setting NetCDF clobber attribute
    clobber = 'w' if kwargs['CLOBBER'] else 'a'
    # opening netCDF file for writing
    fileID = netCDF4.Dataset(kwargs['FILENAME'], clobber, format="NETCDF4")

    # Calculating the number of cos and sin harmonics up to LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm,_,n_geo = output['covariance'].shape

    # Defining the netCDF dimensions
    fileID.createDimension('lm', n_harm)
    fileID.createDimension('geocenter', n_geo)

    # defining the netCDF variables
    nc = {}
    # degree and order
    nc['l'] = fileID.createVariable('l', 'i', ('lm',))
    nc['m'] = fileID.createVariable('m', 'i', ('lm',))
    nc['cs'] = fileID.createVariable('cs', 'i', ('lm',))
    nc['geocenter'] = fileID.createVariable('geocenter', 'i',
        ('geocenter',))
    # spherical harmonics
    nc['covariance'] = fileID.createVariable('covariance', 'd',
        ('lm','lm','geocenter'))

    # filling netCDF variables
    for key,val in output.items():
        nc[key][:] = val.copy()

    # Defining attributes for degree and order
    nc['l'].long_name = 'spherical_harmonic_degree'# SH degree long name
    nc['l'].units = 'Wavenumber'# SH degree units
    nc['m'].long_name = 'spherical_harmonic_order'# SH order long name
    nc['m'].units = 'Wavenumber'# SH order units
    nc['cs'].long_name = 'cosine/sine harmonics'
    # Defining attributes for harmonics
    nc['covariance'].long_name = 'spherical_harmonic_covariance'
    nc['covariance'].units = kwargs['UNITS']
    # defining attributes for geocenter
    nc['geocenter'].units = '1'
    nc['geocenter'].long_name = 'Geocenter'
    nc['geocenter'].flag_meanings = '1: C10, 2: C11, 3: S11'
    nc['geocenter'].flag_values = [1,2,3]
    nc['geocenter'].valid_min = 1
    nc['geocenter'].valid_max = 3

    # global variables of NetCDF file
    fileID.creator_name = 'Tyler C. Sutterley and Isabella Velicogna'
    fileID.creator_email = 'tsutterl@uw.edu and isabella@uci.edu'
    fileID.creator_url = 'https://www.ess.uci.edu/~velicogna/index.html'
    fileID.reference = 'https://doi.org/10.3390/rs11182108'
    fileID.creator_type = 'group'
    fileID.creator_institution = ('University of Washington; '
        'University of California, Irvine')
    fileID.history = 'Created at UC Irvine'
    fileID.source = 'derived'
    fileID.title = 'Geocenter covariance'
    fileID.summary = ('Geocenter covariance coefficients.  Geocenter '
        'coefficients represent the largest-scale variability of '
        'hydrologic, cryospheric, and solid Earth processes.')
    project = []
    project.append('NASA Gravity Recovery And Climate Experiment (GRACE)')
    project.append('GRACE Follow-On (GRACE-FO)')
    fileID.project = ', '.join(project)
    fileID.processing_level = 3
    ack = []
    ack.append('GRACE is a joint mission of NASA (USA) and DLR (Germany)')
    ack.append('GRACE-FO is a joint mission of NASA (USA) and GFZ (Germany)')
    fileID.acknowledgement = ', '.join(ack)
    keywords = []
    keywords.append('GRACE')
    keywords.append('GRACE-FO')
    keywords.append('Spherical Harmonic Model')
    keywords.append('Gravitational Field')
    keywords.append('Geopotential')
    keywords.append('Time Variable Gravity')
    keywords.append('Mass Transport')
    keywords.append('Satellite Geodesy')
    fileID.keywords = ', '.join(keywords)
    vocabulary = 'NASA Global Change Master Directory (GCMD) Science Keywords'
    fileID.keywords_vocabulary = vocabulary
    # add software information
    fileID.software_reference = gravtk.version.project_name
    fileID.software_version = gravtk.version.full_version
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output netCDF structure information
    logging.info(kwargs['FILENAME'])
    logging.info(list(fileID.variables.keys()))

    # Closing the netCDF file
    fileID.close()

# PURPOSE: print a file log for the GRACE degree one analysis
def output_log_file(input_arguments, output_files):
    # format: kernel_degree_one_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'kernel_degree_one_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = pathlib.Path(input_arguments.directory).joinpath('geocenter')
    # create a unique log and open the log file
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the GRACE degree one analysis
def output_error_log_file(input_arguments):
    # format: kernel_degree_one_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'kernel_degree_one_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = pathlib.Path(input_arguments.directory).joinpath('geocenter')
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
        description="""Calculates the sensitivity of geocenter calculations
            for each degree and order
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
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
    parser.add_argument('--kl','-k',
        type=float, default=0.021,
        help='Degree 1 gravitational Load Love number')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for ocean models')
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
    # Output log file for each job in forms
    # kernel_degree_one_run_2002-04-01_PID-00000.log
    # kernel_degree_one_failed_run_2002-04-01_PID-00000.log
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
        # run kernel_degree_one algorithm with parameters
        output_files = kernel_degree_one(
            args.directory,
            args.lmax,
            args.radius,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            LOVE_K1=args.kl,
            SOLVER=args.solver,
            FINGERPRINT=args.fingerprint,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
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
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
