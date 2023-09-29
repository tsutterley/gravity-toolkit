#!/usr/bin/env python
u"""
calc_degree_one.py
Written by Tyler Sutterley (09/2023)

Calculates degree 1 variations using GRACE coefficients of degree 2 and greater,
    and ocean bottom pressure variations from ECCO and OMCT/MPIOM

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
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    -k X, --kl X: Degree 1 Gravitational Load Love number
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
    --slr-c40 X: Replace C40 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --slr-c50 X: Replace C50 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --ocean-model X: Ocean model to use
    -F X, --format X: Input data format for ocean models
        ascii
        netCDF4
        HDF5
    --ocean-file X: Index file for ocean model harmonics
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
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
    -p, --plot: Create output plots for components and iterations
    -C, --copy: Copy output files for distribution and archival
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    scipy: Scientific Tools for Python
        https://scipy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces low-degree harmonics with SLR values (if specified)
    time.py: utilities for calculating time operations
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
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
    Updated 09/2023: output comprehensive netCDF4 files with all components
        simplify I-matrix and G-matrix calculations
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2023: add options for least-squares solver
    Updated 03/2023: place matplotlib import within try/except statement
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 08/2022: set default land-sea mask file in arguments
    Updated 07/2022: set plot tick formatter to not use offsets
    Updated 05/2022: use argparse descriptions within documentation
        use GIA reference and citation output from GIA read program
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: add GSFC low-degree harmonics
        use gravity_toolkit geocenter class for operations
    Updated 10/2021: using python logging for handling verbose output
    Updated 08/2021: updated output filename of netCDF4 iterations file
    Updated 07/2021: fix YAML headers for S11 description
        simplified file imports using wrappers in harmonics
        remove choices for argparse processing centers
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: reference for GFZ oblateness and figure axis solutions
        define int/float precision to prevent deprecation warning
    Updated 04/2021: can replace figure axis and azimuthal dependence with SLR
    Updated 01/2021: harmonics object output from gen_stokes.py
    Updated 12/2020: added more love number options and from gfc for mean files
        updated model parameter to be a string to allow more ECCO model options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: added parameter for make degree of SLF expansion optional
        pass Legendre polynomials up to SLF expansion to functions
        use utilities to define path to load love numbers file
    Updated 07/2020: using spatial class for reading and operating spatial data
    Updated 06/2020: only output public data formats with and without dealiasing
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function.  include output in public format
        using Holmes and Featherstone relation for calculating Legendre polynomials
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: output all iterations to a single netCDF4 file
        can replace C30 with coefficients from satellite laser ranging (SLR)
    Updated 06/2019: added parameter LANDMASK for setting the land-sea mask
    Updated 05/2019: use exact date for calculating GIA drift
    Updated 02/2019: saving metadata to output figure
    Updated 11/2018: add flag if using a GIA model other than AW13 ICE5G
    Updated 10/2018: added sea level equation solver for calculating SLF directly
        read spherical harmonics from ECCO that take into account bathymetry
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 07/2018: can create components plot with OMCT component (no ECCO)
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: use an updated model of ECCO ocean bottom pressure (kf080i)
    Updated 03/2018: added extrapolation of load love numbers if LMAX > 696
    Updated 04/2017: added --mode command line option to set the permissions level
    Updated 03/2017: print number of iterations used to output log file if ITERATIVE
    Updated 06/2016: adjustments to plot similar to Swenson et al. (2008) Figure 1
        edited iteration option with while loop and eps
    Updated 05/2016: using __future__ print function. output with format.
        minor clean up of directories. formatting on mm.w.e.
    Updated 02/2016: compute I-matrix directly rather than for each parameter,
        use getopt parameters to set number of PROCESSES to run in parallel,
            whether or not to output a log file, added help module
    Updated 10/2015: using preliminary geocenter components from Chen (1999) for LWM
    Updated 09/2015: added sea level fingerprint option
    Updated 08/2015: changed sys.exit to a raise exception instance
        Added pole tide and GAE/GAF/GAG correction parameters
    Updated 05/2015: added plot and parameter to recreate Figure 1 of Swenson (2008)
    Updated 03/2015: updated eustatic ratio and updated comments
        added multiprocessing error handling with traceback
    Updated 12/2014: added output in mm.w.e.
    Updated 11/2014: complete code update
        added main definition and multiprocessing
        added updated destriping algorithm
    Forked 06/2013 from calc_deg_one.pro
    Written 09/2012
"""
from __future__ import print_function

import sys
import os
import re
import time
import copy
import shutil
import logging
import pathlib
import argparse
import warnings
import traceback
import collections
import numpy as np
import scipy.linalg
import gravity_toolkit as gravtk

# attempt imports
try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.offsetbox
    import matplotlib.ticker as ticker
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.warn("matplotlib not available", ImportWarning)
try:
    import netCDF4
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.warn("netCDF4 not available", ImportWarning)

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: import GRACE/GRACE-FO GSM files for a given months range
def load_grace_GSM(base_dir, PROC, DREL, START, END, MISSING, LMAX,
    MMAX=None, SLR_C20=None, SLR_21=None, SLR_22=None, SLR_C30=None,
    SLR_C40=None, SLR_C50=None, POLE_TIDE=False):
    # GRACE/GRACE-FO dataset
    DSET = 'GSM'
    # do not import degree 1 coefficients for the GRACE GSM solution
    # 0: No degree 1
    DEG1 = 0
    # reading GRACE/GRACE-FO GSM solutions for input date range
    # replacing low-degree harmonics with SLR values if specified
    # correcting for Pole-Tide if specified
    # atmospheric jumps will be corrected externally if specified
    grace_Ylms = gravtk.grace_input_months(base_dir,
        PROC, DREL, DSET, LMAX, START, END, MISSING, SLR_C20, DEG1,
        MMAX=MMAX, SLR_21=SLR_21, SLR_22=SLR_22, SLR_C30=SLR_C30,
        SLR_C40=SLR_C40, SLR_C50=SLR_C50, POLE_TIDE=POLE_TIDE,
        ATM=False, MODEL_DEG1=False)
    # returning input variables as a harmonics object
    return gravtk.harmonics().from_dict(grace_Ylms)

# PURPOSE: import GRACE/GRACE-FO dealiasing files for a given months range
def load_AOD(base_dir, PROC, DREL, DSET, START, END, MISSING, LMAX):
    # do not replace low degree harmonics for AOD solutions
    SLR_C20 = 'N'
    # do not replace degree 1 coefficients for the GRACE AOD solution
    # 0: No degree 1 replacement
    DEG1 = 0
    # reading GRACE/GRACE-FO AOD solutions for input date range
    grace_Ylms = gravtk.grace_input_months(base_dir,
        PROC, DREL, DSET, LMAX, START, END, MISSING, SLR_C20, DEG1,
        POLE_TIDE=False, ATM=False)
    # returning input variables as a harmonics object
    return gravtk.harmonics().from_dict(grace_Ylms)

# PURPOSE: model the seasonal component of an initial degree 1 model
# using preliminary estimates of annual and semi-annual variations from LWM
# as calculated in Chen et al. (1999), doi:10.1029/1998JB900019
# NOTE: this is to get an accurate assessment of the land water mass for the
# eustatic component (not for the ocean component from GRACE)
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
def calc_degree_one(base_dir, PROC, DREL, MODEL, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
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
    SLR_C40=None,
    SLR_C50=None,
    DATAFORM=None,
    MEAN_FILE=None,
    MEANFORM=None,
    MODEL_INDEX=None,
    ITERATIVE=False,
    SOLVER=None,
    FINGERPRINT=False,
    EXPANSION=None,
    LANDMASK=None,
    PLOT=False,
    COPY=False,
    MODE=0o775):

    # output directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('geocenter')
    # create output directory if non-existent
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output attributes for geocenter netCDF4 files
    attributes = collections.OrderedDict()
    attributes['generating_institute'] = PROC
    attributes['product_release'] = DREL
    attributes['product_name'] = 'GSM'
    attributes['product_type'] = 'gravity_field'
    MISSION = dict(RL05='GRACE', RL06='GRACE/GRACE-FO')
    attributes['title'] = f'{MISSION[DREL]} Geocenter Coefficients'
    attributes['solver'] = SOLVER

    # list object of output files for file logs (full path)
    output_files = []

    # output flag and attributes for using iterative solution
    if ITERATIVE:
        iter_str = '_iter'
        max_iter = 15
        attributes['solution_type'] = 'iterative'
    else:
        iter_str = ''
        max_iter = 1
        attributes['solution_type'] = 'single'
    # output flag and attributes for using sea level fingerprints
    if FINGERPRINT:
        slf_str = '_SLF'
        attributes['eustatic_sea_level'] = 'self_attraction_and_loading'
    else:
        slf_str = ''
        attributes['eustatic_sea_level'] = 'uniform_redistribution'

    # output flag for low-degree harmonic replacements
    if SLR_21 in ('CSR','GFZ','GSFC'):
        C21_str = f'_w{SLR_21}_21'
    else:
        C21_str = ''
    if SLR_22 in ('CSR','GSFC'):
        C22_str = f'_w{SLR_22}_22'
    else:
        C22_str = ''
    if SLR_C30 in ('GSFC',):
        # C30 replacement now default for all solutions
        C30_str = ''
    elif SLR_C30 in ('CSR','GFZ','LARES'):
        C30_str = f'_w{SLR_C30}_C30'
    else:
        C30_str = ''
    if SLR_C40 in ('CSR','GSFC','LARES'):
        C40_str = f'_w{SLR_C40}_C40'
    else:
        C40_str = ''
    if SLR_C50 in ('CSR','GSFC','LARES'):
        C50_str = f'_w{SLR_C50}_C50'
    else:
        C50_str = ''
    # combine satellite laser ranging flags
    slr_str = ''.join([C21_str,C22_str,C30_str,C40_str,C50_str])

    # read load love numbers
    LOVE = gravtk.load_love_numbers(EXPANSION,
        LOVE_NUMBERS=LOVE_NUMBERS, REFERENCE='CF',
        FORMAT='class')
    # add attributes for earth model and love numbers
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # set gravitational load love number to a specific value
    if LOVE_K1:
        LOVE.kl[1] = np.copy(LOVE_K1)

    # maximum spherical harmonic order
    if not MMAX:
        MMAX = np.copy(LMAX)
    # add attributes for LMAX and MMAX
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    rho_e = factors.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e# Average Radius of the Earth [cm]
    l = factors.l
    # Factor for converting to Mass SH
    dfactor = factors.cmwe
    # add attributes for earth parameters
    attributes['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['earth_density'] = f'{factors.rho_e:0.3f} g/cm'
    attributes['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'

    # Read Smoothed Ocean and Land Functions
    # Open the land-sea NetCDF file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK,
        date=False, varname='LSMASK')
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
    ocean_Ylms = gravtk.gen_stokes(ocean_function,
        landsea.lon, landsea.lat, UNITS=1, LMIN=0, LMAX=1,
        LOVE=LOVE, PLM=PLM[:2,:2,:])

    # Gaussian Smoothing (Jekeli, 1981)
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        attributes['smoothing_radius'] = f'{RAD:0.0f} km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))

    # load GRACE/GRACE-FO data
    GSM_Ylms = load_grace_GSM(base_dir, PROC, DREL, START, END, MISSING, LMAX,
        MMAX=MMAX, SLR_C20=SLR_C20, SLR_21=SLR_21, SLR_22=SLR_22,
        SLR_C30=SLR_C30, SLR_C40=SLR_C40, SLR_C50=SLR_C50, POLE_TIDE=POLE_TIDE)
    GAD_Ylms = load_AOD(base_dir, PROC, DREL, 'GAD', START, END, MISSING, LMAX)
    GAC_Ylms = load_AOD(base_dir, PROC, DREL, 'GAC', START, END, MISSING, LMAX)
    # add attributes for input GRACE/GRACE-FO spherical harmonics
    for att_name, att_val in GSM_Ylms.attributes['ROOT'].items():
        attributes[att_name] = att_val
    # use a mean file for the static field to remove
    if MEAN_FILE:
        # read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=MEANFORM, date=False)
        # remove the input mean
        GSM_Ylms.subtract(mean_Ylms)
        attributes['lineage'].append(MEAN_FILE.name)
    else:
        GSM_Ylms.mean(apply=True)
    # remove the mean from the GRACE/GRACE-FO dealiasing data
    GAD_Ylms.mean(apply=True)
    GAC_Ylms.mean(apply=True)
    # convert GAC to geocenter object
    GAC = gravtk.geocenter().from_harmonics(GAC_Ylms)
    # filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        # destriping GRACE GSM and GAD coefficients
        ds_str = '_FL'
        GSM_Ylms = GSM_Ylms.destripe()
        GAD_Ylms = GAD_Ylms.destripe()
        attributes['filtering'] = 'Destriped'
    else:
        # using standard GRACE GSM harmonics
        ds_str = ''
    # GRACE dates
    tdec = np.copy(GSM_Ylms.time)
    months = np.copy(GSM_Ylms.month)
    # number of months considered
    n_files = len(GSM_Ylms)

    # input GIA spherical harmonic datafiles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
    if GIA:
        gia_str = f'_{GIA_Ylms_rate.title}'
        attributes['GIA'] = (str(GIA_Ylms_rate.citation), GIA_FILE.name)
    else:
        gia_str = ''
    # monthly GIA calculated by gia_rate*time elapsed
    # finding change in GIA each month
    GIA_Ylms = GIA_Ylms_rate.drift(GSM_Ylms.time, epoch=2003.3)
    GIA_Ylms.month[:] = np.copy(GSM_Ylms.month)
    # save geocenter coefficients of monthly GIA variability
    gia = gravtk.geocenter().from_harmonics(GIA_Ylms)

    # GRACE GAD degree 1
    GAD = gravtk.geocenter()
    GAD.time = np.copy(GAD_Ylms.time)
    GAD.month = np.copy(GAD_Ylms.month)
    GAD.C10 = np.zeros((n_files))
    GAD.C11 = np.zeros((n_files))
    GAD.S11 = np.zeros((n_files))
    for t in range(0,n_files):
        # converting GAD degree 1 harmonics to mass
        # NOTE: following Swenson (2008): do not use the kl Load Love number
        # to convert the GAD coefficients into coefficients of mass as
        # the GAC and GAD products are computed with a Load Love number of 0
        GAD.C10[t] = rho_e*rad_e*np.squeeze(GAD_Ylms.clm[1,0,t])*(2.0 + 1.0)/3.0
        GAD.C11[t] = rho_e*rad_e*np.squeeze(GAD_Ylms.clm[1,1,t])*(2.0 + 1.0)/3.0
        GAD.S11[t] = rho_e*rad_e*np.squeeze(GAD_Ylms.slm[1,1,t])*(2.0 + 1.0)/3.0
    # removing the mean of the GAD OBP coefficients
    GAD.mean(apply=True)

    # read atmospheric jump corrections from Fagiolini et al. (2015)
    ATM_Ylms = GSM_Ylms.zeros_like()
    ATM_Ylms.time[:] = np.copy(GSM_Ylms.time)
    ATM_Ylms.month[:] = np.copy(GSM_Ylms.month)
    if ATM:
        atm_corr = gravtk.read_ecmwf_corrections(base_dir,
            LMAX, ATM_Ylms.month)
        ATM_Ylms.clm[:,:,:] = np.copy(atm_corr['clm'])
        ATM_Ylms.slm[:,:,:] = np.copy(atm_corr['slm'])
        # removing the mean of the atmospheric jump correction coefficients
        ATM_Ylms.mean(apply=True)
    # truncate to degree and order LMAX/MMAX
    ATM_Ylms = ATM_Ylms.truncate(lmax=LMAX, mmax=MMAX)
    # save geocenter coefficients of the atmospheric jump corrections
    atm = gravtk.geocenter().from_harmonics(ATM_Ylms)

    # read bottom pressure model if applicable
    if MODEL not in ('OMCT','MPIOM'):
        # read input data files for ascii (txt), netCDF4 (nc) or HDF5 (H5)
        MODEL_INDEX = pathlib.Path(MODEL_INDEX).expanduser().absolute()
        OBP_Ylms = gravtk.harmonics().from_index(MODEL_INDEX,
            format=DATAFORM)
        attributes['lineage'].extend([f.name for f in OBP_Ylms.filename])
        # reduce to GRACE/GRACE-FO months and truncate to degree and order
        OBP_Ylms = OBP_Ylms.subset(GSM_Ylms.month).truncate(lmax=LMAX,mmax=MMAX)
        # filter ocean bottom pressure coefficients
        if DESTRIPE:
            OBP_Ylms = OBP_Ylms.destripe()
        # removing the mean of the ecco spherical harmonic coefficients
        OBP_Ylms.mean(apply=True)
        # converting ecco degree 1 harmonics to coefficients of mass
        OBP = gravtk.geocenter.from_harmonics(OBP_Ylms).scale(dfactor[1])

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    m = GSM_Ylms.m[:, np.newaxis]
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

    # get seasonal variations of an initial geocenter correction
    # for use in the land water mass calculation
    seasonal_geocenter = model_seasonal_geocenter(tdec)

    # iterate solutions: if not single iteration
    n_iter = 0
    eps = np.inf
    eps_max = 1e-6

    # Calculating data matrices
    # GRACE Eustatic degree 1 from land variations
    eustatic = gravtk.geocenter()
    eustatic.C10 = np.zeros((n_files))
    eustatic.C11 = np.zeros((n_files))
    eustatic.S11 = np.zeros((n_files))
    # Allocate for G matrix parameters
    # G matrix calculates the GRACE ocean mass variations
    G = gravtk.geocenter()
    G.C10 = np.zeros((n_files))
    G.C11 = np.zeros((n_files))
    G.S11 = np.zeros((n_files))
    # DMAT is the degree one matrix ((C10,C11,S11) x Time) in terms of mass
    DMAT = np.zeros((3,n_files))
    # degree 1 iterations
    iteration = gravtk.geocenter()
    iteration.C10 = np.zeros((n_files,max_iter))
    iteration.C11 = np.zeros((n_files,max_iter))
    iteration.S11 = np.zeros((n_files,max_iter))
    # calculate non-iterated terms for each file (G-matrix parameters)
    for t in range(n_files):
        # calculate geocenter component of ocean mass with GRACE
        # allocate for product of grace and legendre polynomials
        pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
        psin = np.zeros((MMAX+1, nlat))#-[m,lat]
        # Summing product of plms and c/slms over all SH degrees >= 2
        # Removing monthly GIA signal and atmospheric correction
        Ylms = GSM_Ylms.index(t)
        Ylms.subtract(GIA_Ylms.index(t))
        Ylms.subtract(ATM_Ylms.index(t))
        for i in range(0, nlat):
            l = np.arange(2,LMAX+1)
            pcos[:,i] = np.sum(plmout[l,:,i]*Ylms.clm[l,:], axis=0)
            psin[:,i] = np.sum(plmout[l,:,i]*Ylms.slm[l,:], axis=0)
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
            # calculate eustatic component from GRACE (can iterate)
            if (n_iter == 0):
                # for first iteration (will be only iteration if not ITERATIVE):
                # seasonal component of geocenter variation for land water
                GSM_Ylms.clm[1,0,t] = seasonal_geocenter.C10[t]
                GSM_Ylms.clm[1,1,t] = seasonal_geocenter.C11[t]
                GSM_Ylms.slm[1,1,t] = seasonal_geocenter.S11[t]
            else:
                # for all others: use previous iteration of inversion
                # for each of the geocenter solutions (C10, C11, S11)
                GSM_Ylms.clm[1,0,t] = iteration.C10[t,n_iter-1]
                GSM_Ylms.clm[1,1,t] = iteration.C11[t,n_iter-1]
                GSM_Ylms.slm[1,1,t] = iteration.S11[t,n_iter-1]

            # allocate for product of grace and legendre polynomials
            pcos = np.zeros((MMAX+1, nlat))#-[m,lat]
            psin = np.zeros((MMAX+1, nlat))#-[m,lat]
            # Summing product of plms and c/slms over all SH degrees
            # Removing monthly GIA signal and atmospheric correction
            Ylms = GSM_Ylms.index(t)
            Ylms.subtract(GIA_Ylms.index(t))
            Ylms.subtract(ATM_Ylms.index(t))
            for i in range(0, nlat):
                # for land water: use an initial seasonal geocenter estimate
                # from Chen et al. (1999) then the iterative if specified
                l = np.arange(1,LMAX+1)
                pcos[:,i] = np.sum(plmout[l,:,i]*Ylms.clm[l,:], axis=0)
                psin[:,i] = np.sum(plmout[l,:,i]*Ylms.slm[l,:], axis=0)

            # Multiplying by c/s(phi#m) to get surface density in cm w.e. (lonxlat)
            # this will be a spatial field similar to outputs from stokes_combine.py
            # ccos/ssin are mXphi, pcos/psin are mXtheta: resultant matrices are phiXtheta
            # The summation over spherical harmonic order is in this multiplication
            lmass = np.dot(np.transpose(ccos),pcos) + np.dot(np.transpose(ssin),psin)

            # use sea level fingerprints or eustatic from GRACE land components
            if FINGERPRINT:
                # calculate total sea level fingerprint for eustatic component
                # steps to calculate sea level from GRACE land-water change:
                # 1) calculate total land mass at time t (GRACE*land function)
                # NOTE: this is an unscaled GRACE estimate that uses the
                # buffered land function when solving the sea-level equation.
                # possible improvement using scaled estimate with real coastlines
                land_Ylms = gravtk.gen_stokes(lmass*land_function,
                    landsea.lon, landsea.lat, UNITS=1, LMIN=0,
                    LMAX=EXPANSION, PLM=PLM, LOVE=LOVE)
                # 2) calculate sea level fingerprints of land mass at time t
                # use maximum of 3 iterations for computational efficiency
                sea_level = gravtk.sea_level_equation(land_Ylms.clm, land_Ylms.slm,
                    landsea.lon, landsea.lat, land_function, LMAX=EXPANSION,
                    LOVE=LOVE, BODY_TIDE_LOVE=0, FLUID_LOVE=0, ITERATIONS=3,
                    POLAR=True, PLM=PLM, ASTYPE=np.float64, SCALE=1e-32,
                    FILL_VALUE=0)
                # 3) convert sea level fingerprints into spherical harmonics
                slf_Ylms = gravtk.gen_stokes(sea_level, landsea.lon, landsea.lat,
                    UNITS=1, LMIN=0, LMAX=1, PLM=PLM[:2,:2,:], LOVE=LOVE)
                # 4) convert the slf degree 1 harmonics to mass with dfactor
                eustatic.C10[t] = slf_Ylms.clm[1,0]*dfactor[1]
                eustatic.C11[t] = slf_Ylms.clm[1,1]*dfactor[1]
                eustatic.S11[t] = slf_Ylms.slm[1,1]*dfactor[1]
            else:
                # steps to calculate eustatic component from GRACE land-water change:
                # 1) calculate total mass of 1 cm of ocean height (calculated above)
                # 2) calculate total land mass at time t (GRACE*land function)
                # NOTE: possible improvement using the sea-level equation to solve
                # for the spatial pattern of sea level from the land water mass
                land_Ylms = gravtk.gen_stokes(lmass*land_function,
                    landsea.lon, landsea.lat, UNITS=1, LMIN=0, LMAX=1,
                    PLM=PLM[:2,:2,:], LOVE=LOVE)
                # 3) calculate ratio between the total land mass and the total mass
                # of 1 cm of ocean height (negative as positive land = sea level drop)
                # this converts the total land change to ocean height change
                eustatic_ratio = -land_Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
                # 4) scale degree one coefficients of ocean function with ratio
                # and convert the eustatic degree 1 harmonics to mass with dfactor
                scale_factor = eustatic_ratio*dfactor[1]
                eustatic.C10[t] = ocean_Ylms.clm[1,0]*scale_factor
                eustatic.C11[t] = ocean_Ylms.clm[1,1]*scale_factor
                eustatic.S11[t] = ocean_Ylms.slm[1,1]*scale_factor

            # eustatic coefficients of degree 1
            # for OMCT/MPIOM:
            # equal to the eustatic component only as OMCT/MPIOM model is
            # already removed from the GRACE/GRACE-FO GSM coefficients
            CMAT = np.array([eustatic.C10[t],eustatic.C11[t],eustatic.S11[t]])
            # replacing the OBP harmonics of degree 1
            if MODEL not in ('OMCT','MPIOM'):
                # calculate difference between ECCO and GAD as the OMCT/MPIOM
                # model is already removed from the GRACE GSM coefficients
                GADMAT = np.array([GAD.C10[t], GAD.C11[t], GAD.S11[t]])
                OBPMAT = np.array([OBP.C10[t], OBP.C11[t], OBP.S11[t]])
                # effectively adding back OMCT/MPIOM and then removing ECCO
                CMAT += OBPMAT - GADMAT

            # G Matrix for time t
            GMAT = np.array([G.C10[t], G.C11[t], G.S11[t]])
            # calculate degree 1 solution for iteration
            # this is mathematically equivalent to an iterative procedure
            # whereby the initial degree one coefficients are used to update
            # the G Matrix until (C10, C11, S11) converge
            # for OMCT/MPIOM: min(eustatic from land - measured ocean)
            # for ECCO: min((OBP-GAD) + eustatic from land - measured ocean)
            if (SOLVER == 'inv'):
                DMAT[:,t] = np.dot(np.linalg.inv(IMAT), (CMAT-GMAT))
            elif (SOLVER == 'lstsq'):
                DMAT[:,t] = np.linalg.lstsq(IMAT, (CMAT-GMAT), rcond=-1)[0]
            elif SOLVER in ('gelsd', 'gelsy', 'gelss'):
                DMAT[:,t], res, rnk, s = scipy.linalg.lstsq(IMAT, (CMAT-GMAT),
                    lapack_driver=SOLVER)
            # save geocenter for iteration and time t after restoring GIA+ATM
            iteration.C10[t,n_iter] = DMAT[0,t]/dfactor[1]+gia.C10[t]+atm.C10[t]
            iteration.C11[t,n_iter] = DMAT[1,t]/dfactor[1]+gia.C11[t]+atm.C11[t]
            iteration.S11[t,n_iter] = DMAT[2,t]/dfactor[1]+gia.S11[t]+atm.S11[t]
        # remove mean of each solution for iteration
        iteration.C10[:,n_iter] -= iteration.C10[:,n_iter].mean()
        iteration.C11[:,n_iter] -= iteration.C11[:,n_iter].mean()
        iteration.S11[:,n_iter] -= iteration.S11[:,n_iter].mean()
        # calculate difference between original geocenter coefficients and the
        # calculated coefficients for each of the geocenter solutions
        sigma_C10 = np.sum((GSM_Ylms.clm[1,0,:] - iteration.C10[:,n_iter])**2)
        sigma_C11 = np.sum((GSM_Ylms.clm[1,1,:] - iteration.C11[:,n_iter])**2)
        sigma_S11 = np.sum((GSM_Ylms.slm[1,1,:] - iteration.S11[:,n_iter])**2)
        power = GSM_Ylms.clm[1,0,:]**2 + GSM_Ylms.clm[1,1,:]**2 + GSM_Ylms.slm[1,1,:]**2
        eps = np.sqrt(sigma_C10 + sigma_C11 + sigma_S11)/np.sqrt(np.sum(power))
        # add 1 to n_iter counter
        n_iter += 1

    # Convert inverted solutions into fully normalized spherical harmonics
    # restore geocenter variation from glacial isostatic adjustment (GIA)
    # restore atmospheric jump corrections from Fagiolini (2015) if applicable
    # for each of the geocenter solutions (C10, C11, S11)
    # for the iterative case this will be the final iteration
    DEG1 = gravtk.geocenter()
    DEG1.C10 = DMAT[0,:]/dfactor[1] + gia.C10[:] + atm.C10[:]
    DEG1.C11 = DMAT[1,:]/dfactor[1] + gia.C11[:] + atm.C11[:]
    DEG1.S11 = DMAT[2,:]/dfactor[1] + gia.S11[:] + atm.S11[:]
    # remove mean of geocenter for each component
    DEG1.mean(apply=True)
    # calculate geocenter variations with dealiasing restored
    aod = DEG1.copy()
    aod.add(GAC)

    # output degree 1 coefficients with and without dealiasing
    file_format = '{0}_{1}_{2}{3}{4}{5}{6}{7}{8}.{9}'
    output_format = '{0:11.4f}{1:14.6e}{2:14.6e}{3:14.6e} {4:03d}\n'
    # public file format in fully normalized spherical harmonics
    # before and after restoring the atmospheric and oceanic dealiasing
    for AOD in ['','_wAOD']:
        # local version with all descriptor flags
        a1=(PROC,DREL,MODEL,slf_str,iter_str,slr_str,gia_str,AOD,ds_str,'txt')
        FILE1 = DIRECTORY.joinpath(file_format.format(*a1))
        fid1 = FILE1.open(mode='w', encoding='utf8')
        # print headers for cases with and without dealiasing
        print_header(fid1)
        print_harmonic(fid1,LOVE.kl[1])
        print_global(fid1,PROC,DREL,MODEL.replace('_',' '),AOD,GIA_Ylms_rate,
            SLR_C20,SLR_21,GSM_Ylms.month)
        print_variables(fid1,'single precision','fully normalized')
        # for each GRACE/GRACE-FO month
        for t,mon in enumerate(GSM_Ylms.month):
            # geocenter coefficients with and without AOD restored
            if AOD:
                args=(tdec[t],aod.C10[t],aod.C11[t],aod.S11[t],mon)
            else:
                args=(tdec[t],DEG1.C10[t],DEG1.C11[t],DEG1.S11[t],mon)
            # output geocenter coefficients to file
            fid1.write(output_format.format(*args))
        # close the output file
        fid1.close()
        output_files.append(FILE1)
        # set the permissions mode of the output file
        FILE1.chmod(mode=MODE)
        # create public and archival copies of data
        if COPY:
            # create symbolic link for public distribution without flags
            a2=(PROC,DREL,MODEL,slf_str,iter_str,'','',AOD,'','txt')
            FILE2 = DIRECTORY.joinpath(file_format.format(*a2))
            os.symlink(FILE1,FILE2) if not FILE2.exists() else None
            output_files.append(FILE2)
            # create copy of file with date for archiving
            today = time.strftime('_%Y-%m-%d',time.localtime())
            a3=(PROC,DREL,MODEL,slf_str,iter_str,'','',AOD,today,'txt')
            FILE3 = DIRECTORY.joinpath(file_format.format(*a3))
            shutil.copyfile(FILE1,FILE3)
            # copy modification times and permissions for archive file
            shutil.copystat(FILE1,FILE3)
            output_files.append(FILE3)

    # output all degree 1 coefficients as a netCDF4 file
    a4=(PROC,DREL,MODEL,slf_str,iter_str,slr_str,gia_str,'',ds_str,'nc')
    FILE4 = DIRECTORY.joinpath(file_format.format(*a4))
    fileID = netCDF4.Dataset(FILE4, mode='w')
    # Defining the NetCDF4 dimensions
    fileID.createDimension('iteration', n_iter)
    fileID.createDimension('time', n_files)

    # variable attributes
    attrs = dict(time={}, month={}, C10={}, C11={}, S11={})
    attrs['time']['units'] = 'years'
    attrs['time']['long_name'] = 'Date_in_Decimal_Years'
    attrs['month']['long_name'] = 'GRACE_month'
    attrs['month']['units'] = 'months since 2001-12-01'
    attrs['month']['calendar'] = 'standard'
    attrs['C10']['units'] = 'fully_normalized'
    attrs['C10']['long_name'] = 'cosine_spherical_harmonic_of_degree_1,_order_0'
    attrs['C10']['description'] = 'spherical_harmonics'
    attrs['C11']['units'] = 'fully_normalized'
    attrs['C11']['long_name'] = 'cosine_spherical_harmonic_of_degree_1,_order_1'
    attrs['C11']['description'] = 'spherical_harmonics'
    attrs['S11']['units'] = 'fully_normalized'
    attrs['S11']['long_name'] = 'sine_spherical_harmonic_of_degree_1,_order_1'
    attrs['S11']['description'] = 'spherical_harmonics'

    # defining the NetCDF4 variables
    nc = {}
    nc['time'] = fileID.createVariable('time', tdec.dtype, ('time',))
    nc['month'] = fileID.createVariable('month', months.dtype, ('time',))
    # filling NetCDF4 variables
    nc['time'][:] = tdec[:].copy()
    nc['month'][:] = months[:].copy()
    # set attributes for time and month
    for key in ('time','month'):
        for att_name, att_val in attrs[key].items():
            nc[key].setncattr(att_name, att_val)

    # degree 1 coefficients from the iterative solution
    for key in iteration.fields:
        var = iteration.get(key)
        nc[key] = fileID.createVariable(key, var.dtype,
            ('time','iteration',), zlib=True)
        nc[key][:] = var[:,:n_iter]
        for att_name, att_val in attrs[key].items():
            nc[key].setncattr(att_name, att_val)

    # atmospheric and oceanic dealiasing
    nc['AOD'] = {}
    g1 = fileID.createGroup('AOD')
    g1.description = f'Atmospheric and oceanic dealiasing'
    gac = GAC.scale(1.0/dfactor[1])
    for key in gac.fields:
        var = gac.get(key)
        nc['AOD'][key] = g1.createVariable(key, var.dtype,
            ('time',), zlib=True)
        nc['AOD'][key][:] = var[:]
        for att_name, att_val in attrs[key].items():
            nc['AOD'][key].setncattr(att_name, att_val)

    # ocean bottom pressure
    nc['OBP'] = {}
    g2 = fileID.createGroup('OBP')
    g2.description = f'Ocean bottom pressure from {MODEL}'
    if MODEL not in ('OMCT','MPIOM'):
        obp = OBP.scale(1.0/dfactor[1])
    else:
        obp = GAD.scale(1.0/dfactor[1])
    for key in obp.fields:
        var = obp.get(key)
        nc['OBP'][key] = g2.createVariable(key, var.dtype,
            ('time',), zlib=True)
        nc['OBP'][key][:] = var[:]
        for att_name, att_val in attrs[key].items():
            nc['OBP'][key].setncattr(att_name, att_val)

    # GRACE components of ocean water mass
    nc['OWM'] = {}
    g3 = fileID.createGroup('OWM')
    g3.description = f'Ocean water mass from {MISSION[DREL]}'
    owm = G.scale(1.0/dfactor[1])
    for key in owm.fields:
        var = owm.get(key)
        nc['OWM'][key] = g3.createVariable(key, var.dtype,
            ('time',), zlib=True)
        nc['OWM'][key][:] = var[:]
        for att_name, att_val in attrs[key].items():
            nc['OWM'][key].setncattr(att_name, att_val)

    # eustatic sea level from land water mass
    nc['ESL'] = {}
    g4 = fileID.createGroup('ESL')
    g4.description = 'Eustatic sea level from land water mass'
    esl = eustatic.scale(1.0/dfactor[1])
    for key in esl.fields:
        var = esl.get(key)
        nc['ESL'][key] = g4.createVariable(key, var.dtype,
            ('time',), zlib=True)
        nc['ESL'][key][:] = var[:]
        for att_name, att_val in attrs[key].items():
            nc['ESL'][key].setncattr(att_name, att_val)

    # define global attributes
    for att_name, att_val in attributes.items():
        fileID.setncattr(att_name, att_val)
    # define creation date attribute
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
    # close the output file
    fileID.close()
    # set the permissions mode of the output file
    FILE4.chmod(mode=MODE)
    output_files.append(FILE4)

    # create plot similar to Figure 1 of Swenson et al (2008)
    if PLOT:
        # 3 row plot (C10, C11 and S11) with:
        # - ECCO OBP geocenter
        # - OMCT/MPIOM OBP geocenter
        # - eustatic sea level geocenter
        # - G-matrix ocean water mass components
        ax = {}
        fig, (ax[0], ax[1], ax[2]) = plt.subplots(num=1, nrows=3,
            sharex=True, sharey=True, figsize=(6,9))
        ii = np.nonzero((tdec >= 2003.) & (tdec < 2008.))
        # remove means of individual geocenter components
        if MODEL not in ('OMCT','MPIOM'):
            OBP.mean(apply=True,indices=ii)
        G.mean(apply=True,indices=ii)
        GAD.mean(apply=True,indices=ii)
        eustatic.mean(apply=True,indices=ii)
        for i,key in enumerate(G.fields):
            # plot ocean bottom pressure for alternative models
            if MODEL not in ('OMCT','MPIOM'):
                ax[i].plot(tdec, 10.*OBP.get(key), color='#1ed565', lw=2)
            # plot GRACE components
            ax[i].plot(tdec, 10.*G.get(key), color='orange', lw=2)
            # plot OMCT/MPIOM ocean bottom pressure
            ax[i].plot(tdec, 10.*GAD.get(key), color='blue', lw=2)
            # plot eustatic components
            ax[i].plot(tdec, 10.*eustatic.get(key), color='r', lw=2)
            ax[i].set_ylabel('[mm]', fontsize=14)
            # add axis labels and adjust font sizes for axis ticks
            # axis label
            artist = matplotlib.offsetbox.AnchoredText(key, pad=0.,
                prop=dict(size=16,weight='bold'), frameon=False, loc=2)
            ax[i].add_artist(artist)
            # axes tick adjustments
            for tick in ax[i].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax[i].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            # adjust ticks
            ax[i].get_xaxis().set_tick_params(which='both', direction='in')
            ax[i].get_yaxis().set_tick_params(which='both', direction='in')
        # labels and set limits to Swenson range
        ax[2].set_xlabel('Time [Yr]', fontsize=14)
        ax[2].set_xlim(2003, 2007)
        ax[2].set_ylim(-6, 6)
        ax[2].xaxis.set_ticks(np.arange(2003, 2008, 1))
        ax[2].xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
        ax[2].yaxis.set_ticks(np.arange(-6, 8, 2))
        ax[2].xaxis.get_major_formatter().set_useOffset(False)
        # adjust locations of subplots and save to file
        fig.subplots_adjust(left=0.1,right=0.96,bottom=0.06,top=0.98,hspace=0.1)
        args = (PROC,DREL,MODEL,slf_str,iter_str,slr_str,gia_str,ds_str)
        FILE = 'Swenson_Figure_1_{0}_{1}_{2}{3}{4}{5}{6}{7}.pdf'.format(*args)
        PLOT1 = DIRECTORY.joinpath(FILE)
        metadata = {'Title':pathlib.Path(sys.argv[0]).name}
        plt.savefig(PLOT1, format='pdf', metadata=metadata)
        plt.clf()
        # set the permissions mode of the output files
        PLOT1.chmod(mode=MODE)
        output_files.append(PLOT1)

    # if ITERATIVE: create plot showing iteration solutions
    if PLOT and ITERATIVE:
        # 3 row plot (C10, C11 and S11)
        ax = {}
        fig, (ax[0], ax[1], ax[2]) = plt.subplots(num=2, nrows=3,
            sharex=True, figsize=(6,9))
        # show solutions for each iteration
        cmap = copy.copy(cm.rainbow)
        plot_colors = iter(cmap(np.linspace(0,1,n_iter)))
        iteration_mmwe = iteration.scale(10.0*dfactor[1])
        for j in range(n_iter):
            c = next(plot_colors)
            # C10, C11 and S11
            ax[0].plot(GSM_Ylms.month,iteration_mmwe.C10[:,j],c=c)
            ax[1].plot(GSM_Ylms.month,iteration_mmwe.C11[:,j],c=c)
            ax[2].plot(GSM_Ylms.month,iteration_mmwe.S11[:,j],c=c)
        # add axis labels and adjust font sizes for axis ticks
        for i,key in enumerate(iteration_mmwe.fields):
            ax[i].set_ylabel('mm', fontsize=14)
            # axis label
            artist = matplotlib.offsetbox.AnchoredText(key, pad=0.,
                prop=dict(size=16,weight='bold'), frameon=False, loc=2)
            ax[i].add_artist(artist)
            # axes tick adjustments
            for tick in ax[i].xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax[i].yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            # adjust ticks
            ax[i].get_xaxis().set_tick_params(which='both', direction='in')
            ax[i].get_yaxis().set_tick_params(which='both', direction='in')
        # labels and set limits
        ax[2].set_xlabel('Grace Month', fontsize=14)
        xmin = np.floor(GSM_Ylms.month[0]/10.)*10.
        xmax = np.ceil(GSM_Ylms.month[-1]/10.)*10.
        ax[2].set_xlim(xmin,xmax)
        ax[2].xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax[2].xaxis.get_major_formatter().set_useOffset(False)
        # adjust locations of subplots and save to file
        fig.subplots_adjust(left=0.12,right=0.94,bottom=0.06,top=0.98,hspace=0.1)
        args = (PROC,DREL,MODEL,slf_str,slr_str,gia_str,ds_str)
        FILE = 'Geocenter_Iterative_{0}_{1}_{2}{3}{4}{5}{6}.pdf'.format(*args)
        PLOT2 = DIRECTORY.joinpath(FILE)
        metadata = {'Title':pathlib.Path(sys.argv[0]).name}
        plt.savefig(PLOT2, format='pdf', metadata=metadata)
        plt.clf()
        # set the permissions mode of the output files
        PLOT2.chmod(mode=MODE)
        output_files.append(PLOT2)

    # return the list of output files and the number of iterations
    return (output_files, n_iter)

# PURPOSE: print YAML header to top of file
def print_header(fid):
    # print header
    fid.write('{0}:\n'.format('header'))
    # data dimensions
    fid.write('  {0}:\n'.format('dimensions'))
    fid.write('    {0:22}: {1:d}\n'.format('degree',1))
    fid.write('    {0:22}: {1:d}\n'.format('order',1))
    fid.write('\n')

# PURPOSE: print spherical harmonic attributes to YAML header
def print_harmonic(fid,kl):
    # non-standard attributes
    fid.write('  {0}:\n'.format('non-standard_attributes'))
    # load love number
    fid.write('    {0:22}:\n'.format('love_number'))
    long_name = 'Gravitational Load Love Number of Degree 1 (k1)'
    fid.write('      {0:20}: {1}\n'.format('long_name',long_name))
    fid.write('      {0:20}: {1:0.3f}\n'.format('value',kl))
    # data format
    data_format = '(f11.4,3e14.6,i4)'
    fid.write('    {0:22}: {1}\n'.format('formatting_string',data_format))
    fid.write('\n')

# PURPOSE: print global attributes to YAML header
def print_global(fid,PROC,DREL,MODEL,AOD,GIA,SLR,S21,month):
    fid.write('  {0}:\n'.format('global_attributes'))
    MISSION = dict(RL05='GRACE',RL06='GRACE/GRACE-FO')
    title = '{0} Geocenter Coefficients {1} {2}'.format(MISSION[DREL],PROC,DREL)
    fid.write('    {0:22}: {1}\n'.format('title',title))
    summary = []
    summary.append(('Geocenter coefficients derived from {0} mission '
        'measurements and {1} ocean model outputs.').format(MISSION[DREL],MODEL))
    if AOD:
        summary.append(('  These coefficients represent the largest-scale '
            'variability of atmospheric, oceanic, hydrologic, cryospheric, '
            'and solid Earth processes.'))
    else:
        summary.append(('  These coefficients represent the largest-scale '
            'variability of hydrologic, cryospheric, and solid Earth '
            'processes.  In addition, the coefficients represent the '
            'atmospheric and oceanic processes not captured in the {0} {1} '
            'de-aliasing product.').format(MISSION[DREL],DREL))
    # get GIA parameters
    summary.append(('  Glacial Isostatic Adjustment (GIA) estimates from '
        '{0} have been restored.').format(GIA.citation))
    if AOD:
        summary.append(('  Monthly atmospheric and oceanic de-aliasing product '
            'has been restored.'))
    elif (DREL == 'RL05') and not AOD:
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
    if AOD:
        src += ('  Atmospheric and oceanic variation restored using the {2} '
            'de-aliasing product.')
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
    PRODUCT_VERSION = f'Release-{DREL[2:]}'
    fid.write('    {0:22}: {1}\n'.format('product_version',PRODUCT_VERSION))
    fid.write('    {0:22}:\n'.format('references'))
    reference = []
    # geocenter citations
    reference.append(('T. C. Sutterley, and I. Velicogna, "Improved estimates '
        'of geocenter variability from time-variable gravity and ocean model '
        'outputs", Remote Sensing, 11(18), 2108, (2019). '
        'https://doi.org/10.3390/rs11182108'))
    reference.append(('S. C. Swenson, D. P. Chambers, and J. Wahr, "Estimating '
        'geocenter variations from a combination of GRACE and ocean model '
        'output", Journal of Geophysical Research - Solid Earth, 113(B08410), '
        '(2008). https://doi.org/10.1029/2007JB005338'))
    # GIA citation
    reference.append(GIA.reference)
    # ECMWF jump corrections citation
    if (DREL == 'RL05') and not AOD:
        reference.append(('E. Fagiolini, F. Flechtner, M. Horwath, H. Dobslaw, '
            '''"Correction of inconsistencies in ECMWF's operational '''
            '''analysis data during de-aliasing of GRACE gravity models", '''
            'Geophysical Journal International, 202(3), 2150, (2015). '
            'https://doi.org/10.1093/gji/ggv276'))
    # SLR citation for a given solution
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
    # print list of references
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
    # date range and date created
    calendar_year,calendar_month = gravtk.time.grace_to_calendar(month)
    start_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[0],calendar_month[0])
    fid.write('    {0:22}: {1}\n'.format('time_coverage_start', start_time))
    end_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[-1],calendar_month[-1])
    fid.write('    {0:22}: {1}\n'.format('time_coverage_end', end_time))
    today = time.strftime('%Y-%m-%d',time.localtime())
    fid.write('    {0:22}: {1}\n'.format('date_created', today))
    fid.write('\n')

# PURPOSE: print variable descriptions to YAML header
def print_variables(fid,data_precision,data_units):
    # variables
    fid.write('  {0}:\n'.format('variables'))
    # time
    fid.write('    {0:22}:\n'.format('mid-epoch_time'))
    long_name = 'mid-date of each measurement epoch'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
    fid.write('      {0:20}: {1}\n'.format('units', 'decimal-years'))
    fid.write('      {0:20}: {1}\n'.format('comment', '1st column'))
    # C10
    fid.write('    {0:22}:\n'.format('C10'))
    long_name = 'C10 coefficient; cosine coefficient for degree 1 and order 0'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '2nd column'))
    # C11
    fid.write('    {0:22}:\n'.format('C11'))
    long_name = 'C11 coefficient; cosine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '3rd column'))
    # S11
    fid.write('    {0:22}:\n'.format('S11'))
    long_name = 'S11 coefficient; sine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '4th column'))
    # GRACE month
    fid.write('    {0:22}:\n'.format('month'))
    long_name = 'GRACE month of each measurement epoch'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    description = 'months starting 2002-01-01'
    fid.write('      {0:20}: {1}\n'.format('description', description))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
    fid.write('      {0:20}: {1}\n'.format('units', 'month'))
    fid.write('      {0:20}: {1}\n'.format('comment', '5th column'))
    # end of header
    fid.write('\n\n# End of YAML header\n')

# PURPOSE: print a file log for the GRACE degree one analysis
def output_log_file(input_arguments, output_files, n_iter):
    # format: calc_degree_one_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_degree_one_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = pathlib.Path(input_arguments.directory).joinpath('geocenter')
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

# PURPOSE: print a error file log for the GRACE degree one analysis
def output_error_log_file(input_arguments):
    # format: calc_degree_one_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_degree_one_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Calculates degree 1 variations using GRACE/GRACE-FO
            coefficients of degree 2 and greater, and ocean bottom pressure
            variations from ECCO and OMCT/MPIOM
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # start and end GRACE/GRACE-FO months
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
        type=float, default=0.021, nargs='?',
        help='Degree 1 gravitational Load Love number')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # GIA model type list
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
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', default='AW13-ICE6G', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
    # use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    # correct for pole tide drift follow Wahr et al. (2015)
    parser.add_argument('--pole-tide',
        default=False, action='store_true',
        help='Correct for pole tide drift')
    # replace low degree harmonics with values from Satellite Laser Ranging
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
    parser.add_argument('--slr-c40',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C40 coefficients with SLR values')
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    # ocean model list
    choices = []
    choices.append('OMCT')
    choices.append('MPIOM')
    choices.append('ECCO_kf080i')
    choices.append('ECCO_dr080i')
    choices.append('ECCO_V4r3')
    choices.append('ECCO_V4r4')
    choices.append('ECCO_V5alpha')
    parser.add_argument('--ocean-model',
        metavar='MODEL', type=str,
        default='MPIOM', choices=choices,
        help='Ocean model to use')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for ocean models')
    # index file for ocean model harmonics
    parser.add_argument('--ocean-file',
        type=pathlib.Path,
        help='Index file for ocean model harmonics')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=pathlib.Path,
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
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
    # copy output files
    parser.add_argument('--copy','-C',
        default=False, action='store_true',
        help='Copy output files for distribution and archival')
    # Output log file for each job in forms
    # calc_degree_one_run_2002-04-01_PID-00000.log
    # calc_degree_one_failed_run_2002-04-01_PID-00000.log
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
        # run calc_degree_one algorithm with parameters
        output_files,n_iter = calc_degree_one(
            args.directory,
            args.center,
            args.release,
            args.ocean_model,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
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
            SLR_C40=args.slr_c40,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            MODEL_INDEX=args.ocean_file,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            ITERATIVE=args.iterative,
            SOLVER=args.solver,
            FINGERPRINT=args.fingerprint,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            PLOT=args.plot,
            COPY=args.copy,
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
