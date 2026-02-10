#!/usr/bin/env python
u"""
kernel_degree_one_error.py
Written by Tyler Sutterley (06/2024)

Calculates degree 1 errors using GRACE coefficients of degree 2 and greater,
    and ocean bottom pressure variations from OMCT/MPIOM

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
    -F X, --format X: Input/output data format for harmonics files
        ascii
        netCDF4
        HDF5
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
    --fingerprint: Redistribute land-water flux using sea level fingerprints
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
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces low-degree harmonics with SLR values (if specified)
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    time_series.smooth.py: smoothes a time-series using a Loess-type algorithm
    harmonics.py: class for processing GRACE/GRACE-FO spherical harmonic data
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    time.py: utilities for calculating time operations
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
    Updated 06/2024: use wrapper to importlib for optional dependencies
    Updated 10/2023: generalize mission variable to be GRACE/GRACE-FO
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 02/2023: use love numbers class with additional attributes
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 05/2022: use argparse descriptions within documentation
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: calculate uncertainty using meshgrid of harmonic arrays
        add YAML headers to output geocenter uncertainty file
        can use variable loglevels for verbose output
    Written 11/2021
"""
from __future__ import print_function

import sys
import os
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
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

# PURPOSE: calculate the uncertainty in the geocenter time-series
def kernel_degree_one_error(base_dir, PROC, DREL, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    LOVE_K1=None,
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
    ERROR_FILES=[],
    FINGERPRINT=False,
    EXPANSION=None,
    ERROR_TYPE='RMS',
    MODE=0o775):

    # output directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('geocenter')
    # create output directory if non-existent
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # list object of output files for file logs (full path)
    output_files = []

    # GRACE/GRACE-FO dataset
    DSET = 'GSM'
    # do not import degree 1 coefficients
    DEG1 = ''

    # delta coefficients flag for monte carlo run
    delta_str = '_kernel'
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = f'M{MMAX:d}' if MMAX and (MMAX != LMAX) else ''
    # output flag for using sea level fingerprints
    slf_str = '_SLF' if FINGERPRINT else ''
    # filter GRACE/GRACE-FO coefficients
    ds_str = '_FL' if DESTRIPE else ''
    # atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    # ocean model string
    model_str = 'MPIOM' if (DREL == 'RL06') else 'OMCT'
    # output flag for using sea level fingerprints
    slf_str = '_SLF' if FINGERPRINT else ''
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

    # suffix for input ascii, netcdf and HDF5 files
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

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

    # Gaussian Smoothing (Jekeli, 1981)
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # reading GRACE months for input date range
    # replacing low-degree harmonics with SLR values if specified
    # correcting for Pole-Tide drift if specified
    # atmospheric jumps will be corrected externally if specified
    Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX, SLR_21=SLR_21,
        SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C40=SLR_C40, SLR_C50=SLR_C50,
        POLE_TIDE=POLE_TIDE, ATM=False, MODEL_DEG1=False)
    # create harmonics object from GRACE/GRACE-FO data
    GSM_Ylms = gravtk.harmonics().from_dict(Ylms)
    # use a mean file for the static field to remove
    if MEAN_FILE:
        # read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=MEANFORM, date=False)
        # remove the input mean
        GSM_Ylms.subtract(mean_Ylms)
    else:
        GSM_Ylms.mean(apply=True)
    # filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        # destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GSM_Ylms = GSM_Ylms.destripe()
    else:
        # using standard GRACE/GRACE-FO harmonics
        ds_str = ''
    # full path to directory for specific GRACE/GRACE-FO product
    GSM_Ylms.directory = pathlib.Path(Ylms['directory']).expanduser().absolute()
    # GRACE dates
    tdec = np.copy(GSM_Ylms.time)
    months = np.copy(GSM_Ylms.month)
    # number of months considered
    n_files = len(GSM_Ylms.month)

    # read atmospheric jump corrections from Fagiolini et al. (2015)
    ATM_Ylms = GSM_Ylms.zeros_like()
    ATM_Ylms.time[:] = np.copy(GSM_Ylms.time)
    ATM_Ylms.month[:] = np.copy(GSM_Ylms.month)
    if ATM:
        atm_corr = gravtk.read_ecmwf_corrections(base_dir, LMAX, ATM_Ylms.month)
        ATM_Ylms.clm[:,:,:] = np.copy(atm_corr['clm'])
        ATM_Ylms.slm[:,:,:] = np.copy(atm_corr['slm'])
        # removing the mean of the atmospheric jump correction coefficients
        ATM_Ylms.mean(apply=True)
    # truncate to degree and order LMAX/MMAX
    ATM_Ylms = ATM_Ylms.truncate(lmax=LMAX, mmax=MMAX)

    # input spherical harmonic datafiles to be used in monte carlo
    error_Ylms = []
    # for each file to be removed
    for ERROR_FILE in ERROR_FILES:
        # file in ascii, netCDF4 or HDF5 formats
        Ylms = gravtk.harmonics().from_file(ERROR_FILE, format=DATAFORM)
        # truncate to degree and order and append to list
        error_Ylms.append(Ylms.truncate(lmax=LMAX, mmax=MMAX))

    # calculating GRACE/GRACE-FO error (Wahr et al. 2006)
    # output GRACE error file (for both LMAX==MMAX and LMAX != MMAX cases)
    fargs = (PROC,DREL,DSET,LMAX,order_str,ds_str,atm_str,GSM_Ylms.month[0],
        GSM_Ylms.month[-1], suffix[DATAFORM])
    delta_format = '{0}_{1}_{2}_DELTA_CLM_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    DELTA_FILE = GSM_Ylms.directory.joinpath(delta_format.format(*fargs))
    # check full path of the GRACE directory for delta file
    # if file was previously calculated: will read file
    # else: will calculate the GRACE/GRACE-FO error
    if not DELTA_FILE.exists():
        # add output delta file to list object
        output_files.append(DELTA_FILE)

        # Delta coefficients of GRACE time series (Error components)
        delta_Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
        delta_Ylms.clm = np.zeros((LMAX+1, MMAX+1))
        delta_Ylms.slm = np.zeros((LMAX+1, MMAX+1))
        # Smoothing Half-Width (CNES is a 10-day solution)
        # All other solutions are monthly solutions (HFWTH for annual = 6)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            HFWTH = 19
        else:
            HFWTH = 6
        # Equal to the noise of the smoothed time-series
        # for each spherical harmonic order
        for m in range(0,MMAX+1):# MMAX+1 to include MMAX
            # for each spherical harmonic degree
            for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                # Delta coefficients of GRACE time series
                for cs,csharm in enumerate(['clm','slm']):
                    # calculate GRACE Error (Noise of smoothed time-series)
                    # With Annual and Semi-Annual Terms
                    val1 = getattr(GSM_Ylms, csharm)
                    smth = gravtk.time_series.smooth(tdec, val1[l,m,:],
                        HFWTH=HFWTH)
                    # number of smoothed points
                    nsmth = len(smth['data'])
                    tsmth = np.mean(smth['time'])
                    # GRACE/GRACE-FO delta Ylms
                    # variance of data-(smoothed+annual+semi)
                    val2 = getattr(delta_Ylms, csharm)
                    val2[l,m] = np.sqrt(np.sum(smth['noise']**2)/nsmth)

        # attributes for output files
        attributes = {}
        attributes['title'] = 'GRACE/GRACE-FO Spherical Harmonic Errors'
        attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
        # save GRACE/GRACE-FO delta harmonics to file
        delta_Ylms.time = np.copy(tsmth)
        delta_Ylms.month = np.int64(nsmth)
        delta_Ylms.to_file(DELTA_FILE, format=DATAFORM, **attributes)
        # set the permissions mode of the output harmonics file
        DELTA_FILE.chmod(mode=MODE)
        # append delta harmonics file to output files list
        output_files.append(DELTA_FILE)
    else:
        # read GRACE/GRACE-FO delta harmonics from file
        delta_Ylms = gravtk.harmonics().from_file(DELTA_FILE,
            format=DATAFORM)
        # truncate GRACE/GRACE-FO delta clm and slm to d/o LMAX/MMAX
        delta_Ylms = delta_Ylms.truncate(lmax=LMAX, mmax=MMAX)
        tsmth = np.squeeze(delta_Ylms.time)
        nsmth = np.int64(delta_Ylms.month)

    # read covariance file with all descriptor flags
    args = ('COV',LMAX,order_str,gw_str,slf_str,'nc')
    FILE = DIRECTORY.joinpath('{0}_L{1:d}{2}{3}{4}.{5}'.format(*args))
    COV = {}
    # opening netCDF covariance file for reading
    with netCDF4.Dataset(FILE, mode='r') as fileID:
        # copy each netCDF variables
        for key,val in fileID.variables.items():
            COV[key] = val[:].copy()

    # calculate total uncertainty
    if (ERROR_TYPE == 'RMS'):
        # add all the uncertainties in quadrature
        Ylms = delta_Ylms.power(2)
        for eYlms in error_Ylms:
            Ylms.add(eYlms.power(2))
        # calculate the RMS and then smooth harmonics
        Ylms = Ylms.power(0.5).convolve(wt)
    else:
        # add all the uncertainties
        Ylms = delta_Ylms.copy()
        for eYlms in error_Ylms:
            Ylms.add(eYlms)
        # smooth harmonics
        Ylms = Ylms.convolve(wt)

    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    # Switching between Cosine and Sine Stokes
    YLM = np.zeros((n_harm))

    # ii is a counter variable for building the column array
    ii = 0
    for cs,csharm in enumerate(['clm','slm']):
        # get harmonics
        temp = getattr(Ylms, csharm)
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # cosine and sine harmonics
                YLM[ii] = temp[l,m]
                # add 1 to counter
                ii += 1

    # output geocenter errors
    DEG1 = gravtk.geocenter()
    # create a meshgrid for multiplying each harmonic
    rows,cols = np.meshgrid(YLM, YLM)
    # reduce to unique harmonic pairs
    ii,jj = np.tril_indices(n_harm)
    # for each geocenter coefficient
    for d,gc in enumerate(['C10','C11','S11']):
        tmp = np.sum(rows[ii,jj]*cols[ii,jj]*COV['covariance'][ii,jj,d]**2)
        setattr(DEG1,gc,np.sqrt(tmp))

    # output degree 1 coefficient errors
    file_format = '{0}_{1}_{2}{3}{4}{5}{6}.{7}'
    output_format = '{0:11.4f}{1:14.6e}{2:14.6e}{3:14.6e} {4:03d}\n'
    # local version with all descriptor flags
    a1 = (PROC,DREL,model_str,slf_str,slr_str,delta_str,ds_str,'txt')
    FILE1 = DIRECTORY.joinpath(file_format.format(*a1))
    fid1 = FILE1.open(mode='w', encoding='utf8')
    logging.info(str(FILE1))
    # print headers
    print_header(fid1)
    print_harmonic(fid1,LOVE.kl[1])
    print_global(fid1,PROC,DREL,model_str.replace('_',' '),
        SLR_C20,SLR_21,months)
    print_variables(fid1,'single precision','fully normalized')
    # output geocenter coefficients to file
    args=(tdec.mean(),DEG1.C10,DEG1.C11,DEG1.S11,n_files)
    fid1.write(output_format.format(*args))
    # close the output file
    fid1.close()
    # change the permissions mode
    FILE1.chmod(mode=MODE)
    output_files.append(FILE1)

    # return the list of output files
    return output_files

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
def print_global(fid,PROC,DREL,MODEL,SLR,S21,month):
    fid.write('  {0}:\n'.format('global_attributes'))
    MISSION = 'GRACE/GRACE-FO'
    title = '{0} Geocenter Coefficients {1} {2}'.format(MISSION,PROC,DREL)
    fid.write('    {0:22}: {1}\n'.format('title',title))
    summary = []
    summary.append(('Geocenter error coefficients derived from {0} mission '
        'measurements and {1} ocean model outputs.').format(MISSION,MODEL))
    summary.append(('  These coefficients represent the largest-scale '
        'variability of hydrologic, cryospheric, and solid Earth '
        'processes.  In addition, the coefficients represent the '
        'atmospheric and oceanic processes not captured in the {0} {1} '
        'de-aliasing product.').format(MISSION,DREL))
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
    hist = '{0} Level-3 Data created at UC Irvine'.format(MISSION)
    fid.write('    {0:22}: {1}\n'.format('history',hist))
    src = 'An inversion using {0} measurements and {1} ocean model outputs.'
    args = (MISSION,MODEL,DREL)
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
    # ECMWF jump corrections citation
    if (DREL == 'RL05'):
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
    long_name = 'mid-date of measurement epoch'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
    fid.write('      {0:20}: {1}\n'.format('units', 'decimal-years'))
    fid.write('      {0:20}: {1}\n'.format('comment', '1st column'))
    # eC10
    fid.write('    {0:22}:\n'.format('eC10'))
    long_name = 'eC10 uncertainty; cosine coefficient for degree 1 and order 0'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '2nd column'))
    # eC11
    fid.write('    {0:22}:\n'.format('eC11'))
    long_name = 'eC11 uncertainty; cosine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '3rd column'))
    # eS11
    fid.write('    {0:22}:\n'.format('eS11'))
    long_name = 'eS11 uncertainty; sine coefficient for degree 1 and order 1'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
    fid.write('      {0:20}: {1}\n'.format('units', data_units))
    fid.write('      {0:20}: {1}\n'.format('comment', '4th column'))
    # GRACE month
    fid.write('    {0:22}:\n'.format('month'))
    long_name = 'Number of months'
    fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
    description = 'Number of months used in determining satellite errors'
    fid.write('      {0:20}: {1}\n'.format('description', description))
    fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
    fid.write('      {0:20}: {1}\n'.format('units', 'month'))
    fid.write('      {0:20}: {1}\n'.format('comment', '5th column'))
    # end of header
    fid.write('\n\n# End of YAML header\n')

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
        description="""Calculates degree 1 errors using GRACE/GRACE-FO
            coefficients of degree 2 and greater, and ocean bottom pressure
            variations from ECCO and OMCT/MPIOM
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
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
    # input/output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=pathlib.Path,
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    # additional error files to include
    parser.add_argument('--error-file',
        type=pathlib.Path,
        nargs='+', default=[],
        help='Additional error files to include in total uncertainty')
    # run with sea level fingerprints
    parser.add_argument('--fingerprint',
        default=False, action='store_true',
        help='Redistribute land-water flux using sea level fingerprints')
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
        # run kernel_degree_one_error algorithm with parameters
        output_files = kernel_degree_one_error(
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
            LOVE_NUMBERS=args.love,
            LOVE_K1=args.kl,
            ATM=args.atm_correction,
            POLE_TIDE=args.pole_tide,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C40=args.slr_c40,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            ERROR_FILES=args.error_file,
            FINGERPRINT=args.fingerprint,
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
