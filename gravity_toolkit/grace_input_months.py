#!/usr/bin/env python
u"""
grace_input_months.py
Written by Tyler Sutterley (09/2021)
Contributions by Hugo Lecomte and Yara Mohajerani

Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order
    and for a specified date range
Includes degree 1 with with input values (if specified)
Replaces C20 with SLR values (if specified)
Replaces C21/S21/C22/S22/C30/C50 with SLR values for months 179+ (if specified)
Corrects for ECMWF atmospheric "jumps" using the GAE, GAF and GAG files
Corrects for Pole Tide drift following Wahr et al. (2015)

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data
    PROC: Data processing center or satellite mission
        CSR: University of Texas Center for Space Research
        GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        JPL: Jet Propulsion Laboratory
        CNES: French Centre National D'Etudes Spatiales
        GRAZ: Institute of Geodesy from GRAZ University of Technology
        COSTG: Combination Service for Time-variable Gravity Fields
        Swarm: Time-variable gravity data from Swarm satellites
    DREL: GRACE/GRACE-FO/Swarm data release
    DSET: GRACE/GRACE-FO/Swarm data product
        GAA: non-tidal atmospheric correction
        GAB: non-tidal oceanic correction
        GAC: combined non-tidal atmospheric and oceanic correction
        GAD: ocean bottom pressure product
        GSM: monthly static field product
    LMAX: Upper bound of Spherical Harmonic Degrees (e.g. 60)
    start_mon: starting month to consider in analysis
    end_mon: ending month to consider in analysis
    missing: missing months to not consider in analysis
    SLR_C20: Replaces C20 with SLR values
        N: use original values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GFZ: use values from GFZ
        GSFC: use values from GSFC (TN-14)
    DEG1: Use Degree 1 coefficients
        None: No degree 1
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
            https://grace.jpl.nasa.gov/data/get-data/geocenter/
        SLR: satellite laser ranging coefficients from CSR
            http://download.csr.utexas.edu/pub/slr/geocenter/
        SLF: Sutterley and Velicogna coefficients, Remote Sensing (2019)
            https://doi.org/10.6084/m9.figshare.7388540
        Swenson: GRACE-derived coefficients from Sean Swenson
            https://doi.org/10.1029/2007JB005338
        GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
            http://gravis.gfz-potsdam.de/corrections

OUTPUTS:
    clm: GRACE/GRACE-FO cosine spherical harmonic to degree/order LMAX and MMAX
    slm: GRACE/GRACE-FO sine spherical harmonic to degree/order LMAX and MMAX
    eclm: GRACE/GRACE-FO uncalibrated cosine spherical harmonic errors
    eslm: GRACE/GRACE-FO uncalibrated sine spherical harmonic errors
    time: time of each GRACE/GRACE-FO measurement (mid-month)
    month: GRACE/GRACE-FO months of input datasets
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    title: string denoting low degree zonals, geocenter and corrections
    directory: directory of exact GRACE/GRACE-FO product

OPTIONS:
    MMAX: Upper bound of Spherical Harmonic Orders (default=LMAX)
    SLR_21: replaces C21 and S21 with SLR values
        None: use original values
        CSR: use values from CSR
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC
    SLR_22: replaces C22 and S22 with SLR values
        None: use original values
        CSR: use values from CSR
    SLR_C30: replaces C30 with SLR values
        None: use original values
        CSR: use values from CSR (5x5 with 6,1)
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC (TN-14)
    SLR_C50: replaces C50 with SLR values
        None: use original values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    POLE_TIDE: correct GSM data with pole tides following Wahr et al (2015)
    ATM: correct data with ECMWF "jump" corrections GAE, GAF and GAG
    MODEL_DEG1: least-squares model missing degree 1 coefficients (True/False)
    DEG1_GIA: GIA-correction used when calculating degree 1 coefficients

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    grace_date.py: reads GRACE index file and calculates dates for each month
    read_SLR_C20.py: reads C20 files from satellite laser ranging (CSR or GSFC)
    read_SLR_C30.py: reads C30 files from satellite laser ranging (CSR or GSFC)
    read_tellus_geocenter.py: reads degree 1 files from JPL Tellus
    read_swenson_geocenter.py: reads degree 1 files from Sean Swenson
    read_SLR_geocenter.py: reads degree 1 files from Satellite Laser Ranging
    read_GRACE_geocenter.py: reads degree 1 files from Sutterley et al. (2019)
    read_GRACE_harmonics.py: read spherical harmonic data from SHM files
    read_gfc_harmonics.py: reads spherical harmonic data from gfc files

UPDATE HISTORY:
    Updated 09/2021: added time-variable gravity data from GRAZ and Swarm
    Updated 08/2021: fix spherical harmonic errors for SLR C21,S21,C22,S22
    Updated 07/2021: fix inputs to AOD-corrected SLR geocenter coefficients
        output uncalibrated spherical harmonic errors (eclm and eslm)
        test if GRACE directory for product exists at start of program
    Updated 06/2021: can use SLR figure axis harmonics produced by GSFC
        read GRACE/GRACE-FO fields before reading replacement values
    Updated 05/2021: can use SLR low-degree harmonic values produced by GFZ
        define int/float precision to prevent deprecation warning
    Updated 04/2021: can replace figure axis and azimuthal dependence with SLR
    Updated 12/2020: updated SLR geocenter for new solutions from Minkang Cheng
    Updated 11/2020: set regress_model RELATIVE option to 2003.3 to match others
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 06/2020: set relative time to mean of input within regress_model
    Updated 03/2020: for public release.  output degree and order in dict
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: edit exception for degree 1 if no matching month is found
    Updated 08/2019: added LARES C30 from John Ries (C30_LARES_filtered.txt)
        added exceptions for C20, C30 and degree 1 if no matching month is found
    Updated 07/2019: added option for reading SLR C30 data (TN-14)
        added option DEG1_GIA to set the GIA used when calculating geocenter
    Updated 11/2018: split atmospheric jump corrections into a separate function
    Updated 10/2018: decode gzip read for python3 Compatibility
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 05/2018: set filename for release 6 C20 coefficients from SLR
    Updated 03/2018: include a flag if using ECMWF "jump" corrections
    Updated 12/2015: added ECMWF GAG "jump" correction. New TN-09 GAF correction
    Updated 09/2015: added sea level fingerprint geocenter option
    Updated 08/2015: added ECMWF "jump" corrections GAE and GAF
        as described in the AOD1B processing document (pages 21-26)
        ftp://podaac.jpl.nasa.gov/allData/grace/docs/AOD1B_20140520.pdf
        Made pole tide correction optional as POLE_TIDE (Y/N)
    Updated 05/2015: added Pole Tide correction from Wahr et al. (2015):
            The Pole Tide and its Effect on GRACE Time-Variable Gravity
            Measurements: Implications for Estimates of Surface Mass Variations
        Added MMAX option for new data files with a lower truncation
    Updated 12/2014: added portion for using iterated geocenter coefficients
    Updated 05/2014: can remove different means from an input file
    Updated 02/2014: minor update to if statements.
    Updated 01/2014: output directory of exact GRACE product
    Updated 09/2013: added option to least-squares model the degree 1
        coefficients for missing months
    Updated 07/2013: can use different geocenter solutions
    Written 05/2013
"""
from __future__ import print_function, division

import os
import re
import gzip
import copy
import numpy as np
from gravity_toolkit.grace_date import grace_date
from gravity_toolkit.read_SLR_C20 import read_SLR_C20
from gravity_toolkit.read_SLR_CS2 import read_SLR_CS2
from gravity_toolkit.read_SLR_C30 import read_SLR_C30
from gravity_toolkit.read_SLR_C50 import read_SLR_C50
from gravity_toolkit.read_swenson_geocenter import read_swenson_geocenter
from gravity_toolkit.read_tellus_geocenter import read_tellus_geocenter
from gravity_toolkit.read_SLR_geocenter import aod_corrected_SLR_geocenter
from gravity_toolkit.read_gravis_geocenter import read_gravis_geocenter
from read_GRACE_geocenter.read_GRACE_geocenter import read_GRACE_geocenter
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics

def grace_input_months(base_dir, PROC, DREL, DSET, LMAX, start_mon, end_mon,
    missing, SLR_C20, DEG1, MMAX=None, SLR_21='', SLR_22='',  SLR_C30='',
    SLR_C50='', MODEL_DEG1=False, DEG1_GIA='', ATM=False, POLE_TIDE=False):
    """
    Reads GRACE/GRACE-FO files for a spherical harmonic degree and order
        and a date range
    Can include geocenter values for degree 1 coefficients
    Can replace C20 with SLR values for all months
    Can replace C21,S21,C22,S22,C30,C50 with SLR values for months 179+
    Can correct for ECMWF atmospheric "jumps" using GAE/GAF/GAG files
    Can correct for Pole Tide drift following Wahr et al. (2015)

    Arguments
    ---------
    base_dir: Working data directory for GRACE/GRACE-FO data
    PROC: (CSR/CNES/JPL/GFZ) data processing center
    DREL: (RL01/RL02/RL03/RL04/RL05/RL06) data release
    DSET: (GAA/GAB/GAC/GAD/GSM) data product
    LMAX: Upper bound of Spherical Harmonic Degrees
    start_mon: starting month to consider in analysis
    end_mon: ending month to consider in analysis
    missing: missing months to not consider in analysis
    SLR_C20: Replaces C20 with SLR values
        N: use original values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GFZ: use values from GFZ
        GSFC: use values from GSFC (TN-14)
    DEG1: Use Degree 1 coefficients
        None: No degree 1
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        SLR: satellite laser ranging coefficients from CSR
        SLF: Sutterley and Velicogna coefficients, Remote Sensing (2019)
        Swenson: GRACE-derived coefficients from Sean Swenson
        GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS

    Keyword arguments
    -----------------
    MMAX: Upper bound of Spherical Harmonic Orders
    SLR_21: replaces C21 and S21 with SLR values
        None: use original values
        CSR: use values from CSR
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC
    SLR_22: replaces C22 and S22 with SLR values
        None: use original values
        CSR: use values from CSR
    SLR_C30: replaces C30 with SLR values
        None: use original values
        CSR: use values from CSR (5x5 with 6,1)
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC (TN-14)
    SLR_C50: replaces C50 with SLR values
        None: use original values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    POLE_TIDE: correct GSM data with pole tides following Wahr et al (2015)
    ATM: correct data with ECMWF "jump" corrections GAE, GAF and GAG
    MODEL_DEG1: least-squares model missing degree 1 coefficients (True/False)
    DEG1_GIA: GIA-correction used when calculating degree 1 coefficients

    Returns
    -------
    clm: GRACE/GRACE-FO cosine spherical harmonics
    slm: GRACE/GRACE-FO sine spherical harmonics
    eclm: GRACE/GRACE-FO uncalibrated cosine spherical harmonic errors
    eslm: GRACE/GRACE-FO uncalibrated sine spherical harmonic errors
    time: time of each GRACE/GRACE-FO measurement (mid-month)
    month: GRACE/GRACE-FO months of input datasets
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    title: string denoting low degree zonals replacement, geocenter usage and corrections
    directory: directory of exact GRACE/GRACE-FO product
    """

    #-- Directory of exact GRACE product
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)
    #-- test if GRACE product directory exists
    if not os.access(grace_dir, os.F_OK):
        raise FileNotFoundError(grace_dir)

    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX

    #-- Range of months from start_mon to end_mon (end_mon+1 to include end_mon)
    #-- Removing the missing months and months not to consider
    months = sorted(set(np.arange(start_mon,end_mon+1)) - set(missing))
    #-- number of months to consider in analysis
    n_cons = len(months)

    #-- Initializing input data matrices
    grace_Ylms = {}
    grace_Ylms['clm'] = np.zeros((LMAX+1,MMAX+1,n_cons))
    grace_Ylms['slm'] = np.zeros((LMAX+1,MMAX+1,n_cons))
    grace_Ylms['eclm'] = np.zeros((LMAX+1,MMAX+1,n_cons))
    grace_Ylms['eslm'] = np.zeros((LMAX+1,MMAX+1,n_cons))
    grace_Ylms['time'] = np.zeros((n_cons))
    grace_Ylms['month'] = np.zeros((n_cons),dtype=np.int64)
    #-- output dimensions
    grace_Ylms['l'] = np.arange(LMAX+1)
    grace_Ylms['m'] = np.arange(MMAX+1)
    #-- input directory for product
    grace_Ylms['directory'] = copy.copy(grace_dir)

    #-- associate GRACE/GRACE-FO files with each GRACE/GRACE-FO month
    grace_files=grace_date(base_dir,PROC=PROC,DREL=DREL,DSET=DSET,OUTPUT=False)

    #-- importing data from GRACE/GRACE-FO files
    for i,grace_month in enumerate(months):
        #-- read spherical harmonic data products
        infile = grace_files[grace_month]
        if PROC in ('GRAZ','Swarm'):
            #-- Degree 2 zonals will be converted to a tide free state
            Ylms = read_gfc_harmonics(infile, TIDE='tide_free')
        else:
            #-- Effects of Pole tide drift will be compensated if specified
            Ylms = read_GRACE_harmonics(infile, LMAX, MMAX=MMAX,
                POLE_TIDE=POLE_TIDE)
        #-- truncate harmonics to degree and order
        grace_Ylms['clm'][:,:,i] = Ylms['clm'][0:LMAX+1,0:MMAX+1]
        grace_Ylms['slm'][:,:,i] = Ylms['slm'][0:LMAX+1,0:MMAX+1]
        #-- truncate harmonic errors to degree and order
        grace_Ylms['eclm'][:,:,i] = Ylms['eclm'][0:LMAX+1,0:MMAX+1]
        grace_Ylms['eslm'][:,:,i] = Ylms['eslm'][0:LMAX+1,0:MMAX+1]
        #-- copy date variables
        grace_Ylms['time'][i] = np.copy(Ylms['time'])
        grace_Ylms['month'][i] = np.int64(grace_month)
    #-- single accelerometer months
    single_acc_months = np.copy(grace_Ylms['month'][grace_Ylms['month'] > 176])

    #-- SLR low-degree harmonic, geocenter and correction flags
    FLAGS = []

    #-- Replacing C2,0 with SLR C2,0
    #-- Running function read_SLR_C20.py
    #-- reading SLR C2,0 file for given release if specified
    if (SLR_C20 == 'CSR'):
        if (DREL == 'RL04'):
            SLR_file = os.path.join(base_dir,'TN-05_C20_SLR.txt')
        elif (DREL == 'RL05'):
            SLR_file = os.path.join(base_dir,'TN-07_C20_SLR.txt')
        elif (DREL == 'RL06'):
            # SLR_file = os.path.join(base_dir,'TN-11_C20_SLR.txt')
            SLR_file = os.path.join(base_dir,'C20_RL06.txt')
        C20_input = read_SLR_C20(SLR_file)
        FLAGS.append('_wCSR_C20')
    elif (SLR_C20 == 'GFZ'):
        SLR_file=os.path.join(base_dir,'GFZ_{0}_C20_SLR.dat'.format(DREL))
        C20_input = read_SLR_C20(SLR_file)
        FLAGS.append('_wGFZ_C20')
    elif (SLR_C20 == 'GSFC'):
        SLR_file=os.path.join(base_dir,'TN-14_C30_C20_GSFC_SLR.txt')
        C20_input = read_SLR_C20(SLR_file)
        FLAGS.append('_wGSFC_C20')

    #-- Replacing C2,1/S2,1 with SLR
    #-- Running function read_SLR_CS2.py
    if (SLR_21 == 'CSR'):
        SLR_file = os.path.join(base_dir,'C21_S21_{0}.txt'.format(DREL))
        C21_input = read_SLR_CS2(SLR_file)
        FLAGS.append('_wCSR_21')
    elif (SLR_21 == 'GFZ'):
        GravIS_file = 'GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat'
        SLR_file = os.path.join(base_dir,GravIS_file)
        C21_input = read_SLR_CS2(SLR_file)
        FLAGS.append('_wGFZ_21')
    elif (SLR_21 == 'GSFC'):
        #-- calculate monthly averages from 7-day arcs
        # SLR_file = os.path.join(base_dir,'GSFC_C21_S21.txt')
        SLR_file = os.path.join(base_dir,'gsfc_slr_5x5c61s61.txt')
        C21_input = read_SLR_CS2(SLR_file, DATE=grace_Ylms['time'], ORDER=1)
        FLAGS.append('_wGSFC_21')

    #-- Replacing C2,2/S2,2 with SLR
    #-- Running function read_SLR_CS2.py
    if (SLR_22 == 'CSR'):
        SLR_file = os.path.join(base_dir,'C22_S22_{0}.txt'.format(DREL))
        C22_input = read_SLR_CS2(SLR_file)
        FLAGS.append('_wCSR_22')
    elif (SLR_22 == 'GSFC'):
        SLR_file = os.path.join(base_dir,'gsfc_slr_5x5c61s61.txt')
        C22_input = read_SLR_CS2(SLR_file, DATE=grace_Ylms['time'], ORDER=2)
        FLAGS.append('_wGSFC_22')

    #-- Replacing C3,0 with SLR C3,0
    #-- Running function read_SLR_C30.py
    if (SLR_C30 == 'CSR'):
        SLR_file=os.path.join(base_dir,'CSR_Monthly_5x5_Gravity_Harmonics.txt')
        C30_input = read_SLR_C30(SLR_file)
        FLAGS.append('_wCSR_C30')
    elif (SLR_C30 == 'LARES'):
        SLR_file=os.path.join(base_dir,'C30_LARES_filtered.txt')
        C30_input = read_SLR_C30(SLR_file)
        FLAGS.append('_wLARES_C30')
    elif (SLR_C30 == 'GSFC'):
        SLR_file=os.path.join(base_dir,'TN-14_C30_C20_GSFC_SLR.txt')
        C30_input = read_SLR_C30(SLR_file)
        FLAGS.append('_wGSFC_C30')
    elif (SLR_C30 == 'GFZ'):
        GravIS_file = 'GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat'
        SLR_file = os.path.join(base_dir,GravIS_file)
        C30_input = read_SLR_C30(SLR_file)
        FLAGS.append('_wGFZ_C30')

    #-- Replacing C5,0 with SLR C5,0
    #-- Running function read_SLR_C50.py
    if (SLR_C50 == 'CSR'):
        SLR_file=os.path.join(base_dir,'CSR_Monthly_5x5_Gravity_Harmonics.txt')
        C50_input = read_SLR_C50(SLR_file)
        FLAGS.append('_wCSR_C50')
    elif (SLR_C50 == 'LARES'):
        SLR_file=os.path.join(base_dir,'C50_LARES_filtered.txt')
        C50_input = read_SLR_C50(SLR_file)
        FLAGS.append('_wLARES_C50')
    elif (SLR_C50 == 'GSFC'):
        # SLR_file=os.path.join(base_dir,'GSFC_SLR_C20_C30_C50_GSM_replacement.txt')
        SLR_file = os.path.join(base_dir,'gsfc_slr_5x5c61s61.txt')
        C50_input = read_SLR_C50(SLR_file, DATE=grace_Ylms['time'])
        FLAGS.append('_wGSFC_C50')

    #-- Correcting for Degree 1 (geocenter variations)
    #-- reading degree 1 file for given release if specified
    if (DEG1 == 'Tellus'):
        #-- Tellus (PO.DAAC) degree 1
        if DREL in ('RL04','RL05'):
            DEG1_file = os.path.join(base_dir,'geocenter',
                'deg1_coef_{0}.txt'.format(DREL))
            JPL = False
        else:
            DEG1_file = os.path.join(base_dir,'geocenter',
                'TN-13_GEOC_{0}_{1}.txt'.format(PROC,DREL))
            JPL = True
        #-- Running function read_tellus_geocenter.py
        DEG1_input = read_tellus_geocenter(DEG1_file,JPL=JPL)
        FLAGS.append('_w{0}_DEG1'.format(DEG1))
    elif (DEG1 == 'SLR'):
        #-- CSR Satellite Laser Ranging (SLR) degree 1
        # #-- SLR-derived degree-1 mass variations
        # #-- ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
        # DEG1_file=os.path.join(base_dir,'geocenter','GCN_{0}.txt'.format(DREL))
        # DEG1_input=aod_corrected_SLR_geocenter(DEG1_file,DREL,HEADER=16,
        #     COLUMNS=['time','X','Y','Z','X_sigma','Y_sigma','Z_sigma'])

        # #-- new CF-CM file of degree-1 mass variations
        # #-- https://cddis.nasa.gov/lw20/docs/2016/papers/14-Ries_paper.pdf
        # #-- http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
        # DEG1_file = os.path.join(base_dir,'geocenter','GCN_L1_L2_30d_CF-CM.txt')
        # DEG1_input = aod_corrected_SLR_geocenter(DEG1_file,DREL,HEADER=111,
        #     COLUMNS=['time','X','Y','Z','X_sigma','Y_sigma','Z_sigma'])

        #-- new file of degree-1 mass variations from Minkang Cheng
        #-- http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s
        DEG1_file = os.path.join(base_dir,'geocenter','gct2est.220_5s')
        DEG1_input = aod_corrected_SLR_geocenter(DEG1_file,DREL,HEADER=15,
            RADIUS=6.378136e9,COLUMNS=['MJD','time','X','Y','Z','XM','YM','ZM',
            'X_sigma','Y_sigma','Z_sigma','XM_sigma','YM_sigma','ZM_sigma'])
        FLAGS.append('_w{0}_DEG1'.format(DEG1))
    elif (DEG1 == 'SLF'):
        #-- read iterated degree one files from Sutterley and Velicogna (2019)
        #-- that includes self-attraction and loading effects
        #-- include flag for datasets with different GIA corrections
        MODEL = dict(RL04='OMCT', RL05='OMCT', RL06='MPIOM')
        args = (PROC,DREL,MODEL[DREL],'SLF_iter',DEG1_GIA)
        DEG1_file = os.path.join(base_dir,'geocenter',
            '{0}_{1}_{2}_{3}{4}.txt'.format(*args))
        DEG1_input = read_GRACE_geocenter(DEG1_file)
        FLAGS.append('_w{0}_DEG1'.format(DEG1))
    elif (DEG1 == 'Swenson'):
        #-- degree 1 coefficients provided by Sean Swenson in mm w.e.
        DEG1_file = os.path.join(base_dir,'geocenter',
            'gad_gsm.{0}.txt'.format(DREL))
        #-- Running function read_swenson_geocenter.py
        DEG1_input = read_swenson_geocenter(DEG1_file)
        FLAGS.append('_w{0}_DEG1'.format(DEG1))
    elif (DEG1 == 'GFZ'):
        #-- degree 1 coefficients provided by GFZ GravIS
        #-- http://gravis.gfz-potsdam.de/corrections
        DEG1_file = os.path.join(base_dir,'geocenter',
            'GRAVIS-2B_GFZOP_GEOCENTER_0002.dat')
        #-- Running function read_gravis_geocenter.py
        DEG1_input = read_gravis_geocenter(DEG1_file)
        FLAGS.append('_w{0}_DEG1'.format(DEG1))

    #-- atmospheric flag if correcting ECMWF "jumps" (using GAE/GAF/GAG files)
    if ATM:
        FLAGS.append('_wATM')
    #-- pole tide flag if correcting for pole tide drift (Wahr et al. 2015)
    if POLE_TIDE:
        FLAGS.append('_wPT')
    #-- full output string (SLR, geocenter and correction flags)
    grace_Ylms['title'] = ''.join(FLAGS)

    #-- Replace C20 with SLR coefficients
    if SLR_C20 in ('CSR','GFZ','GSFC'):
        #-- verify that there are replacement C20 months for specified range
        months_test = sorted(set(months) - set(C20_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C20 Months ({0})'.format(gm))
        #-- replace C20 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C20_input['month'] == grace_month)
            if (count != 0):
                k, = np.nonzero(C20_input['month'] == grace_month)
                grace_Ylms['clm'][2,0,i] = np.copy(C20_input['data'][k])
                grace_Ylms['eclm'][2,0,i] = np.copy(C20_input['error'][k])

    #-- Replace C21/S21 with SLR coefficients for single-accelerometer months
    if SLR_21 in ('CSR','GFZ','GSFC'):
        #-- verify that there are replacement C21/S21 months for specified range
        months_test = sorted(set(single_acc_months) - set(C21_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C21/S21 Months ({0})'.format(gm))
        #-- replace C21/S21 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C21_input['month'] == grace_month)
            if (count != 0) and (grace_month > 176):
                k, = np.nonzero(C21_input['month'] == grace_month)
                grace_Ylms['clm'][2,1,i] = np.copy(C21_input['C2m'][k])
                grace_Ylms['slm'][2,1,i] = np.copy(C21_input['S2m'][k])
                grace_Ylms['eclm'][2,1,i] = np.copy(C21_input['eC2m'][k])
                grace_Ylms['eslm'][2,1,i] = np.copy(C21_input['eS2m'][k])

    #-- Replace C22/S22 with SLR coefficients for single-accelerometer months
    if SLR_22 in ('CSR','GSFC'):
        #-- verify that there are replacement C22/S22 months for specified range
        months_test = sorted(set(single_acc_months) - set(C22_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C22/S22 Months ({0})'.format(gm))
        #-- replace C22/S22 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C22_input['month'] == grace_month)
            if (count != 0) and (grace_month > 176):
                k, = np.nonzero(C22_input['month'] == grace_month)
                grace_Ylms['clm'][2,2,i] = np.copy(C22_input['C2m'][k])
                grace_Ylms['slm'][2,2,i] = np.copy(C22_input['S2m'][k])
                grace_Ylms['eclm'][2,2,i] = np.copy(C22_input['eC2m'][k])
                grace_Ylms['eslm'][2,2,i] = np.copy(C22_input['eS2m'][k])

    #-- Replace C30 with SLR coefficients for single-accelerometer months
    if SLR_C30 in ('CSR','GFZ','GSFC','LARES'):
        #-- verify that there are replacement C30 months for specified range
        months_test = sorted(set(single_acc_months) - set(C30_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C30 Months ({0})'.format(gm))
        #-- replace C30 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C30_input['month'] == grace_month)
            if (count != 0) and (grace_month > 176):
                k, = np.nonzero(C30_input['month'] == grace_month)
                grace_Ylms['clm'][3,0,i] = np.copy(C30_input['data'][k])
                grace_Ylms['eclm'][3,0,i] = np.copy(C30_input['error'][k])

    #-- Replace C50 with SLR coefficients for single-accelerometer months
    if SLR_C50 in ('CSR','GSFC','LARES'):
        #-- verify that there are replacement C50 months for specified range
        months_test = sorted(set(single_acc_months) - set(C50_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C50 Months ({0})'.format(gm))
        #-- replace C50 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C50_input['month'] == grace_month)
            if (count != 0) and (grace_month > 176):
                k, = np.nonzero(C50_input['month'] == grace_month)
                grace_Ylms['clm'][5,0,i] = np.copy(C50_input['data'][k])
                grace_Ylms['eclm'][5,0,i] = np.copy(C50_input['error'][k])

    #-- Use Degree 1 coefficients
    #-- Tellus: Tellus Degree 1 (PO.DAAC following Sun et al., 2016)
    #-- SLR: CSR Satellite Laser Ranging (SLR) Degree 1 - GRACE AOD
    #-- SLF: OMCT/MPIOM coefficients with Sea Level Fingerprint land-water mass
    #-- Swenson: GRACE-derived coefficients from Sean Swenson
    #-- GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    if DEG1 in ('Tellus','SLR','SLF','Swenson','GFZ'):
        #-- check if modeling degree 1 or if all months are available
        if MODEL_DEG1:
            #-- least-squares modeling the degree 1 coefficients
            #-- fitting annual, semi-annual, linear and quadratic terms
            C10_model = regress_model(DEG1_input['time'], DEG1_input['C10'],
                grace_Ylms['time'], ORDER=2, CYCLES=[0.5,1.0], RELATIVE=2003.3)
            C11_model = regress_model(DEG1_input['time'], DEG1_input['C11'],
                grace_Ylms['time'], ORDER=2, CYCLES=[0.5,1.0], RELATIVE=2003.3)
            S11_model = regress_model(DEG1_input['time'], DEG1_input['S11'],
                grace_Ylms['time'], ORDER=2, CYCLES=[0.5,1.0], RELATIVE=2003.3)
        else:
            #-- check that all months are available for a given geocenter
            months_test = sorted(set(months) - set(DEG1_input['month']))
            if months_test:
                gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
                raise IOError('No Matching Geocenter Months ({0})'.format(gm))
        #-- for each considered date
        for i,grace_month in enumerate(months):
            k, = np.nonzero(DEG1_input['month'] == grace_month)
            count = np.count_nonzero(DEG1_input['month'] == grace_month)
            #-- Degree 1 is missing for particular month
            if (count == 0) and MODEL_DEG1:
                #-- using least-squares modeled coefficients from
                #-- lsq_model_degree_one.py
                grace_Ylms['clm'][1,0,i] = np.copy(C10_model[i])
                grace_Ylms['clm'][1,1,i] = np.copy(C11_model[i])
                grace_Ylms['slm'][1,1,i] = np.copy(S11_model[i])
            else:#-- using coefficients from data file
                grace_Ylms['clm'][1,0,i] = np.copy(DEG1_input['C10'][k])
                grace_Ylms['clm'][1,1,i] = np.copy(DEG1_input['C11'][k])
                grace_Ylms['slm'][1,1,i] = np.copy(DEG1_input['S11'][k])

    #-- read and add/remove the GAE and GAF atmospheric correction coefficients
    if ATM:
        #-- read ECMWF correction files from Fagiolini et al. (2015)
        atm_corr = read_ecmwf_corrections(base_dir,LMAX,months,MMAX=MMAX)
        #-- Removing GAE/GAF/GAG from RL05 GSM Products
        if (DSET == 'GSM'):
            for m in range(0,MMAX+1):#-- MMAX+1 to include l
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    grace_Ylms['clm'][l,m,:] -= atm_corr['clm'][l,m,:]
                    grace_Ylms['slm'][l,m,:] -= atm_corr['slm'][l,m,:]
        #-- Adding GAE/GAF/GAG to RL05 Atmospheric Products (GAA,GAC)
        elif DSET in ('GAC','GAA'):
            for m in range(0,MMAX+1):#-- MMAX+1 to include l
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    grace_Ylms['clm'][l,m,:] += atm_corr['clm'][l,m,:]
                    grace_Ylms['slm'][l,m,:] += atm_corr['slm'][l,m,:]

    #-- return the harmonic solutions with possible low-degree replacements
    #-- return the harmonic dimensions (spectral and temporal)
    #-- return string specifying processing and correction flags
    #-- return directory of exact GRACE/GRACE-FO product
    return grace_Ylms

#-- PURPOSE: read atmospheric jump corrections from Fagiolini et al. (2015)
def read_ecmwf_corrections(base_dir, LMAX, months, MMAX=None):
    """
    Read atmospheric jump corrections from Fagiolini et al. (2015)

    Arguments
    ---------
    base_dir: Working data directory for GRACE/GRACE-FO data
    LMAX: Upper bound of Spherical Harmonic Degrees
    months: list of GRACE/GRACE-FO months

    Keyword arguments
    -----------------
    MMAX: Upper bound of Spherical Harmonic orders

    Returns
    -------
    clm: atmospheric correction cosine spherical harmonics
    slm: atmospheric correction sine spherical harmonics
    """
    #-- correction files
    corr_file = {}
    corr_file['GAE'] = 'TN-08_GAE-2_2006032-2010031_0000_EIGEN_G---_0005.gz'
    corr_file['GAF'] = 'TN-09_GAF-2_2010032-2015131_0000_EIGEN_G---_0005.gz'
    corr_file['GAG'] = 'TN-10_GAG-2_2015132-2099001_0000_EIGEN_G---_0005.gz'
    #-- atmospheric correction coefficients
    atm_corr_clm = {}
    atm_corr_slm = {}
    #-- number of months to consider in analysis
    n_cons = len(months)
    #-- set maximum order if not equal to maximum degree
    MMAX = LMAX if (MMAX is None) else MMAX
    #-- iterate through python dictionary keys (GAE, GAF, GAG)
    for key, val in corr_file.items():
        #-- allocate for clm and slm of atmospheric corrections
        atm_corr_clm[key] = np.zeros((LMAX+1,MMAX+1))
        atm_corr_slm[key] = np.zeros((LMAX+1,MMAX+1))
        #-- GRACE correction files are compressed gz files
        with gzip.open(os.path.join(base_dir, val),'rb') as f:
            file_contents = f.read().decode('ISO-8859-1').splitlines()
        #-- for each line in the GRACE correction file
        for line in file_contents:
            #-- find if line starts with GRCOF2
            if bool(re.match(r'GRCOF2',line)):
                #-- split the line into individual components
                line_contents = line.split()
                #-- degree and order for the line
                l1 = np.int64(line_contents[1])
                m1 = np.int64(line_contents[2])
                #-- if degree and order are below the truncation limits
                if ((l1 <= LMAX) and (m1 <= MMAX)):
                    atm_corr_clm[key][l1,m1] = np.float64(line_contents[3])
                    atm_corr_slm[key][l1,m1] = np.float64(line_contents[4])

    #-- create output atmospheric corrections to be removed/added to data
    atm_corr = {}
    atm_corr['clm'] = np.zeros((LMAX+1,LMAX+1,n_cons))
    atm_corr['slm'] = np.zeros((LMAX+1,LMAX+1,n_cons))
    #-- for each considered date
    for i,grace_month in enumerate(months):
        #-- remove correction based on dates
        if (grace_month >= 50) & (grace_month <= 97):
            atm_corr['clm'][:,:,i] = atm_corr_clm['GAE'][:,:]
            atm_corr['slm'][:,:,i] = atm_corr_slm['GAE'][:,:]
        elif (grace_month >= 98) & (grace_month <= 161):
            atm_corr['clm'][:,:,i] = atm_corr_clm['GAF'][:,:]
            atm_corr['slm'][:,:,i] = atm_corr_slm['GAF'][:,:]
        elif (grace_month > 161):
            atm_corr['clm'][:,:,i] = atm_corr_clm['GAG'][:,:]
            atm_corr['slm'][:,:,i] = atm_corr_slm['GAG'][:,:]

    #-- return the atmospheric corrections
    return atm_corr

#-- PURPOSE: calculate a regression model for extrapolating values
def regress_model(t_in, d_in, t_out, ORDER=2, CYCLES=None, RELATIVE=0.0):
    """
    Calculates a regression model for extrapolating values

    Arguments
    ---------
    t_in: input time array
    d_in: input data array
    t_out: time array for output regressed values

    Keyword arguments
    -----------------
    ORDER: maximum polynomial order for regression model
    CYCLES: list of cyclical terms to include in regression model
    RELATIVE: relative time for polynomial coefficients in fit

    Returns
    -------
    d_out: output regressed value data array
    """

    #-- remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    t_out = np.squeeze(t_out)
    #-- check dimensions of output
    t_out = np.atleast_1d(t_out)
    #-- set relative to mean of input time
    if not RELATIVE:
        RELATIVE = np.mean(t_in)
    #-- create design matrix based on polynomial order and harmonics
    DMAT = []
    MMAT = []
    #-- add polynomial orders (0=constant, 1=linear, 2=quadratic)
    for o in range(ORDER+1):
        DMAT.append((t_in-RELATIVE)**o)
        MMAT.append((t_out-RELATIVE)**o)
    #-- add cyclical terms (0.5=semi-annual, 1=annual)
    for c in CYCLES:
        DMAT.append(np.sin(2.0*np.pi*t_in/np.float64(c)))
        DMAT.append(np.cos(2.0*np.pi*t_in/np.float64(c)))
        MMAT.append(np.sin(2.0*np.pi*t_out/np.float64(c)))
        MMAT.append(np.cos(2.0*np.pi*t_out/np.float64(c)))
    #-- Calculating Least-Squares Coefficients
    #-- Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(np.transpose(DMAT), d_in, rcond=-1)[0]
    #-- return modeled time-series
    return np.dot(np.transpose(MMAT),beta_mat)
