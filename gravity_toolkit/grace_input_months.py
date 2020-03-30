#!/usr/bin/env python
u"""
grace_input_months.py
Written by Tyler Sutterley (03/2020)

Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order
    and for a specified date range
Replaces Degree 1 with with input values (if specified)
Replaces C20 with SLR values (if specified)
Replaces C30 with SLR values for months 179+ (if specified)
Corrects for ECMWF atmospheric "jumps" using the GAE, GAF and GAG files
Corrects for Pole Tide drift following Wahr et al. (2015)

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data
    PROC: (CSR/CNES/JPL/GFZ) data processing center
    DREL: (RL01,RL02,RL03,RL04,RL05,RL06) data release
    DSET: (GAA/GAB/GAC/GAD/GSM) data product
    LMAX: Upper bound of Spherical Harmonic Degrees (e.g. 60)
    start_mon: starting month to consider in analysis
    end_mon: ending month to consider in analysis
    missing: missing months to not consider in analysis
    SLR_C20: Replaces C20 with SLR values
        N: use original values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GSFC: use values from GSFC (TN-14)
    DEG1: Use Degree 1 coefficients
        None: No degree 1
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
            https://grace.jpl.nasa.gov/data/get-data/geocenter/
        SLR: satellite laser ranging coefficients from CSR
            ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
        SLF: Sutterley and Velicogna coefficients, Remote Sensing (2019)
            https://doi.org/10.6084/m9.figshare.7388540

OUTPUTS:
    clm: GRACE/GRACE-FO cosine spherical harmonic to degree/order LMAX and MMAX
    slm: GRACE/GRACE-FO sine spherical harmonic to degree/order LMAX and MMAX
    time: time of each GRACE/GRACE-FO measurement (mid-month)
    month: GRACE/GRACE-FO months of input datasets
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    title: string denoting low degree zonals, geocenter and corrections
    directory: directory of exact GRACE/GRACE-FO product

OPTIONS:
    MMAX: Upper bound of Spherical Harmonic Orders (default=LMAX)
    SLR_C30: replaces C30 with SLR values
        None: use original values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC (TN-14)
    POLE_TIDE: correct GSM data with pole tides following Wahr et al (2015)
    ATM: correct data with ECMWF "jump" corrections GAE, GAF and GAG
    MODEL_DEG1: least-squares model missing degree 1 coefficients (True/False)
    DEG1_GIA: GIA-correction used when calculating degree 1 coefficients

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)
    scipy: Scientific Tools for Python (http://www.scipy.org/)
    PyYAML: YAML parser and emitter for Python (https://github.com/yaml/pyyaml)

PROGRAM DEPENDENCIES:
    grace_date.py: reads GRACE index file and calculates dates for each month
    read_SLR_C20.py: reads C20 files from satellite laser ranging (CSR or GSFC)
    read_SLR_C30.py: reads C30 files from satellite laser ranging (CSR or GSFC)
    convert_julian.py: Return the calendar date and time given Julian date
    read_tellus_geocenter.py: reads PO.DAAC degree 1 files
    read_SLR_geocenter.py: reads degree 1 files from Satellite Laser Ranging
    read_GRACE_geocenter.py: reads degree 1 files from Sutterley et al. (2019)
    read_GRACE_harmonics.py: reads an input GRACE data file and calculates date

UPDATE HISTORY:
    Updated 03/2020 for public release.  output degree and order in dict
"""
from __future__ import print_function, division

import os
import re
import gzip
import numpy as np
from gravity_toolkit.grace_date import grace_date
from gravity_toolkit.read_SLR_C20 import read_SLR_C20
from gravity_toolkit.read_SLR_C30 import read_SLR_C30
from gravity_toolkit.read_tellus_geocenter import read_tellus_geocenter
from gravity_toolkit.read_SLR_geocenter import aod_corrected_SLR_geocenter
from read_GRACE_geocenter.read_GRACE_geocenter import read_GRACE_geocenter
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics

def grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
    start_mon, end_mon, missing, SLR_C20, DEG1, MMAX=None, SLR_C30='',
    MODEL_DEG1=False, DEG1_GIA='', ATM=False, POLE_TIDE=False):

    #-- Directory of exact GRACE product
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)

    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX

    #-- Replacing C2,0 with SLR C2,0
    #-- Running function read_SLR_C20.py
    #-- reading SLR C2,0 file for given release if specified
    if (SLR_C20 == 'CSR'):
        if (DREL == 'RL04'):
            SLR_file = os.path.join(base_dir,'TN-05_C20_SLR.txt')
        elif (DREL == 'RL05'):
            SLR_file = os.path.join(base_dir,'TN-07_C20_SLR.txt')
        elif (DREL == 'RL06'):
            SLR_file = os.path.join(base_dir,'TN-11_C20_SLR.txt')
        C20_input = read_SLR_C20(SLR_file)
        C20_str = '_wCSR_C20'
    elif (SLR_C20 == 'GSFC'):
        SLR_file=os.path.join(base_dir,'TN-14_C30_C20_GSFC_SLR.txt')
        C20_input = read_SLR_C20(SLR_file)
        C20_str = '_wGSFC_C20'
    else:
        C20_str = ''

    #-- Replacing C3,0 with SLR C3,0
    #-- Running function read_SLR_C30.py
    if (SLR_C30 == 'CSR'):
        SLR_file=os.path.join(base_dir,'CSR_Monthly_5x5_Gravity_Harmonics.txt')
        C30_input = read_SLR_C30(SLR_file)
        C30_str = '_wCSR_C30'
    elif (SLR_C30 == 'LARES'):
        SLR_file=os.path.join(base_dir,'C30_LARES_filtered.txt')
        C30_input = read_SLR_C30(SLR_file)
        C30_str = '_wLARES_C30'
    elif (SLR_C30 == 'GSFC'):
        SLR_file=os.path.join(base_dir,'TN-14_C30_C20_GSFC_SLR.txt')
        C30_input = read_SLR_C30(SLR_file)
        C30_str = '_wGSFC_C30'
    else:
        C30_str = ''

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
        DEG1_str = '_w{0}_DEG1'.format(DEG1)
    elif (DEG1 == 'SLR'):
        #-- CSR Satellite Laser Ranging (SLR) degree 1
        # #-- SLR-derived degree-1 mass variations
        # #-- ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
        # DEG1_file=os.path.join(base_dir,'geocenter','GCN_{0}.txt'.format(DREL))
        # DEG1_input=aod_corrected_slr_deg1(DEG1_file,DREL,skiprows=16)

        #-- new CF-CM file of degree-1 mass variations
        #-- https://cddis.nasa.gov/lw20/docs/2016/papers/14-Ries_paper.pdf
        #-- ftp://ftp.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
        DEG1_file = os.path.join(base_dir,'geocenter','GCN_L1_L2_30d_CF-CM.txt')
        DEG1_input = aod_corrected_slr_deg1(DEG1_file,DREL,skiprows=111)
        DEG1_str = '_w{0}_DEG1'.format(DEG1)
    elif (DEG1 == 'SLF'):
        #-- read iterated degree one files from Sutterley and Velicogna (2019)
        #-- that includes self-attraction and loading effects
        #-- include flag for datasets with different GIA corrections
        MODEL = dict(RL04='OMCT', RL05='OMCT', RL06='MPIOM')
        args = (PROC,DREL,MODEL[DREL],'SLF_iter',DEG1_GIA)
        DEG1_file = os.path.join(base_dir,'geocenter',
            '{0}_{1}_{2}_{3}{4}.txt'.format(*args))
        DEG1_input = read_GRACE_geocenter(DEG1_file)
        DEG1_str = '_w{0}_DEG1'.format(DEG1)
    else:#-- not using a degree 1 file (non-GSM or only using degree 2+)
        DEG1_str = ''

    #-- atmospheric flag if correcting ECMWF "jumps" (using GAE/GAF/GAG files)
    atm_str = '_wATM' if ATM else ''
    #-- pole tide flag if correcting for pole tide drift (Wahr et al. 2015)
    pt_str = '_wPT' if POLE_TIDE else ''
    #-- full output string (C20, C30, geocenter and atmospheric flags)
    out_str = C20_str + C30_str + DEG1_str + atm_str + pt_str

    #-- Range of months from start_mon to end_mon (end_mon+1 to include end_mon)
    #-- Removing the missing months and months not to consider
    months = sorted(set(np.arange(start_mon,end_mon+1)) - set(missing))
    #-- number of months to consider in analysis
    n_cons = len(months)

    #-- Initializing input data matrices
    grace_clm = np.zeros((LMAX+1,MMAX+1,n_cons))
    grace_slm = np.zeros((LMAX+1,MMAX+1,n_cons))
    tdec = np.zeros((n_cons))
    mon = np.zeros((n_cons),dtype=np.int)
    #-- output dimensions
    lout = np.arange(LMAX+1)
    mout = np.arange(MMAX+1)

    #-- associate GRACE/GRACE-FO files with each GRACE/GRACE-FO month
    grace_files=grace_date(base_dir,PROC=PROC,DREL=DREL,DSET=DSET,OUTPUT=False)

    #-- importing data from GRACE/GRACE-FO files
    for i,grace_month in enumerate(months):
        #-- Effects of Pole tide drift will be compensated if soecified
        infile = grace_files[grace_month]
        Ylms = read_GRACE_harmonics(infile,LMAX,MMAX=MMAX,POLE_TIDE=POLE_TIDE)
        grace_clm[:,:,i] = Ylms['clm'][0:LMAX+1,0:MMAX+1]
        grace_slm[:,:,i] = Ylms['slm'][0:LMAX+1,0:MMAX+1]
        tdec[i] = Ylms['time']
        mon[i] = np.int(grace_month)

    #-- Replace C20 with SLR coefficients
    if SLR_C20 in ('CSR','GSFC'):
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
                grace_clm[2,0,i] = C20_input['data'][k]

    #-- Replace C30 with SLR coefficients for single-accelerometer months
    if SLR_C30 in ('CSR','GSFC','LARES'):
        #-- verify that there are replacement C30 months for specified range
        months_test = sorted(set(mon[mon > 176]) - set(C30_input['month']))
        if months_test:
            gm = ','.join('{0:03d}'.format(gm) for gm in months_test)
            raise IOError('No Matching C30 Months ({0})'.format(gm))
        #-- replace C30 with SLR coefficients
        for i,grace_month in enumerate(months):
            count = np.count_nonzero(C30_input['month'] == grace_month)
            if (count != 0) and (grace_month > 176):
                k, = np.nonzero(C30_input['month'] == grace_month)
                grace_clm[3,0,i] = C30_input['data'][k]

    #-- Use Degree 1 coefficients
    #-- Tellus: Tellus Degree 1 (PO.DAAC following Sun et al., 2016)
    #-- SLR: CSR Satellite Laser Ranging (SLR) Degree 1 - GRACE AOD
    #-- SLF: OMCT/MPIOM coefficients with Sea Level Fingerprint land-water mass
    if DEG1 in ('Tellus','SLR','SLF'):
        #-- check if modeling degree 1 or if all months are available
        if MODEL_DEG1:
            #-- least-squares modeling the degree 1 coefficients
            #-- fitting annual, semi-annual, linear and quadratic terms
            C10_model = regress_model(DEG1_input['time'], DEG1_input['C10'],
                tdec, ORDER=2, CYCLES=[0.5,1.0])
            C11_model = regress_model(DEG1_input['time'], DEG1_input['C11'],
                tdec, ORDER=2, CYCLES=[0.5,1.0])
            S11_model = regress_model(DEG1_input['time'], DEG1_input['S11'],
                tdec, ORDER=2, CYCLES=[0.5,1.0])
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
                grace_clm[1,0,i] = C10_model[i]
                grace_clm[1,1,i] = C11_model[i]
                grace_slm[1,1,i] = S11_model[i]
            else:#-- using coefficients from data file
                grace_clm[1,0,i] = DEG1_input['C10'][k]
                grace_clm[1,1,i] = DEG1_input['C11'][k]
                grace_slm[1,1,i] = DEG1_input['S11'][k]

    #-- read and add/remove the GAE and GAF atmospheric correction coefficients
    if ATM:
        #-- read ECMWF correction files from Fagiolini et al. (2015)
        atm_corr = read_ecmwf_corrections(base_dir,LMAX,months,MMAX=MMAX)
        #-- Removing GAE/GAF/GAG from RL05 GSM Products
        if (DSET == 'GSM'):
            for m in range(0,MMAX+1):#-- MMAX+1 to include l
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    grace_clm[l,m,:] -= atm_corr['clm'][l,m,:]
                    grace_slm[l,m,:] -= atm_corr['slm'][l,m,:]
        #-- Adding GAE/GAF/GAG to RL05 Atmospheric Products (GAA,GAC)
        elif DSET in ('GAC','GAA'):
            for m in range(0,MMAX+1):#-- MMAX+1 to include l
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    grace_clm[l,m,:] += atm_corr['clm'][l,m,:]
                    grace_slm[l,m,:] += atm_corr['slm'][l,m,:]

    return {'clm':grace_clm, 'slm':grace_slm, 'time':tdec, 'month':mon,
        'l':lout, 'm':mout, 'title':out_str, 'directory':grace_dir}

#-- PURPOSE: read atmospheric jump corrections from Fagiolini et al. (2015)
def read_ecmwf_corrections(base_dir, LMAX, months, MMAX=None):
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
            if bool(re.match('GRCOF2',line)):
                #-- split the line into individual components
                line_contents = line.split()
                #-- degree and order for the line
                l1 = np.int(line_contents[1])
                m1 = np.int(line_contents[2])
                #-- if degree and order are below the truncation limits
                if ((l1 <= LMAX) and (m1 <= MMAX)):
                    atm_corr_clm[key][l1,m1] = np.float(line_contents[3])
                    atm_corr_slm[key][l1,m1] = np.float(line_contents[4])

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
def regress_model(t_in, d_in, t_out, ORDER=2, CYCLES=None, RELATIVE=None):

    #-- remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    t_out = np.squeeze(t_out)
    #-- check dimensions of output
    if (np.ndim(t_out) == 0):
        t_out = np.array([t_out])

    #-- CREATING DESIGN MATRIX FOR REGRESSION
    DMAT = []
    MMAT = []
    #-- add polynomial orders (0=constant, 1=linear, 2=quadratic)
    for o in range(ORDER+1):
        DMAT.append((t_in-RELATIVE)**o)
        MMAT.append((t_out-RELATIVE)**o)
    #-- add cyclical terms (0.5=semi-annual, 1=annual)
    for c in CYCLES:
        DMAT.append(np.sin(2.0*np.pi*t_in/np.float(c)))
        DMAT.append(np.cos(2.0*np.pi*t_in/np.float(c)))
        MMAT.append(np.sin(2.0*np.pi*t_out/np.float(c)))
        MMAT.append(np.cos(2.0*np.pi*t_out/np.float(c)))

    #-- Calculating Least-Squares Coefficients
    #-- Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(np.transpose(DMAT), d_in, rcond=-1)[0]

    #-- return modeled time-series
    return np.dot(np.transpose(MMAT),beta_mat)
