#!/usr/bin/env python
u"""
read_GRACE_harmonics.py
Written by Tyler Sutterley (05/2021)

Reads GRACE files and extracts spherical harmonic data and drift rates (RL04)
Adds drift rates to clm and slm for release 4 harmonics
Correct GSM data for drift in pole tide following Wahr et al. (2015)
Parses date of GRACE/GRACE-FO data from filename

INPUTS:
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    LMAX: Maximum degree of spherical harmonics (degree of truncation)

OPTIONS:
    MMAX: Maximum order of spherical harmonics (order of truncation)
        default is the maximum spherical harmonic degree
    POLE_TIDE: correct GSM data for pole tide drift following Wahr et al. (2015)

OUTPUTS:
    time: mid-month date in year-decimal
    start: start date of range as Julian day
    end: end date of range as Julian day
    clm: cosine spherical harmonics of input data (LMAX,MMAX)
    slm: sine spherical harmonics of input data (LMAX,MMAX)
    eclm: cosine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)
    eslm: sine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 12/2020: using utilities from time module
    Updated 08/2020: flake8 compatible regular expression strings
        input file can be "diskless" bytesIO object
    Updated 07/2020: added function docstrings
    Updated 08/2019: specify yaml loader (PyYAML yaml.load(input) Deprecation)
    Updated 07/2019: replace colons in yaml header if within quotations
    Updated 11/2018: decode gzip read with ISO-8859-1 for python3 compatibility
    Updated 05/2018: updates to file name structure with release 6 and GRACE-FO
        output file headers and parse new YAML headers for RL06 and GRACE-FO
    Written 10/2017 for public release
"""
import os
import re
import io
import gzip
import yaml
import numpy as np
import gravity_toolkit.time

#-- PURPOSE: read Level-2 GRACE and GRACE-FO spherical harmonic files
def read_GRACE_harmonics(input_file, LMAX, MMAX=None, POLE_TIDE=False):
    """
    Extracts spherical harmonic coefficients from GRACE/GRACE-FO files
    Adds drift rates to spherical harmonics for Release 4 data
    Correct data prior to Release 6 for pole tide drift
    Parses date of GRACE/GRACE-FO data from filename

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    LMAX: Maximum degree of spherical harmonics (degree of truncation)

    Keyword arguments
    -----------------
    MMAX: Maximum order of spherical harmonics
    POLE_TIDE: correct for pole tide drift following Wahr et al. (2015)

    Returns
    -------
    time: mid-month date in year-decimal
    start: start date of range as Julian day
    end: end date of range as Julian day
    clm: cosine spherical harmonics coefficients
    slm: sine spherical harmonics coefficients
    eclm: cosine spherical harmonic uncalibrated standard deviations
    eslm: sine spherical harmonic uncalibrated standard deviations
    """

    #-- parse filename
    PFX,SD,ED,N,PRC,F1,DRL,F2,SFX = parse_file(input_file)
    file_contents = extract_file(input_file, (SFX=='.gz'))

    #-- JPL Mascon solutions
    if PRC in ('JPLMSC'):
        DSET = 'GSM'
        DREL = np.int64(DRL)
        FLAG = r'GRCOF2'
    #-- Kusche et al. (2009) DDK filtered solutions 10.1007/s00190-009-0308-3
    elif PFX.startswith('kfilter_DDK'):
        DSET = 'GSM'
        DREL = np.int64(DRL)
        FLAG = r'gfc'
    #-- Standard GRACE solutions
    else:
        DSET = PFX
        DREL = np.int64(DRL)
        FLAG = r'GRCOF2'

    #-- output python dictionary with GRACE data and date information
    grace_L2_input = {}
    #-- extract GRACE date information from input file name
    start_yr = np.float64(SD[:4])
    end_yr = np.float64(ED[:4])
    start_day = np.float64(SD[4:])
    end_day = np.float64(ED[4:])
    #-- calculate mid-month date taking into account if measurements are
    #-- on different years
    dpy = gravity_toolkit.time.calendar_days(start_yr).sum()

    #-- For data that crosses years (end_yr - start_yr should be at most 1)
    end_cyclic = ((end_yr - start_yr)*dpy+end_day)
    #-- Calculate mid-month value
    mid_day = np.mean([start_day, end_cyclic])
    #-- Calculating the mid-month date in decimal form
    grace_L2_input['time'] = start_yr + mid_day/dpy
    #-- Calculating the Julian dates of the start and end date
    grace_L2_input['start'] = 2400000.5 + \
        gravity_toolkit.time.convert_calendar_dates(start_yr,1.0,start_day,
        epoch=(1858,11,17,0,0,0))
    grace_L2_input['end'] = 2400000.5 + \
        gravity_toolkit.time.convert_calendar_dates(end_yr,1.0,end_day,
        epoch=(1858,11,17,0,0,0))

    #-- set maximum spherical harmonic order
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX
    #-- Spherical harmonic coefficient matrices to be filled from data file
    grace_L2_input['clm'] = np.zeros((LMAX+1,MMAX+1))
    grace_L2_input['slm'] = np.zeros((LMAX+1,MMAX+1))
    #-- spherical harmonic uncalibrated standard deviations
    grace_L2_input['eclm'] = np.zeros((LMAX+1,MMAX+1))
    grace_L2_input['eslm'] = np.zeros((LMAX+1,MMAX+1))
    if ((DREL == 4) and (DSET == 'GSM')):
        #-- clm and slm drift rates for RL04
        drift_c = np.zeros((LMAX+1,MMAX+1))
        drift_s = np.zeros((LMAX+1,MMAX+1))

    #-- extract GRACE and GRACE-FO file headers
    #-- replace colons in header if within quotations
    head = [re.sub(r'\"(.*?)\:\s(.*?)\"',r'"\1, \2"',l) for l in file_contents
        if not re.match(r'{0}|GRDOTA'.format(FLAG),l)]
    if ((N == 'GRAC') and (DREL >= 6)) or (N == 'GRFO'):
        #-- parse the YAML header for RL06 or GRACE-FO (specifying yaml loader)
        grace_L2_input.update(yaml.load('\n'.join(head),Loader=yaml.BaseLoader))
    else:
        #-- save lines of the GRACE file header removing empty lines
        grace_L2_input['header'] = [l.rstrip() for l in head if l]

    #-- for each line in the GRACE/GRACE-FO file
    for line in file_contents:
        #-- find if line starts with data marker flag (e.g. GRCOF2)
        if bool(re.match(FLAG,line)):
            #-- split the line into individual components
            line_contents = line.split()
            #-- degree and order for the line
            l1 = np.int64(line_contents[1])
            m1 = np.int64(line_contents[2])
            #-- if degree and order are below the truncation limits
            if ((l1 <= LMAX) and (m1 <= MMAX)):
                grace_L2_input['clm'][l1,m1] = np.float64(line_contents[3])
                grace_L2_input['slm'][l1,m1] = np.float64(line_contents[4])
                grace_L2_input['eclm'][l1,m1] = np.float64(line_contents[5])
                grace_L2_input['eslm'][l1,m1] = np.float64(line_contents[6])
        #-- find if line starts with drift rate flag
        elif bool(re.match(r'GRDOTA',line)):
            #-- split the line into individual components
            line_contents = line.split()
            l1 = np.int64(line_contents[1])
            m1 = np.int64(line_contents[2])
            #-- Reading Drift rates for low degree harmonics
            drift_c[l1,m1] = np.float64(line_contents[3])
            drift_s[l1,m1] = np.float64(line_contents[4])

    #-- Adding drift rates to clm and slm for RL04
    #-- if drift rates exist at any time, will add to harmonics
    #-- Will convert the secular rates into a stokes contribution
    #-- Currently removes 2003.3 to get the temporal average close to 0.
    #-- note: += means grace_xlm = grace_xlm + drift_x
    if ((DREL == 4) and (DSET == 'GSM')):
        #-- time since 2003.3
        dt = (grace_L2_input['time']-2003.3)
        grace_L2_input['clm'][:,:] += dt*drift_c[:,:]
        grace_L2_input['slm'][:,:] += dt*drift_s[:,:]

    #-- Correct Pole Tide following Wahr et al. (2015) 10.1002/2015JB011986
    if POLE_TIDE and (DSET == 'GSM'):
        #-- time since 2000.0
        dt = (grace_L2_input['time']-2000.0)
        #-- CSR and JPL Pole Tide Correction
        if PRC in ('UTCSR','JPLEM','JPLMSC'):
            #-- values for IERS mean pole [2010]
            if (grace_L2_input['time'] < 2010.0):
                a = np.array([0.055974,1.8243e-3,1.8413e-4,7.024e-6])
                b = np.array([-0.346346,-1.7896e-3,1.0729e-4,0.908e-6])
            elif (grace_L2_input['time'] >= 2010.0):
                a = np.array([0.023513,7.6141e-3,0.0,0.0])
                b = np.array([-0.358891,0.6287e-3,0.0,0.0])
            #-- calculate m1 and m2 values
            m1 = np.copy(a[0])
            m2 = np.copy(b[0])
            for x in range(1,4):
                m1 += a[x]*dt**x
                m2 += b[x]*dt**x
            #-- pole tide values for CSR and JPL
            #-- CSR and JPL both remove the IERS mean pole from m1 and m2
            #-- before computing their harmonic solutions
            C21_PT = -1.551e-9*(m1 - 0.62e-3*dt) - 0.012e-9*(m2 + 3.48e-3*dt)
            S21_PT = 0.021e-9*(m1 - 0.62e-3*dt) - 1.505e-9*(m2 + 3.48e-3*dt)
            #-- correct GRACE spherical harmonics for pole tide
            #-- note: -= means grace_xlm = grace_xlm - PT
            grace_L2_input['clm'][2,1] -= C21_PT
            grace_L2_input['slm'][2,1] -= S21_PT
        #-- GFZ Pole Tide Correction
        elif PRC in ('EIGEN','GFZOP'):
            #-- pole tide values for GFZ
            #-- GFZ removes only a constant pole position
            C21_PT = -1.551e-9*(-0.62e-3*dt) - 0.012e-9*(3.48e-3*dt)
            S21_PT = 0.021e-9*(-0.62e-3*dt) - 1.505e-9*(3.48e-3*dt)
            #-- correct GRACE spherical harmonics for pole tide
            #-- note: -= means grace_xlm = grace_xlm - PT
            grace_L2_input['clm'][2,1] -= C21_PT
            grace_L2_input['slm'][2,1] -= S21_PT

    #-- return the GRACE data, GRACE date (mid-month in decimal), and the
    #-- start and end days as Julian dates
    return grace_L2_input

#-- PURPOSE: extract parameters from filename
def parse_file(input_file):
    """
    Extract parameters from filename

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    """
    #-- compile numerical expression operator for parameters from files
    #-- UTCSR: The University of Texas at Austin Center for Space Research
    #-- EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
    #-- GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
    #-- JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
    #-- JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    regex_pattern = (r'(.*?)-2_(\d+)-(\d+)_(.*?)_({0})_(.*?)_(\d+)(.*?)'
        r'(\.gz|\.gfc)?$').format('UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC')
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- extract parameters from input filename
    if isinstance(input_file, io.IOBase):
        return rx.findall(input_file.filename).pop()
    else:
        return rx.findall(os.path.basename(input_file)).pop()

#-- PURPOSE: read input file and extract contents
def extract_file(input_file, compressed):
    """
    Read input file and extract contents

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    compressed: denotes if the file is compressed
    """
    #-- tilde expansion of input file if not byteIO object
    if not isinstance(input_file, io.IOBase):
        input_file = os.path.expanduser(input_file)
    #-- check if file is uncompressed byteIO object
    if isinstance(input_file, io.IOBase) and not compressed:
        #-- extract spherical harmonic coefficients
        return input_file.read().decode('ISO-8859-1').splitlines()
    else:
        #-- check if file is compressed (read with gzip if gz)
        file_opener = gzip.open if compressed else open
        #-- opening data file to extract spherical harmonic coefficients
        with file_opener(input_file,'rb') as f:
            return f.read().decode('ISO-8859-1').splitlines()
