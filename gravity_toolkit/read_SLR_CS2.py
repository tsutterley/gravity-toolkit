#!/usr/bin/env python
u"""
read_SLR_CS2.py
Written by Hugo Lecomte and Tyler Sutterley (04/2022)

Reads monthly degree 2,m (figure axis and azimuthal dependence)
    spherical harmonic data files from satellite laser ranging (SLR)

Dataset distributed by CSR
    http://download.csr.utexas.edu/pub/slr/degree_2/
        C21_S21_RL06.txt or C22_S22_RL06.txt
Dataset distributed by GFZ
    ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
        GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr

CALLING SEQUENCE:
    SLR_2m = read_SLR_CS2(SLR_file)

INPUTS:
    SLR_file:
        CSR 2,1: C21_S21_RL06.txt
        CSR 2,2: C22_S22_RL06.txt
        GFZ: GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat
        GSFC: GSFC_C21_S21.txt

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    DATE: mid-point of monthly solution for calculating 28-day arc averages

OUTPUTS:
    C2m: SLR degree 2 order m cosine stokes coefficients
    S2m: SLR degree 2 order m sine stokes coefficients
    eC2m: SLR degree 2 order m cosine stokes coefficient error
    eS2m: SLR degree 2 order m sine stokes coefficient error
    month: GRACE/GRACE-FO month of measurement (Apr. 2002 = 004)
    time: date of SLR measurement

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_SLR_harmonics.py: low-degree spherical harmonic coefficients from SLR

REFERENCES:
    Cheng et al., " Variations of the Earth's figure axis from satellite
        laser ranging and GRACE", Journal of Geophysical Research,
        116, B01409, (2011). https://doi.org/10.1029/2010JB000850
    Dahle et al., "The GFZ GRACE RL06 Monthly Gravity Field Time Series:
        Processing Details, and Quality Assessment", Remote Sensing,
        11(18), 2116, (2019). https://doi.org/10.3390/rs11182116
    Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential
        GSM Coefficients GFZ RL06 (Level-2B Product)."
        V. 0002. GFZ Data Services, (2019).
        https://doi.org/10.5880/GFZ.GRAVIS_06_L2B
    Chen el al., "Assessment of degree-2 order-1 gravitational changes
        from GRACE and GRACE Follow-on, Earth rotation, satellite laser
        ranging, and models", Journal of Geodesy, 95(38), (2021).
        https://doi.org/10.1007/s00190-021-01492-x

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 11/2021: reader for new weekly 5x5+6,1 fields from NASA GSFC
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 08/2021: output empty spherical harmonic errors for GSFC
    Updated 06/2021: added GSFC 7-day SLR figure axis solutions
    Updated 05/2021: added GFZ GravIS GRACE/SLR low degree solutions
    Updated 04/2021: use adjust_months function to fix special months cases
    Written 11/2020
"""
import os
import re
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.read_SLR_harmonics

#-- PURPOSE: read Degree 2,m data from Satellite Laser Ranging (SLR)
def read_SLR_CS2(SLR_file, ORDER=1, DATE=None, HEADER=True):
    """
    Reads CS2,m spherical harmonic coefficients from SLR measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    ORDER: int, default 1
        Spherical harmonic order to extract from low-degree fields
    DATE: float or NoneType, default None
        Mid-point of monthly solution for calculating 28-day arc averages
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    C2m: float
        SLR degree 2 order m cosine stokes coefficients
    S2m: float
        SLR degree 2 order m sine stokes coefficients
    eC2m: float
        SLR degree 2 order m cosine stokes coefficient error
    eS2m: float
        SLR degree 2 order m sine stokes coefficient error
    month: int
        GRACE/GRACE-FO month of measurement
    time: float
        date of SLR measurement
    """

    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')
    #-- output dictionary with input data
    dinput = {}

    if bool(re.search(r'GSFC_C2(\d)_S2(\d)',SLR_file,re.I)):
        #-- 7-day arc SLR file produced by GSFC
        #-- input variable names and types
        dtype = {}
        dtype['names'] = ('time','C2','S2')
        dtype['formats'] = ('f','f8','f8')
        #-- read SLR 2,1 file from GSFC
        #-- Column 1: Approximate mid-point of 7-day solution (years)
        #-- Column 2: Solution from SLR (normalized)
        #-- Column 3: Solution from SLR (normalized)
        content = np.loadtxt(os.path.expanduser(SLR_file),dtype=dtype)
        #-- duplicate time and harmonics
        tdec = np.repeat(content['time'],7)
        c2m = np.repeat(content['C2'],7)
        s2m = np.repeat(content['S2'],7)
        #-- calculate daily dates to use in centered moving average
        tdec += (np.mod(np.arange(len(tdec)),7) - 3.5)/365.25
        #-- number of dates to use in average
        n_neighbors = 28
        #-- calculate 28-day moving-average solution from 7-day arcs
        dinput['time'] = np.zeros_like(DATE)
        dinput['C2m'] = np.zeros_like(DATE,dtype='f8')
        dinput['S2m'] = np.zeros_like(DATE,dtype='f8')
        #-- no estimated spherical harmonic errors
        dinput['eC2m'] = np.zeros_like(DATE,dtype='f8')
        dinput['eS2m'] = np.zeros_like(DATE,dtype='f8')
        for i,D in enumerate(DATE):
            isort = np.argsort((tdec - D)**2)[:n_neighbors]
            dinput['time'][i] = np.mean(tdec[isort])
            dinput['C2m'][i] = np.mean(c2m[isort])
            dinput['S2m'][i] = np.mean(s2m[isort])
        #-- GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    elif bool(re.search(r'gsfc_slr_5x5c61s61',SLR_file,re.I)):
        #-- read 5x5 + 6,1 file from GSFC and extract coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        #-- duplicate time and harmonics
        tdec = np.repeat(Ylms['time'],7)
        c2m = np.repeat(Ylms['clm'][2,ORDER],7)
        s2m = np.repeat(Ylms['slm'][2,ORDER],7)
        #-- calculate daily dates to use in centered moving average
        tdec += (np.mod(np.arange(len(tdec)),7) - 3.5)/365.25
        #-- number of dates to use in average
        n_neighbors = 28
        #-- calculate 28-day moving-average solution from 7-day arcs
        dinput['time'] = np.zeros_like(DATE)
        dinput['C2m'] = np.zeros_like(DATE,dtype='f8')
        dinput['S2m'] = np.zeros_like(DATE,dtype='f8')
        #-- no estimated spherical harmonic errors
        dinput['eC2m'] = np.zeros_like(DATE,dtype='f8')
        dinput['eS2m'] = np.zeros_like(DATE,dtype='f8')
        for i,D in enumerate(DATE):
            isort = np.argsort((tdec - D)**2)[:n_neighbors]
            dinput['time'][i] = np.mean(tdec[isort])
            dinput['C2m'][i] = np.mean(c2m[isort])
            dinput['S2m'][i] = np.mean(s2m[isort])
        #-- GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    elif bool(re.search(r'C2(\d)_S2(\d)_(RL\d{2})',SLR_file,re.I)):
        #-- SLR RL06 file produced by CSR
        #-- input variable names and types
        dtype = {}
        dtype['names'] = ('time','C2','S2','eC2','eS2',
            'C2aod','S2aod','start','end')
        dtype['formats'] = ('f','f8','f8','f','f','f','f','f','f')
        #-- read SLR 2,1 or 2,2 RL06 file from CSR
        #-- header text is commented and won't be read
        #-- Column 1: Approximate mid-point of monthly solution (years)
        #-- Column 2: Solution from SLR (normalized)
        #-- Column 3: Solution from SLR (normalized)
        #-- Column 4: Solution sigma (1E-10)
        #-- Column 5: Solution sigma (1E-10)
        #-- Column 6: Mean value of Atmosphere-Ocean De-aliasing model (1E-10)
        #-- Column 7: Mean value of Atmosphere-Ocean De-aliasing model (1E-10)
        #-- Columns 8-9: Start and end dates of data used in solution
        content = np.loadtxt(os.path.expanduser(SLR_file),dtype=dtype)
        #-- date and GRACE/GRACE-FO month
        dinput['time'] = content['time'].copy()
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
        #-- remove the monthly mean of the AOD model
        dinput['C2m'] = content['C2'] - content['C2aod']*10**-10
        dinput['S2m'] = content['S2'] - content['S2aod']*10**-10
        #-- scale SLR solution sigmas
        dinput['eC2m'] = content['eC2']*10**-10
        dinput['eS2m'] = content['eS2']*10**-10
    elif bool(re.search(r'GRAVIS-2B_GFZOP',SLR_file,re.I)):
        #-- Combined GRACE/SLR solution file produced by GFZ
        #-- Column  1: MJD of BEGINNING of solution data span
        #-- Column  2: Year and fraction of year of BEGINNING of solution span
        #-- Column  9: Replacement C(2,1)
        #-- Column 10: Replacement C(2,1) - mean C(2,1) (1.0E-10)
        #-- Column 11: C(2,1) formal standard deviation (1.0E-12)
        #-- Column 12: Replacement S(2,1)
        #-- Column 13: Replacement S(2,1) - mean S(2,1) (1.0E-10)
        #-- Column 14: S(2,1) formal standard deviation (1.0E-12)
        with open(os.path.expanduser(SLR_file), mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        #-- number of lines contained in the file
        file_lines = len(file_contents)

        #-- counts the number of lines in the header
        count = 0
        #-- Reading over header text
        while HEADER:
            #-- file line at count
            line = file_contents[count]
            #-- find PRODUCT: within line to set HEADER flag to False when found
            HEADER = not bool(re.match(r'PRODUCT:+',line))
            #-- add 1 to counter
            count += 1

        #-- number of months within the file
        n_mon = file_lines - count
        #-- date and GRACE/GRACE-FO month
        dinput['time'] = np.zeros((n_mon))
        dinput['month'] = np.zeros((n_mon),dtype=int)
        #-- monthly spherical harmonic replacement solutions
        dinput['C2m'] = np.zeros((n_mon))
        dinput['S2m'] = np.zeros((n_mon))
        #-- monthly spherical harmonic formal standard deviations
        dinput['eC2m'] = np.zeros((n_mon))
        dinput['eS2m'] = np.zeros((n_mon))
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall(r'[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            count = len(line_contents)
            #-- check for empty lines
            if (count > 0):
                #-- reading decimal year for start of span
                dinput['time'][t] = np.float64(line_contents[1])
                #-- Spherical Harmonic data for line
                dinput['C2m'][t] = np.float64(line_contents[8])
                dinput['eC2m'][t] = np.float64(line_contents[10])*1e-10
                dinput['S2m'][t] = np.float64(line_contents[11])
                dinput['eS2m'][t] = np.float64(line_contents[13])*1e-10
                #-- GRACE/GRACE-FO month of SLR solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                #-- add to t count
                t += 1
        #-- truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])

    #-- return the SLR-derived degree 2 solutions
    return dinput
