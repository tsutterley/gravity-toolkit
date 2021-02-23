#!/usr/bin/env python
u"""
read_SLR_geocenter.py
Written by Tyler Sutterley (02/2021)

Reads monthly geocenter files from satellite laser ranging provided by CSR
    http://download.csr.utexas.edu/pub/slr/geocenter/
    RL04: GCN_RL04.txt
    RL05: GCN_RL05.txt
New CF-CM geocenter dataset to reflect the true degree-1 mass variations
    http://download.csr.utexas.edu/pub/slr/geocenter/geocenter/README_L1_L2
    http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
New geocenter solutions from Minkang Cheng
    http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s

CALLING SEQUENCE:
    geocenter = read_SLR_geocenter(geocenter_file)

INPUTS:
    geocenter_file: degree 1 file

OPTIONS:
    RADIUS: Earth's radius for calculating spherical harmonics from SLR data
    HEADER: rows of data to skip when importing data
    COLUMNS: column names of ascii file
        time: date in decimal-years
        X: X-component of geocenter variation
        Y: Y-component of geocenter variation
        Z: Z-component of geocenter variation
        X_sigma: X-component uncertainty
        Y_sigma: Y-component uncertainty
        Z_sigma: Z-component uncertainty

OUTPUTS:
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    eC10: cosine d1/o0 spherical harmonic coefficient error
    eC11: cosine d1/o1 spherical harmonic coefficient error
    eS11: sine d1/o1 spherical harmonic coefficient error
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    geocenter.py: converts between spherical harmonics and geocenter variations
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: added option COLUMNS to generalize the ascii data format
        replaced numpy loadtxt with generic read using regular expressions
        using utilities from time module for operations
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 08/2019: add catch to verify input geocenter file exists
    Updated 06/2019: added option RADIUS for setting the Earth's radius
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 04/2017: parallels updates to geocenter function INVERSE option
        use enumerate to iterate over dates.  added option HEADER for headers
    Updated 06/2016: using __future__ print function
    Updated 05/2016: use geocenter files from 6-hour AOD1b glo Ylms calculated
        in aod1b_geocenter.py
    Updated 09/2015: added second function for AOD corrected geocenter values
    Updated 04/2015: using Julian dates to calculate GRACE/GRACE-FO month
    Written 08/2013
"""
from __future__ import print_function

import os
import re
import time
import numpy as np
import gravity_toolkit.time
from gravity_toolkit.geocenter import geocenter

#-- PURPOSE: read geocenter data from Satellite Laser Ranging (SLR)
def read_SLR_geocenter(geocenter_file, RADIUS=None, HEADER=0,
    COLUMNS=['time','X','Y','Z','X_sigma','Y_sigma','Z_sigma']):
    """
    Reads monthly geocenter files from satellite laser ranging

    Arguments
    ---------
    geocenter_file: Satellite Laser Ranging file

    Keyword arguments
    -----------------
    RADIUS: Earth's radius for calculating spherical harmonics from SLR data
    HEADER: rows of data to skip when importing data
    COLUMNS: column names of ascii file
        time: date in decimal-years
        X: X-component of geocenter variation
        Y: Y-component of geocenter variation
        Z: Z-component of geocenter variation
        X_sigma: X-component uncertainty
        Y_sigma: Y-component uncertainty
        Z_sigma: Z-component uncertainty

    Returns
    -------
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    eC10: cosine d1/o0 spherical harmonic coefficient error
    eC11: cosine d1/o1 spherical harmonic coefficient error
    eS11: sine d1/o1 spherical harmonic coefficient error
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal
    """

    #-- check that geocenter file exists
    if not os.access(os.path.expanduser(geocenter_file), os.F_OK):
        raise IOError('Geocenter file not found in file system')

    #-- Input geocenter file and split lines
    with open(os.path.expanduser(geocenter_file),'r') as f:
        file_contents = f.read().splitlines()
    ndate = len(file_contents)-HEADER

    #-- compile regular expression operator to find numerical instances
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- initializing output data
    #-- Degree 1 Stokes Coefficients
    C10 = np.zeros((ndate))
    C11 = np.zeros((ndate))
    S11 = np.zeros((ndate))
    #-- Degree 1 Stokes Coefficient Errors
    eC10 = np.zeros((ndate))
    eC11 = np.zeros((ndate))
    eS11 = np.zeros((ndate))
    #-- Date information
    date = np.zeros((ndate))
    JD = np.zeros((ndate))
    mon = np.zeros((ndate), dtype=np.int32)

    #-- for each date
    for t,file_line in enumerate(file_contents[HEADER:]):
        #-- find numerical instances in line
        #-- replacing fortran double precision exponential
        line_contents = rx.findall(file_line.replace('D','E'))
        #-- extract date
        date[t] = np.float(line_contents[COLUMNS.index('time')])
        #-- extract geocenter variations
        X = np.float(line_contents[COLUMNS.index('X')])
        Y = np.float(line_contents[COLUMNS.index('Y')])
        Z = np.float(line_contents[COLUMNS.index('Z')])
        X_sigma = np.float(line_contents[COLUMNS.index('X_sigma')])
        Y_sigma = np.float(line_contents[COLUMNS.index('Y_sigma')])
        Z_sigma = np.float(line_contents[COLUMNS.index('Z_sigma')])
        #-- converting from geocenter into spherical harmonics
        CS1 = geocenter(X=X,Y=Y,Z=Z,RADIUS=RADIUS,INVERSE=True)
        dCS1 = geocenter(X=X_sigma,Y=Y_sigma,Z=Z_sigma,
            RADIUS=RADIUS, INVERSE=True)
        #-- output harmonics
        C10[t],C11[t],S11[t] = (CS1['C10'], CS1['C11'], CS1['S11'])
        eC10[t],eC11[t],eS11[t] = (dCS1['C10'], dCS1['C11'], dCS1['S11'])

        #-- calendar year of date
        year = np.floor(date[t])
        #-- check if year is a leap year
        dpy = 366.0 if ((year % 4) == 0) else 365.0
        #-- calculation of day of the year (with decimals for fraction of day)
        DofY = dpy*(date[t] % 1)
        #-- Calculation of the Julian date from year and DofY
        JD[t] = np.float(367.0*year -
            np.floor(7.0*(year + np.floor(10.0/12.0))/4.0) -
            np.floor(3.0*(np.floor((year - 8.0/7.0)/100.0) + 1.0)/4.0) +
            np.floor(275.0/9.0) + DofY + 1721028.5)
        #-- convert the julian date into calendar dates (hour, day, month, year)
        cal_date = gravity_toolkit.time.convert_julian(JD[t])
        #-- calculate the GRACE/GRACE-FO month (Apr02 == 004)
        #-- https://grace.jpl.nasa.gov/data/grace-months/
        mon[t] = 12*(cal_date['year']-2002) + cal_date['month']

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    mon = gravity_toolkit.time.adjust_months(mon)

    return {'C10':C10, 'C11':C11, 'S11':S11, 'eC10':eC10, 'eC11':eC11,
        'eS11':eS11, 'month':mon, 'time':date}

#-- special function for outputting AOD corrected SLR geocenter values
#-- need to run aod1b_geocenter.py to calculate the monthly geocenter dealiasing
def aod_corrected_SLR_geocenter(geocenter_file, DREL, RADIUS=None, HEADER=0,
    COLUMNS=[]):
    """
    Reads monthly geocenter files from satellite laser ranging corrected
    for non-tidal ocean and atmospheric variation

    Arguments
    ---------
    geocenter_file: Satellite Laser Ranging file
    DREL: GRACE/GRACE-FO data release

    Keyword arguments
    -----------------
    RADIUS: Earth's radius for calculating spherical harmonics from SLR data
    HEADER: rows of data to skip when importing data
    COLUMNS: column names of ascii file
        time: date in decimal-years
        X: X-component of geocenter variation
        Y: Y-component of geocenter variation
        Z: Z-component of geocenter variation
        X_sigma: X-component uncertainty
        Y_sigma: Y-component uncertainty
        Z_sigma: Z-component uncertainty

    Returns
    -------
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    eC10: cosine d1/o0 spherical harmonic coefficient error
    eC11: cosine d1/o1 spherical harmonic coefficient error
    eS11: sine d1/o1 spherical harmonic coefficient error
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal
    """
    #-- directory setup for AOD1b data starting with input degree 1 file
    #-- this will verify that the input paths work
    AOD1B_dir = os.path.abspath(os.path.join(geocenter_file,os.path.pardir,
        os.path.pardir,'AOD1B',DREL,'geocenter'))

    #-- Input geocenter file and split lines
    with open(os.path.expanduser(geocenter_file),'r') as f:
        file_contents = f.read().splitlines()
    ndate = len(file_contents)-HEADER

    #-- compile regular expression operator to find numerical instances
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- initializing output data
    #-- Degree 1 Stokes Coefficients
    C10 = np.zeros((ndate))
    C11 = np.zeros((ndate))
    S11 = np.zeros((ndate))
    #-- Degree 1 Stokes Coefficient Errors
    eC10 = np.zeros((ndate))
    eC11 = np.zeros((ndate))
    eS11 = np.zeros((ndate))
    #-- Date information
    date = np.zeros((ndate))
    JD = np.zeros((ndate))
    mon = np.zeros((ndate), dtype=np.int32)

    #-- for each date
    for t,file_line in enumerate(file_contents[HEADER:]):
        #-- find numerical instances in line
        #-- replacing fortran double precision exponential
        line_contents = rx.findall(file_line.replace('D','E'))
        #-- extract date
        date[t] = np.float(line_contents[COLUMNS.index('time')])
        #-- extract geocenter variations
        X = np.float(line_contents[COLUMNS.index('X')])
        Y = np.float(line_contents[COLUMNS.index('Y')])
        Z = np.float(line_contents[COLUMNS.index('Z')])
        X_sigma = np.float(line_contents[COLUMNS.index('X_sigma')])
        Y_sigma = np.float(line_contents[COLUMNS.index('Y_sigma')])
        Z_sigma = np.float(line_contents[COLUMNS.index('Z_sigma')])
        #-- converting from geocenter into spherical harmonics
        CS1 = geocenter(X=X,Y=Y,Z=Z,RADIUS=RADIUS,INVERSE=True)
        dCS1 = geocenter(X=X_sigma,Y=Y_sigma,Z=Z_sigma,
            RADIUS=RADIUS, INVERSE=True)

        #-- calendar year of date
        year = np.floor(date[t])
        #-- check if year is a leap year
        dpy = 366.0 if ((year % 4) == 0) else 365.0
        #-- calculation of day of the year (with decimals for fraction of day)
        DofY = dpy*(date[t] % 1)
        #-- Calculation of the Julian date from year and DofY
        JD[t] =np.float(367.*year - np.floor(7.*(year + np.floor(10./12.))/4.) -
            np.floor(3.0*(np.floor((year - 8.0/7.0)/100.0) + 1.0)/4.0) +
            np.floor(275.0/9.0) + DofY + 1721028.5)
        #-- convert the julian date into calendar dates (hour, day, month, year)
        cal_date = gravity_toolkit.time.convert_julian(JD[t], ASTYPE=np.int)
        #-- full path to AOD geocenter for month (using glo coefficients)
        args = (DREL, 'glo', cal_date['year'], cal_date['month'])
        AOD1B_file = 'AOD1B_{0}_{1}_{2:4d}_{3:02d}.txt'.format(*args)
        Ylms = read_AOD1b_geocenter(os.path.join(AOD1B_dir,AOD1B_file),
            cal_date['month'])
        #-- remove AOD from output harmonics
        C10[t] = CS1['C10'] - Ylms['C10']
        C11[t] = CS1['C11'] - Ylms['C11']
        S11[t] = CS1['S11'] - Ylms['S11']
        eC10[t],eC11[t],eS11[t] = (dCS1['C10'], dCS1['C11'], dCS1['S11'])
        #-- calculate the GRACE/GRACE-FO month (Apr02 == 004)
        #-- https://grace.jpl.nasa.gov/data/grace-months/
        mon[t] = 12*(cal_date['year']-2002) + cal_date['month']

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    mon = gravity_toolkit.time.adjust_months(mon)

    return {'C10':C10, 'C11':C11, 'S11':S11, 'eC10':eC10, 'eC11':eC11,
        'eS11':eS11, 'month':mon, 'time':date}

#-- PURPOSE: read AOD1b geocenter for month and calculate the mean harmonics
#-- need to run aod1b_geocenter.py to write these monthly geocenter files
def read_AOD1b_geocenter(AOD1B_file, calendar_month):
    """
    Reads monthly non-tidal ocean and atmospheric variation geocenter files

    Arguments
    ---------
    AOD1B_file: de-aliasing data file
    calendar_month: calendar month of data
    """
    #-- check that file exists
    if not os.access(AOD1B_file, os.F_OK):
        raise IOError('AOD1b File {0} not in File System'.format(AOD1B_file))
    #-- read AOD1b geocenter skipping over commented header text
    with open(AOD1B_file, 'r') as f:
        file_contents=[i for i in f.read().splitlines() if not re.match(r'#',i)]
    #-- extract X,Y,Z from each line in the file
    #-- first column: ISO-formatted date and time
    #-- second-fourth columns: X, Y and Z geocenter variations
    n_lines = len(file_contents)
    X = np.zeros((n_lines))
    Y = np.zeros((n_lines))
    Z = np.zeros((n_lines))
    month = np.zeros((n_lines),dtype=np.int)
    for i,line in enumerate(file_contents):
        line_contents = line.split()
        AOD1B_time = time.strptime(line_contents[0],'%Y-%m-%dT%H:%M:%S')
        month[i] = AOD1B_time.tm_mon
        X[i],Y[i],Z[i] = np.array(line_contents[1:],dtype=np.float)
    #-- use only dates within month (should be all)
    ii, = np.nonzero(month == calendar_month)
    #-- convert mean X,Y,Z into spherical harmonics
    return geocenter(X=X[ii].mean(),Y=Y[ii].mean(),Z=Z[ii].mean(),INVERSE=True)
