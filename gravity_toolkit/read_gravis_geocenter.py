#!/usr/bin/env python
u"""
read_gravis_geocenter.py
Written by Tyler Sutterley (05/2021)

Reads monthly geocenter spherical harmonic data files from
    GFZ GravIS calculated using GRACE/GRACE-FO measurements
    and Ocean Models of degree 1

Dataset distributed by GFZ
    ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
        GRAVIS-2B_GFZOP_GEOCENTER_0002.dat

CALLING SEQUENCE:
    geocenter = read_gravis_geocenter(geocenter_file)

INPUTS:
    geocenter_file: degree 1 file

OUTPUTS:
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

REFERENCES:
    Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential
        GSM Coefficients GFZ RL06 (Level-2B Product)."
        V. 0002. GFZ Data Services, (2019).
        http://doi.org/10.5880/GFZ.GRAVIS_06_L2B

UPDATE HISTORY:
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Written 05/2021
"""
import os
import re
import numpy as np
import gravity_toolkit.time

#-- PURPOSE: read geocenter data from GFZ GravIS SLR/GRACE solutions
def read_gravis_geocenter(geocenter_file, HEADER=True):
    """
    Reads monthly geocenter spherical harmonic data files from
        GFZ GravIS calculated using GRACE/GRACE-FO measurements
        and Ocean Models of degree 1

    Arguments
    ---------
    geocenter_file: degree 1 file

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped (default: True)

    Returns
    -------
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal
    """

    #-- check that geocenter file exists
    if not os.access(os.path.expanduser(geocenter_file), os.F_OK):
        raise FileNotFoundError('Geocenter file not found in file system')
    #-- output dictionary with input data
    dinput = {}

    #-- Combined GRACE/SLR geocenter solution file produced by GFZ GravIS
    #-- Column  1: MJD of BEGINNING of solution data span
    #-- Column  2: Year and fraction of year of BEGINNING of solution data span
    #-- Column  3: Coefficient C(1,0)
    #-- Column  4: Coefficient C(1,0) - mean C(1,0) (1.0E-10)
    #-- Column  5: C(1,0) uncertainty (1.0E-10)
    #-- Column  6: Coefficient C(1,1)
    #-- Column  7: Coefficient C(1,1) - mean C(1,1) (1.0E-10)
    #-- Column  8: C(1,1) uncertainty (1.0E-10)
    #-- Column  9: Coefficient S(1,1)
    #-- Column 10: Coefficient S(1,1) - mean S(1,1) (1.0E-10)
    #-- Column 11: S(1,1) uncertainty (1.0E-10)

    with open(os.path.expanduser(geocenter_file),'r') as f:
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
    dinput['C10'] = np.zeros((n_mon))
    dinput['C11'] = np.zeros((n_mon))
    dinput['S11'] = np.zeros((n_mon))
    #-- monthly spherical harmonic formal standard deviations
    dinput['eC10'] = np.zeros((n_mon))
    dinput['eC11'] = np.zeros((n_mon))
    dinput['eS11'] = np.zeros((n_mon))
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
            dinput['C10'][t] = np.float64(line_contents[2])
            dinput['C11'][t] = np.float64(line_contents[5])
            dinput['S11'][t] = np.float64(line_contents[8])
            #-- monthly spherical harmonic formal standard deviations
            dinput['eC10'][t] = np.float64(line_contents[4])*1e-10
            dinput['eC11'][t] = np.float64(line_contents[7])*1e-10
            dinput['eS11'][t] = np.float64(line_contents[10])*1e-10
            #-- GRACE/GRACE-FO month of geocenter solutions
            dinput['month'][t] = 1+np.round((dinput['time'][t]-2002.)*12.)
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

    #-- return the GFZ GravIS geocenter solutions
    return dinput
