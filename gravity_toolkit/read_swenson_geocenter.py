#!/usr/bin/env python
u"""
read_swenson_geocenter.py
Written by Tyler Sutterley (05/2021)

Reads monthly geocenter coefficients from GRACE measurements and
    Ocean Models of Degree 1 provided by Sean Swenson in mm w.e.

Swenson, S., D. Chambers, and J. Wahr, "Estimating geocenter variations
    from a combination of GRACE and ocean model output", J. Geophys. Res.,
    113(B08410), 2008.  doi:10.1029/2007JB005338

CALLING SEQUENCE:
    geocenter = read_swenson_geocenter(geocenter_file)

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

UPDATE HISTORY:
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use file not found exceptions
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: using utilities from time module
    Updated 08/2019: add catch to verify input geocenter file exists
    Updated 10/2018: verify integers for python3 compatibility
    Updated 03/2017: added catch for HEADER flag
    Updated 06/2016: added option HEADER for files that do not have header text
    Updated 10/2015: checks if months are included as the last column
        and will import.  else will calculate GRACE month from the date
        Updated regex operator to include pure integers as well
    Updated 04/2015: using Julian date to determine GRACE month
    Updated 03/2015: minor update to read and regular expression
    Updated 10/2014: rewrote with general code updates
        using regular expressions to extract data
        using better algorithm to find grace month
    Written 11/2012
"""
import os
import re
import numpy as np
import gravity_toolkit.time

#-- PURPOSE: read geocenter data from Sean Swenson
def read_swenson_geocenter(geocenter_file, HEADER=True):
    """
    Reads monthly geocenter files computed by Sean Swenson using
    GRACE/GRACE-FO measurements and Ocean Models of degree 1

    Arguments
    ---------
    geocenter_file: degree 1 file

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped

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

    #-- read degree 1 file and get contents
    with open(os.path.expanduser(geocenter_file),'r') as f:
        file_contents = f.read().splitlines()
    #-- number of lines contained in the file
    file_lines = len(file_contents)

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text
    while HEADER and (count < file_lines):
        #-- file line at count
        line = file_contents[count]
        #-- find Time within line to set HEADER flag to False when found
        HEADER = not bool(re.search(r"Time",line))
        #-- add 1 to counter
        count += 1

    #-- catch to see if HEADER flag was not set to false
    if HEADER:
        raise IOError('Data lines not found in file {0}'.format(geocenter_file))

    #-- number of months within the file
    n_mon = np.int64(file_lines - count)
    C10 = np.zeros((n_mon))
    C11 = np.zeros((n_mon))
    S11 = np.zeros((n_mon))
    date = np.zeros((n_mon))
    JD = np.zeros((n_mon))
    mon = np.zeros((n_mon), dtype=np.int64)

    #-- Average Density of the Earth [g/cm^3]
    rho_e = 5.517
    #-- Average Radius of the Earth [cm]
    rad_e = 6.371e8

    #-- compile numerical expression operator
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- for every other line:
    for t, line in enumerate(file_contents[count:]):
        #-- find numerical instances in line including integers, exponents,
        #-- decimal points and negatives
        line_contents = rx.findall(line)

        #-- extacting data
        date[t]=np.float64(line_contents[0])
        #-- Convert from mm w.e. into cm w.e. (0.1)
        #-- then convert from cm w.e. into normalized geoid coeffs (rho_e*rad_e)
        C10[t]=0.1*np.float64(line_contents[1])/(rho_e*rad_e)
        C11[t]=0.1*np.float64(line_contents[2])/(rho_e*rad_e)
        S11[t]=0.1*np.float64(line_contents[3])/(rho_e*rad_e)

        #-- calculate the GRACE months
        #-- calendar year of date
        year = np.floor(date[t])

        #-- check if year is a leap year
        if ((year % 4) == 0):
            #-- Leap Year
            dpy = 366
        else:
            #-- Standard Year
            dpy = 365
        #-- calculation of day of the year
        DofY = dpy*(date[t] % 1)
        #-- Calculation of the Julian date from year and DofY
        JD[t] = np.float64(367.0*year - \
            np.floor(7.0*(year + np.floor(10.0/12.0))/4.0) - \
            np.floor(3.0*(np.floor((year - 8.0/7.0)/100.0) + 1.0)/4.0) + \
            np.floor(275.0/9.0) + DofY + 1721028.5)

        #-- months are included as last column
        if (len(line_contents) == 5):
            mon[t] = np.int64(line_contents[4])
        else:
            #-- months to be calculated from date
            #-- convert the julian date into calendar dates (day, month, year)
            cal_date = gravity_toolkit.time.convert_julian(JD[t])
            #-- calculate the GRACE month (Apr02 == 004)
            #-- https://grace.jpl.nasa.gov/data/grace-months/
            #-- Notes on special months (e.g. 119, 120) below
            mon[t] = gravity_toolkit.time.calendar_to_grace(cal_date['year'],
                cal_date['month'])

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    mon = gravity_toolkit.time.adjust_months(mon)

    return {'C10':C10, 'C11':C11, 'S11':S11, 'month':mon, 'date':date}
