#!/usr/bin/env python
u"""
read_tellus_geocenter.py
Written by Tyler Sutterley (08/2020)

Reads monthly geocenter spherical harmonic data files from GRACE Tellus
    Technical Notes (TN-13) calculated using GRACE/GRACE-FO measurements and
    Ocean Models of Degree 1

Datasets distributed by NASA PO.DAAC
https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1

Swenson, S., D. Chambers, and J. Wahr, "Estimating geocenter variations
    from a combination of GRACE and ocean model output", J. Geophys. Res.,
    113(B08410), 2008.  doi:10.1029/2007JB005338

Sun, Y., R. Riva, and P. Ditmar, "Observed changes in the Earth's dynamic
    oblateness from GRACE data and geophysical models", J. Geodesy.,
    90(1), 81-89, 2016. doi:10.1007/s00190-015-0852-y

CALLING SEQUENCE:
    geocenter = read_tellus_geocenter(file)

INPUTS:
    geocenter_file: degree 1 file

OUTPUTS:
    C10: cosine d1/o0 spherical harmonic coefficients
    C11: cosine d1/o1 spherical harmonic coefficients
    S11: sine d1/o1 spherical harmonic coefficients
    eC10: cosine d1/o0 spherical harmonic coefficient error
    eC11: cosine d1/o1 spherical harmonic coefficient error
    eS11: sine d1/o1 spherical harmonic coefficient error
    month: GRACE/GRACE-FO month
    time: date of each month in year-decimal

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    JPL: use JPL TN-13 geocenter files with self-attraction and loading

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 08/2019: add catch to verify input geocenter file exists
    Updated 07/2019: month adjustments for new TN-13 geocenter files
        calculate GRACE/GRACE-FO month based on mean time for JPL TN-13 data files
    Updated 06/2019: can use the new JPL TN-13 geocenter files from Tellus
    Updated 10/2018: using future division for python3 Compatibility
    Updated 06/2016: added option HEADER for files that do not have header text
    Updated 04/2015: added time output with convert_calendar_decimal
    Updated 03/2015: minor update to read and regular expression
    Updated 10/2014: rewrote with general code updates.
        using regular expressions to extract data
    Updated 03/2013: changed outputs to be C10, C11, S11 instead of C1, S1
"""
from __future__ import print_function, division

import os
import re
import numpy as np
from gravity_toolkit.convert_calendar_decimal import convert_calendar_decimal

#-- PURPOSE: read geocenter data from PO.DAAC
def read_tellus_geocenter(geocenter_file, HEADER=True, JPL=False):
    """
    Reads monthly geocenter files computed by JPL Tellus using
    GRACE/GRACE-FO measurements and Ocean Models of degree 1

    Arguments
    ---------
    geocenter_file: degree 1 file

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped
    JPL: use JPL TN-13 geocenter files with self-attraction and loading

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

    #-- read degree 1 file and get contents
    with open(os.path.expanduser(geocenter_file),'r') as f:
        file_contents = f.read().splitlines()
    #-- number of lines contained in the file
    file_lines = len(file_contents)

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text
    header_flag = r"end\sof\sheader" if JPL else r"'\(a6,"
    while HEADER:
        #-- file line at count
        line = file_contents[count]
        #-- find header_flag within line to set HEADER flag to False when found
        HEADER = not bool(re.match(header_flag,line))
        #-- add 1 to counter
        count += 1

    #-- number of months within the file
    n_mon = (file_lines - count)//2
    #-- calendar dates
    year = np.zeros((n_mon))
    month = np.zeros((n_mon))
    #-- grace month of data line
    mon = np.zeros((n_mon), dtype=np.int)
    tdec = np.zeros((n_mon))
    #-- spherical harmonic data
    C1 = np.zeros((n_mon,2))
    S1 = np.zeros((n_mon,2))
    eC1 = np.zeros((n_mon,2))
    eS1 = np.zeros((n_mon,2))

    #-- compile numerical expression operator
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- time count
    t = 0
    #-- for every other line:
    for line in file_contents[count:]:
        #-- find numerical instances in line including integers, exponents,
        #-- decimal points and negatives
        line_contents = rx.findall(line)
        #-- calendar year and month
        if JPL:
            #-- start and end dates of month
            start_yr = np.float(line_contents[7][0:4])
            start_mon = np.float(line_contents[7][4:6])
            start_day = np.float(line_contents[7][6:8])
            end_yr = np.float(line_contents[8][0:4])
            end_mon = np.float(line_contents[8][4:6])
            end_day = np.float(line_contents[8][6:8])
            #-- convert date to year decimal
            t_start = convert_calendar_decimal(start_yr,start_mon,DAY=start_day)
            t_end = convert_calendar_decimal(end_yr,end_mon,DAY=end_day)
            #-- calculate mean time
            tdec[t] = np.mean([t_start,t_end])
            year[t] = np.floor(tdec[t])
            month[t] = np.int(12*(tdec[t] % 1) + 1)
        else:
            year[t] = np.float(line_contents[0][0:4])
            month[t] = np.float(line_contents[0][4:6])
            #-- convert date to year decimal
            tdec[t], = convert_calendar_decimal(year[t],month[t])
        #-- grace month
        mon[t] = np.int(12.0*(year[t] - 2002.) + month[t])
        #-- Accelerometer shutoffs complicate the relations between month number
        if (mon[t] == mon[t-1]) and (mon[t-1] in (118,123,160,169,185,205)):
            mon[t] = mon[t-1] + 1
        #-- spherical harmonic order
        m = np.int(line_contents[2])
        #-- extract spherical harmonic data
        C1[t,m] = np.float(line_contents[3])
        S1[t,m] = np.float(line_contents[4])
        eC1[t,m] = np.float(line_contents[5])
        eS1[t,m] = np.float(line_contents[6])
        #-- will only advance in time after reading the
        #-- order 1 coefficients (t+0=t)
        t += m

    #-- reforming outputs to be individual variables
    C10 = np.squeeze(C1[:,0])
    C11 = np.squeeze(C1[:,1])
    S11 = np.squeeze(S1[:,1])
    eC10 = np.squeeze(eC1[:,0])
    eC11 = np.squeeze(eC1[:,1])
    eS11 = np.squeeze(eS1[:,1])

    return {'month':mon, 'C10':C10, 'C11':C11, 'S11':S11, \
        'eC10':eC10, 'eC11':eC11, 'eS11':eS11, 'time':tdec}
