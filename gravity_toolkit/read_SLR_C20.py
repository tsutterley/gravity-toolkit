#!/usr/bin/env python
u"""
read_SLR_C20.py
Written by Tyler Sutterley (07/2020)

Reads in C20 spherical harmonic coefficients derived from SLR measurements

Dataset distributed by NASA PO.DAAC
    https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/grace/docs
        TN-05_C20_SLR.txt
        TN-07_C20_SLR.txt
        TN-11_C20_SLR.txt
        TN-14_C30_C30_GSFC_SLR.txt
Additional dataset distributed by UTCSR
    ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL05.txt

REFERENCE:
    Cheng, M. and Tapley, B. D., "Variations in the Earth's oblateness during
        the past 28 years", Journal of Geophysical Research: Solid Earth,
        109(B9), B09402, 2004. 10.1029/2004JB003028

CALLING SEQUENCE:
    SLR_C20 = read_SLR_C20(SLR_file)

INPUTS:
    SLR_file:
        RL04: TN-05_C20_SLR.txt
        RL05: TN-07_C20_SLR.txt
        RL06: TN-11_C20_SLR.txt
        CSR: C20_RL05.txt

OUTPUTS:
    data: SLR degree 2 order 0 cosine stokes coefficients (C20)
    error: SLR degree 2 order 0 cosine stokes coefficient error (eC20)
    month: GRACE/GRACE-FO month of measurement (Apr. 2002 = 004)
    date: date of SLR measurement

OPTIONS:
    AOD: remove background De-aliasing product from the SLR solution (for CSR)
    HEADER: file contains header text to be skipped (default: True)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_calendar_decimal.py: Return the decimal year for a calendar date

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 08/2019: add catch to verify input SLR file exists
    Updated 07/2019: added tilde-expansion of input SLR file
    Updated 06/2019: added new GRACE-FO special month (October 2018)
    Updated 11/2018: new TN-11 files only list GRACE months available
    Updated 06/2016: added option HEADER for files that do not have header text
    Updated 05/2016: added option AOD to not remove the AOD correction
    Updated 03/2016: minor update to read PO.DAAC
    Updated 05/2015: minor change to file determination (only regular expressions)
    Updated 02/2015: updated UT/CSR portion and comments
    Updated 09/2014: rewrite of the TN-07 read program
        using regular expressions and convert_calendar_decimal
    Updated 01/2014: updated to use UT/CSR monthly time-series
        as an alternative to PO.DAAC as it is updated more regularly
    Updated 05/2013: adapted for python
    Updated 09/2012: Changed month scheme to output.
        Used to remove the GRACE missing months in this program by feeding in the GRACE months
        BUT, as the new SLR files start with an earlier date, decided to parallel
        the degree-1 read program, and remove the missing months in the read_grace program
    Updated 06/2012: OVERHAUL of dating and modification for 'special' GRACE months
        Initiated from an incorrect date tag in the SLR data file
        New dating will convert from the MJD file into date fraction
        Some GRACE 'months' have the accelerometer turned off
            for half the month to preserve battery power
        These months use half of the prior month in the GRACE global gravity solution
        For these months the SLR file has a second dataline for the modified period
        Will use these marked (*) data to replace the GRACE C2,0
        ALSO converted the mon and slrdate inputs into options
    Updated 01/2012: Updated to feed in SLR file from outside
        Update makes this program universal for each computer
        Won't have to update file on each computer pointing to the SLR file
        Will accommodate upcoming GRACE RL05, which will use different SLR files
    Written 12/2011
"""
import os
import re
import numpy as np
from gravity_toolkit.convert_julian import convert_julian
from gravity_toolkit.convert_calendar_decimal import convert_calendar_decimal

#-- PURPOSE: read oblateness data from Satellite Laser Ranging (SLR)
def read_SLR_C20(SLR_file, HEADER=True, AOD=True):
    """
    Reads C20 spherical harmonic coefficients from SLR measurements

    Arguments
    ---------
    SLR_file: Satellite Laser Ranging file

    Keyword arguments
    -----------------
    AOD: remove background De-aliasing product from the SLR solution
    HEADER: file contains header text to be skipped

    Returns
    -------
    data: SLR degree 2 order 0 cosine stokes coefficients
    error: SLR degree 2 order 0 cosine stokes coefficient error
    month: GRACE/GRACE-FO month of measurement
    date: date of SLR measurement
    """

    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise IOError('SLR file not found in file system')

    #-- determine if imported file is from PO.DAAC or CSR
    if bool(re.search('C20_RL\d+',SLR_file)):
        #-- SLR C20 file from CSR
        #-- Just for checking new months when TN series isn't up to date as the
        #-- SLR estimates always use the full set of days in each calendar month.
        #-- format of the input file (note 64 bit floating point for C20)
        #-- Column 1: Approximate mid-point of monthly solution (years)
        #-- Column 2: C20 from SLR (normalized)
        #-- Column 3: Delta C20 relative to a mean value of -4.841694723127E-4 (1E-10)
        #-- Column 4: Solution sigma (1E-10)
        #-- Column 5: Mean value of Atmosphere-Ocean De-aliasing model (1E-10)
        #-- Columns 6-7: Start and end dates of data used in solution
        dtype = {}
        dtype['names'] = ('time','C20','delta','sigma','AOD','start','end')
        dtype['formats'] = ('f','f8','f','f','f','f','f')
        #-- header text is commented and won't be read
        file_input = np.loadtxt(os.path.expanduser(SLR_file),dtype=dtype)
        #-- date and GRACE/GRACE-FO month
        tdec = file_input['time']
        grace_month = 1 + np.floor((tdec-2002.)*12.)
        C20 = file_input['C20']
        eC20 = file_input['sigma']*1e-10
        #-- Background gravity model includes solid earth and ocean tides, solid
        #-- earth and ocean pole tides, and the Atmosphere-Ocean De-aliasing
        #-- product. The monthly mean of the AOD model has been restored.
        if AOD:#-- Removing AOD product that was restored in the solution
            C20 -= file_input['AOD']*1e-10
    elif bool(re.search('TN-(11|14)',SLR_file)):
        #-- SLR C20 RL06 file from PO.DAAC
        with open(os.path.expanduser(SLR_file),'r') as f:
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
            HEADER = not bool(re.match('PRODUCT:+',line,re.IGNORECASE))
            #-- add 1 to counter
            count += 1

        #-- number of months within the file
        n_mon = file_lines - count
        date_conv = np.zeros((n_mon))
        C20_input = np.zeros((n_mon))
        eC20_input = np.zeros((n_mon))
        mon = np.zeros((n_mon),dtype=np.int)
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall('[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            #-- check for empty lines as there are
            #-- slight differences in RL04 TN-05_C20_SLR.txt
            #-- with blanks between the PRODUCT: line and the data
            count = len(line_contents)
            #-- if count is greater than 0
            if (count > 0):
                #-- modified julian date for line
                MJD = np.float(line_contents[0])
                #-- converting from MJD into month, day and year
                YY,MM,DD,hh,mm,ss = convert_julian(MJD+2400000.5,FORMAT='tuple')
                #-- converting from month, day, year into decimal year
                date_conv[t] = convert_calendar_decimal(YY, MM, DAY=DD, HOUR=hh)
                #-- Spherical Harmonic data for line
                C20_input[t] = np.float(line_contents[2])
                eC20_input[t] = np.float(line_contents[4])*1e-10
                #-- GRACE/GRACE-FO month of SLR solutions
                mon[t] = 1 + np.round((date_conv[t]-2002.)*12.)
                #-- The GRACE/GRACE-FO 'Special Months'
                #-- (November 2011, December 2012, April 2012, October 2019)
                #-- Accelerometer shutoffs make the relation between month number
                #-- and date more complicated as days from other months are used
                #-- Nov11 (month 119) is centered in Oct11 (118)
                #-- May15 (month 161) is centered in Apr15 (160)
                #-- Oct18 (month 202) is centered in Nov18 (203)
                if (mon[t] == mon[t-1]) and (mon[t-1] == 118):
                    mon[t] = mon[t-1] + 1
                elif (mon[t] == mon[t-1]) and (mon[t-1] == 121):
                    mon[t-1] = mon[t] - 1
                elif (mon[t] == mon[t-1]) and (mon[t-1] == 160):
                    mon[t] = mon[t-1] + 1
                elif (mon[t] == mon[t-1]) and (mon[t-1] == 203):
                    mon[t-1] = mon[t] - 1
                #-- add to t count
                t += 1
        #-- convert to output variables and truncate if necessary
        tdec = date_conv[:t]
        C20 = C20_input[:t]
        eC20 = eC20_input[:t]
        grace_month = mon[:t]
    else:
        #-- SLR C20 file from PO.DAAC
        with open(os.path.expanduser(SLR_file),'r') as f:
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
            HEADER = not bool(re.match('PRODUCT:+',line))
            #-- add 1 to counter
            count += 1

        #-- number of months within the file
        n_mon = file_lines - count
        date_conv = np.zeros((n_mon))
        C20_input = np.zeros((n_mon))
        eC20_input = np.zeros((n_mon))
        slr_flag = np.zeros((n_mon),dtype=np.bool)
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall('[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            #-- check for empty lines as there are
            #-- slight differences in RL04 TN-05_C20_SLR.txt
            #-- with blanks between the PRODUCT: line and the data
            count = len(line_contents)
            #-- if count is greater than 0
            if (count > 0):
                #-- modified julian date for line
                MJD = np.float(line_contents[0])
                #-- converting from MJD into month, day and year
                YY,MM,DD,hh,mm,ss = convert_julian(MJD+2400000.5,FORMAT='tuple')
                #-- converting from month, day, year into decimal year
                date_conv[t] = convert_calendar_decimal(YY, MM, DAY=DD, HOUR=hh)
                #-- Spherical Harmonic data for line
                C20_input[t] = np.float(line_contents[2])
                eC20_input[t] = np.float(line_contents[4])*1e-10
                #-- line has * flag
                if bool(re.search('\*',line)):
                    slr_flag[t] = True
                #-- add to t count
                t += 1

        #-- truncate for RL04 if necessary
        date_conv = date_conv[:t]
        C20_input = C20_input[:t]
        eC20_input = eC20_input[:t]
        slr_flag = slr_flag[:t]

        #-- GRACE/GRACE-FO month of SLR solutions
        mon = 1 + np.round((date_conv-2002.)*12.)
        #-- number of unique months
        grace_month = np.unique(mon)
        n_uniq = len(grace_month)
        #-- Removing overlapping months to use the data for
        #-- months with limited GRACE accelerometer use
        tdec = np.zeros((n_uniq))
        C20 = np.zeros((n_uniq))
        eC20 = np.zeros((n_uniq))
        #-- New SLR datasets have * flags for the modified GRACE periods
        #-- these GRACE months use half of a prior month in their solution
        #-- this will find these months (marked above with slr_flag)
        for t in range(n_uniq):
            count = np.count_nonzero(mon == grace_month[t])
            #-- there is only one solution for the month
            if (count == 1):
                i = np.nonzero(mon == grace_month[t])
                tdec[t] = date_conv[i]
                C20[t] = C20_input[i]
                eC20[t] = eC20_input[i]
            #-- there is a special solution for the month
            #-- will the solution flagged with slr_flag
            elif (count == 2):
                i = np.nonzero((mon == grace_month[t]) & slr_flag)
                tdec[t] = date_conv[i]
                C20[t] = C20_input[i]
                eC20[t] = eC20_input[i]

    return {'data':C20, 'error':eC20, 'month':grace_month, 'time':tdec}
