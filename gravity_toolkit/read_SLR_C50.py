#!/usr/bin/env python
u"""
read_SLR_C50.py
Written by Yara Mohajerani and Tyler Sutterley (09/2021)

This program reads in C50 spherical harmonic coefficients from SLR measurements
    https://neptune.gsfc.nasa.gov/gngphys/index.php?section=519

Dataset distributed by NASA PO.DAAC
    https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/gracefo/docs
    C50 file sent by Bryan Loomis
        GSFC_SLR_C20_C30_C50_GSM_replacement.txt
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt

CALLING SEQUENCE:
    SLR_C50 = read_SLR_C50(SLR_file)

INPUTS:
    SLR_file:
        GSFC: GSFC_SLR_C20_C30_C50_GSM_replacement.txt
        CSR: CSR_Monthly_5x5_Gravity_Harmonics.txt

OUTPUTS:
    data: SLR degree 5 order 0 cosine stokes coefficients (C50)
    error: SLR degree 5 order 0 cosine stokes coefficient error (eC50)
    month: GRACE/GRACE-FO month of measurement (April 2002 = 004)
    time: date of SLR measurement

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    C50_MEAN: mean C50 to add to LARES C50 anomalies

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_SLR_monthly_6x1.py: reads monthly 5x5 spherical harmonic coefficients

UPDATE HISTORY:
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: simplified program similar to other SLR readers
        define int/float precision to prevent deprecation warning
    Updated 04/2021: using utilities from time module
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Written 11/2019
"""
import os
import re
import numpy as np
import gravity_toolkit.time
from gravity_toolkit.read_SLR_monthly_6x1 import read_SLR_monthly_6x1

#-- PURPOSE: read Degree 5 zonal data from Satellite Laser Ranging (SLR)
def read_SLR_C50(SLR_file, HEADER=True, C50_MEAN=0.):
    """
    Reads C50 spherical harmonic coefficients from SLR measurements

    Arguments
    ---------
    SLR_file: Satellite Laser Ranging file

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped (default: True)
    C50_MEAN: mean C50 to add to LARES C50 anomalies

    Returns
    -------
    data: SLR degree 5 order 0 cosine stokes coefficients
    error: SLR degree 5 order 0 cosine stokes coefficient error
    month: GRACE/GRACE-FO month of measurement
    time: date of SLR measurement
    """

    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')
    #-- output dictionary with input data
    dinput = {}

    if bool(re.search(r'GSFC_SLR_C(20)_C(30)_C(50)',SLR_file)):

        #-- SLR C50 RL06 file from GSFC
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
            HEADER = not bool(re.match(r'Product:+',line))
            #-- add 1 to counter
            count += 1

        #-- number of months within the file
        n_mon = file_lines - count
        #-- date and GRACE/GRACE-FO month
        dinput['time'] = np.zeros((n_mon))
        dinput['month'] = np.zeros((n_mon),dtype=int)
        #-- monthly spherical harmonic replacement solutions
        dinput['data'] = np.zeros((n_mon))
        #-- monthly spherical harmonic formal standard deviations
        dinput['error'] = np.zeros((n_mon))
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall(r'[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            count = len(line_contents)
            #-- only read lines where C50 data exists (don't read NaN lines)
            if (count > 7):
                #-- modified julian date for line
                MJD = np.float64(line_contents[0])
                #-- converting from MJD into month, day and year
                YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
                    MJD+2400000.5, FORMAT='tuple')
                #-- converting from month, day, year into decimal year
                dinput['time'][t] = gravity_toolkit.time.convert_calendar_decimal(
                    YY, MM, day=DD, hour=hh)
                #-- Spherical Harmonic data for line
                dinput['data'][t] = np.float64(line_contents[10])
                dinput['error'][t] = np.float64(line_contents[12])*1e-10
                #-- GRACE/GRACE-FO month of SLR solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                #-- add to t count
                t += 1
        #-- verify that there imported C50 solutions
        if (t == 0):
            raise Exception('No GSFC C50 data imported')
        #-- truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]
    elif bool(re.search(r'C50_LARES',SLR_file)):
        #-- read LARES filtered values
        LARES_input = np.loadtxt(SLR_file,skiprows=1)
        dinput['time'] = LARES_input[:,0].copy()
        #-- convert C50 from anomalies to absolute
        dinput['data'] = 1e-10*LARES_input[:,1] + C50_MEAN
        #-- filtered data does not have errors
        dinput['error'] = np.zeros_like(LARES_input[:,1])
        #-- calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    else:
        #-- CSR 5x5 + 6,1 file from CSR and extract C5,0 coefficients
        Ylms = read_SLR_monthly_6x1(SLR_file, HEADER=True)
        #-- extract dates, C50 harmonics and errors
        dinput['time'] = Ylms['time'].copy()
        dinput['data'] = Ylms['clm'][5,0,:].copy()
        dinput['error'] = Ylms['error']['clm'][5,0,:].copy()
        #-- converting from MJD into month, day and year
        YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
            Ylms['MJD']+2400000.5, FORMAT='tuple')
        #-- calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(YY,MM)

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])

    #-- return the SLR-derived degree 5 zonal solutions
    return dinput
