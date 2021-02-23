#!/usr/bin/env python
u"""
read_SLR_C30.py
Written by Yara Mohajerani and Tyler Sutterley (12/2020)

Reads monthly degree 3 zonal spherical harmonic data files from SLR
    https://neptune.gsfc.nasa.gov/gngphys/index.php?section=519

Dataset distributed by NASA PO.DAAC
    https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/gracefo/docs
        TN-14_C30_C30_GSFC_SLR.txt
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt

REFERENCE:
    Loomis, B. D., Rachlin, K. E., and Luthcke, S. B., "Improved Earth
        Oblateness Rate Reveals Increased Ice Sheet Losses and Mass-Driven Sea
        Level Rise", Geophysical Research Letters, 46(12), 6910-6917, 2019.
        https://doi.org/10.1029/2019GL082929

CALLING SEQUENCE:
    SLR_C30 = read_SLR_C30(SLR_file)

INPUTS:
    SLR_file:
        GSFC: TN-14_C30_C30_GSFC_SLR.txt
        CSR: CSR_Monthly_5x5_Gravity_Harmonics.txt
        LARES: C30_LARES_filtered.txt

OUTPUTS:
    data: SLR degree 3 order 0 cosine stokes coefficients (C30)
    error: SLR degree 3 order 0 cosine stokes coefficient error (eC30)
    month: GRACE/GRACE-FO month of measurement (Apr. 2002 = 004)
    time: date of SLR measurement

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    C30_MEAN: mean C30 to add to LARES C30 anomalies

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_CSR_monthly_6x1.py: reads monthly 5x5 spherical harmonic coefficients

UPDATE HISTORY:
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: using utilities from time module
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 08/2019: new GSFC format with more columns
        add catch to verify input SLR file exists
        added LARES filtered C30 files from John Ries (C30_LARES_filtered.txt)
        add C30 mean (9.5717395773300e-07) to LARES solutions
    Updated 07/2019: added SLR C3,0 files from PO.DAAC (GSFC)
        read CSR monthly 5x5 file and extract C3,0 coefficients
    Written 05/2019
"""
import os
import re
import numpy as np
import gravity_toolkit.time
from gravity_toolkit.read_CSR_monthly_6x1 import read_CSR_monthly_6x1

#-- PURPOSE: read Degree 3 zonal data from Satellite Laser Ranging (SLR)
def read_SLR_C30(SLR_file, HEADER=True, C30_MEAN=9.5717395773300e-07):
    """
    Reads C30 spherical harmonic coefficients from SLR measurements

    Arguments
    ---------
    SLR_file: Satellite Laser Ranging file

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped (default: True)
    C30_MEAN: mean C30 to add to LARES C30 anomalies

    Returns
    -------
    data: SLR degree 3 order 0 cosine stokes coefficients
    error: SLR degree 3 order 0 cosine stokes coefficient error
    month: GRACE/GRACE-FO month of measurement
    time: date of SLR measurement
    """

    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise IOError('SLR file not found in file system')
    #-- output dictionary with input data
    dinput = {}

    if bool(re.search('TN-(14)',SLR_file)):

        #-- SLR C30 RL06 file from PO.DAAC
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
        date_conv = np.zeros((n_mon))
        C30_input = np.zeros((n_mon))
        eC30_input = np.zeros((n_mon))
        mon = np.zeros((n_mon),dtype=np.int)
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall(r'[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            count = len(line_contents)
            #-- only read lines where C30 data exists (don't read NaN lines)
            if (count > 7):
                #-- modified julian date for line
                MJD = np.float(line_contents[0])
                #-- converting from MJD into month, day and year
                YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
                    MJD+2400000.5, FORMAT='tuple')
                #-- converting from month, day, year into decimal year
                date_conv[t] = gravity_toolkit.time.convert_calendar_decimal(
                    YY, MM, day=DD, hour=hh)
                #-- Spherical Harmonic data for line
                C30_input[t] = np.float(line_contents[5])
                eC30_input[t] = np.float(line_contents[7])*1e-10
                #-- GRACE/GRACE-FO month of SLR solutions
                mon[t] = 1 + np.round((date_conv[t]-2002.)*12.)
                #-- add to t count
                t += 1
        #-- verify that there imported C30 solutions
        #-- (TN-14 data format has changed in the past)
        if (t == 0):
            raise Exception('No GSFC C30 data imported')
        #-- convert to output variables and truncate if necessary
        dinput['time'] = date_conv[:t]
        dinput['data'] = C30_input[:t]
        dinput['error'] = eC30_input[:t]
        dinput['month'] = mon[:t]
    elif bool(re.search('C30_LARES',SLR_file)):
        #-- read LARES filtered values
        LARES_input = np.loadtxt(SLR_file,skiprows=1)
        dinput['time'] = LARES_input[:,0].copy()
        #-- convert C30 from anomalies to absolute
        dinput['data'] = 1e-10*LARES_input[:,1] + C30_MEAN
        #-- filtered data does not have errors
        dinput['error'] = np.zeros_like(LARES_input[:,1])
        #-- calculate GRACE/GRACE-FO month
        dinput['month'] = 1 + np.array(12.0*(LARES_input[:,0]-2002.0),dtype='i')
    else:
        #-- CSR 5x5 + 6,1 file from CSR and extract C3,0 coefficients
        Ylms = read_CSR_monthly_6x1(SLR_file, HEADER=True)
        #-- extract dates, C30 harmonics and errors
        dinput['time'] = Ylms['time'].copy()
        dinput['data'] = Ylms['clm'][3,0,:].copy()
        dinput['error'] = Ylms['error']['clm'][3,0,:].copy()
        #-- converting from MJD into month, day and year
        YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
            Ylms['MJD']+2400000.5, FORMAT='tuple')
        #-- calculate GRACE/GRACE-FO month
        dinput['month'] = np.array(12.0*(YY - 2002.) + MM, dtype=np.int)

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])

    #-- return the input C30 data, year-decimal date, and GRACE/GRACE-FO month
    return dinput
