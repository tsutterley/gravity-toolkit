#!/usr/bin/env python
u"""
C50.py
Written by Yara Mohajerani and Tyler Sutterley (01/2023)

Reads monthly degree 5 zonal spherical harmonic data files from SLR

Dataset distributed by CSR
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr
        gsfc_slr_5x5c61s61.txt

CALLING SEQUENCE:
    SLR_C50 = gravity_toolkit.SLR.C50(SLR_file)

INPUTS:
    SLR_file:
        GSFC: gsfc_slr_5x5c61s61.txt
        CSR: CSR_Monthly_5x5_Gravity_Harmonics.txt

OUTPUTS:
    data: SLR degree 5 order 0 cosine stokes coefficients (C50)
    error: SLR degree 5 order 0 cosine stokes coefficient error (eC50)
    month: GRACE/GRACE-FO month of measurement (April 2002 = 004)
    time: date of SLR measurement

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    C50_MEAN: mean C50 to add to LARES C50 anomalies
    DATE: mid-point of monthly solution for calculating 28-day arc averages

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_SLR_harmonics.py: low-degree spherical harmonic coefficients from SLR

UPDATE HISTORY:
    Updated 01/2023: refactored satellite laser ranging read functions
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: use function for converting from 7-day arcs
    Updated 11/2021: reader for new weekly 5x5+6,1 fields from NASA GSFC
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
import gravity_toolkit.read_SLR_harmonics

# PURPOSE: read Degree 5 zonal data from Satellite Laser Ranging (SLR)
def C50(SLR_file, C50_MEAN=0.0, DATE=None, HEADER=True):
    """
    Reads C50 spherical harmonic coefficients from SLR measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    C50_MEAN: float, default 0.0
        Mean C50 to add to LARES C50 anomalies
    DATE: float or NoneType, default None
        Mid-point of monthly solution for calculating 28-day arc averages
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    data: float
        SLR degree 5 order 0 cosine stokes coefficients
    error: float
        SLR degree 5 order 0 cosine stokes coefficient error
    month: int
        GRACE/GRACE-FO month of measurement
    time: float
        date of SLR measurement
    """

    # check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')
    # output dictionary with input data
    dinput = {}

    if bool(re.search(r'GSFC_SLR_C(20)_C(30)_C(50)',SLR_file,re.I)):

        # SLR C50 RL06 file from GSFC
        with open(os.path.expanduser(SLR_file), mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        # number of lines contained in the file
        file_lines = len(file_contents)

        # counts the number of lines in the header
        count = 0
        # Reading over header text
        while HEADER:
            # file line at count
            line = file_contents[count]
            # find PRODUCT: within line to set HEADER flag to False when found
            HEADER = not bool(re.match(r'Product:+',line))
            # add 1 to counter
            count += 1

        # number of months within the file
        n_mon = file_lines - count
        # date and GRACE/GRACE-FO month
        dinput['time'] = np.zeros((n_mon))
        dinput['month'] = np.zeros((n_mon),dtype=int)
        # monthly spherical harmonic replacement solutions
        dinput['data'] = np.zeros((n_mon))
        # monthly spherical harmonic formal standard deviations
        dinput['error'] = np.zeros((n_mon))
        # time count
        t = 0
        # for every other line:
        for line in file_contents[count:]:
            # find numerical instances in line including exponents,
            # decimal points and negatives
            line_contents = re.findall(r'[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            count = len(line_contents)
            # only read lines where C50 data exists (don't read NaN lines)
            if (count > 7):
                # modified julian date for line
                MJD = np.float64(line_contents[0])
                # converting from MJD into month, day and year
                YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
                    MJD+2400000.5, format='tuple')
                # converting from month, day, year into decimal year
                dinput['time'][t] = gravity_toolkit.time.convert_calendar_decimal(
                    YY, MM, day=DD, hour=hh)
                # Spherical Harmonic data for line
                dinput['data'][t] = np.float64(line_contents[10])
                dinput['error'][t] = np.float64(line_contents[12])*1e-10
                # GRACE/GRACE-FO month of SLR solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                # add to t count
                t += 1
        # verify that there imported C50 solutions
        if (t == 0):
            raise Exception('No GSFC C50 data imported')
        # truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]
    elif bool(re.search(r'gsfc_slr_5x5c61s61',SLR_file,re.I)):
        # read 5x5 + 6,1 file from GSFC and extract coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        # calculate 28-day moving-average solution from 7-day arcs
        dinput.update(gravity_toolkit.convert_weekly(Ylms['time'],
            Ylms['clm'][5,0,:], DATE=DATE, NEIGHBORS=28))
        # no estimated spherical harmonic errors
        dinput['error'] = np.zeros_like(DATE,dtype='f8')
    elif bool(re.search(r'C50_LARES',SLR_file,re.I)):
        # read LARES filtered values
        LARES_input = np.loadtxt(SLR_file,skiprows=1)
        dinput['time'] = LARES_input[:,0].copy()
        # convert C50 from anomalies to absolute
        dinput['data'] = 1e-10*LARES_input[:,1] + C50_MEAN
        # filtered data does not have errors
        dinput['error'] = np.zeros_like(LARES_input[:,1])
        # calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    else:
        # read 5x5 + 6,1 file from CSR and extract C5,0 coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        # extract dates, C50 harmonics and errors
        dinput['time'] = Ylms['time'].copy()
        dinput['data'] = Ylms['clm'][5,0,:].copy()
        dinput['error'] = Ylms['error']['clm'][5,0,:].copy()
        # converting from MJD into month, day and year
        YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(
            Ylms['MJD']+2400000.5, format='tuple')
        # calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(YY,MM)

    # The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    # Accelerometer shutoffs make the relation between month number
    # and date more complicated as days from other months are used
    # For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    # For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    # For all: May 2015 (161) is centered in Apr 2015 (160)
    # For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])

    # return the SLR-derived degree 5 zonal solutions
    return dinput
