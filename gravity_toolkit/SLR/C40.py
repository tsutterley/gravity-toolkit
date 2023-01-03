#!/usr/bin/env python
u"""
C40.py
Written by Tyler Sutterley (01/2023)

Reads monthly degree 4 zonal spherical harmonic data files from SLR

Dataset distributed by CSR
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr
        gsfc_slr_5x5c61s61.txt

CALLING SEQUENCE:
    SLR_C40 = gravity_toolkit.SLR.C40(SLR_file)

INPUTS:
    SLR_file:
        GSFC: gsfc_slr_5x5c61s61.txt
        CSR: CSR_Monthly_5x5_Gravity_Harmonics.txt

OUTPUTS:
    data: SLR degree 4 order 0 cosine stokes coefficients (C40)
    error: SLR degree 4 order 0 cosine stokes coefficient error (eC40)
    month: GRACE/GRACE-FO month of measurement (April 2002 = 004)
    time: date of SLR measurement

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    C40_MEAN: mean C40 to add to LARES C40 anomalies
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
    Written 09/2022
"""
import os
import re
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.read_SLR_harmonics

# PURPOSE: read Degree 4 zonal data from Satellite Laser Ranging (SLR)
def C40(SLR_file, C40_MEAN=0.0, DATE=None, **kwargs):
    """
    Reads C40 spherical harmonic coefficients from SLR measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    C40_MEAN: float, default 0.0
        Mean C40 to add to LARES C40 anomalies
    DATE: float or NoneType, default None
        Mid-point of monthly solution for calculating 28-day arc averages

    Returns
    -------
    data: float
        SLR degree 4 order 0 cosine stokes coefficients
    error: float
        SLR degree 4 order 0 cosine stokes coefficient error
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

    if bool(re.search(r'gsfc_slr_5x5c61s61',SLR_file,re.I)):
        # read 5x5 + 6,1 file from GSFC and extract coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        # calculate 28-day moving-average solution from 7-day arcs
        dinput.update(gravity_toolkit.convert_weekly(Ylms['time'],
            Ylms['clm'][4,0,:], DATE=DATE, NEIGHBORS=28))
        # no estimated spherical harmonic errors
        dinput['error'] = np.zeros_like(DATE,dtype='f8')
    elif bool(re.search(r'C40_LARES',SLR_file,re.I)):
        # read LARES filtered values
        LARES_input = np.loadtxt(SLR_file,skiprows=1)
        dinput['time'] = LARES_input[:,0].copy()
        # convert C40 from anomalies to absolute
        dinput['data'] = 1e-10*LARES_input[:,1] + C40_MEAN
        # filtered data does not have errors
        dinput['error'] = np.zeros_like(LARES_input[:,1])
        # calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    else:
        # read 5x5 + 6,1 file from CSR and extract C4,0 coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        # extract dates, C40 harmonics and errors
        dinput['time'] = Ylms['time'].copy()
        dinput['data'] = Ylms['clm'][4,0,:].copy()
        dinput['error'] = Ylms['error']['clm'][4,0,:].copy()
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

    # return the SLR-derived degree 4 zonal solutions
    return dinput
