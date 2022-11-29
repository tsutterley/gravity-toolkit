#!/usr/bin/env python
u"""
read_SLR_geocenter.py
Written by Tyler Sutterley (04/2022)

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
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: function deprecated. merged with gravity_toolkit.geocenter
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use file not found exceptions
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

import warnings
import gravity_toolkit.geocenter

# PURPOSE: read geocenter data from Satellite Laser Ranging (SLR)
def read_SLR_geocenter(geocenter_file, RADIUS=None, HEADER=0,
    COLUMNS=['time', 'X', 'Y', 'Z', 'X_sigma', 'Y_sigma', 'Z_sigma']):
    """
    Reads monthly geocenter files from satellite laser ranging

    Parameters
    ----------
    geocenter_file: str
        Satellite Laser Ranging file
    RADIUS: float or NoneType, default None
        Earth's radius for calculating spherical harmonics from SLR data
    HEADER: int, default 0
        Rows of data to skip when importing data
    COLUMNS: list
        Column names of ascii file

            - ``'time'``: date in decimal-years
            - ``'X'``: X-component of geocenter variation
            - ``'Y'``: Y-component of geocenter variation
            - ``'Z'``: Z-component of geocenter variation
            - ``'X_sigma'``: X-component uncertainty
            - ``'Y_sigma'``: Y-component uncertainty
            - ``'Z_sigma'``: Z-component uncertainty

    Returns
    -------
    C10: float
        cosine d1/o0 spherical harmonic coefficients
    C11: float
        cosine d1/o1 spherical harmonic coefficients
    S11: float
        sine d1/o1 spherical harmonic coefficients
    eC10: float
        cosine d1/o0 spherical harmonic coefficient error
    eC11: float
        cosine d1/o1 spherical harmonic coefficient error
    eS11: float
        sine d1/o1 spherical harmonic coefficient error
    month: int
        GRACE/GRACE-FO month
    time: float
        date of each month in year-decimal
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.harmonics instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    DEG1 = gravity_toolkit.geocenter(radius=RADIUS).from_SLR(geocenter_file,
        AOD=False, header=HEADER, columns=COLUMNS)
    # return the SLR-derived geocenter solutions
    return DEG1.to_dict()


# special function for outputting AOD corrected SLR geocenter values
# need to run aod1b_geocenter.py to calculate the monthly geocenter dealiasing
def aod_corrected_SLR_geocenter(geocenter_file, DREL, RADIUS=None, HEADER=0,
    COLUMNS=[]):
    """
    Reads monthly geocenter files from satellite laser ranging corrected
    for non-tidal ocean and atmospheric variation

    Parameters
    ----------
    geocenter_file: Satellite Laser Ranging file
    DREL: GRACE/GRACE-FO/Swarm data release

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
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.geocenter instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    DEG1 = gravity_toolkit.geocenter(radius=RADIUS).from_SLR(geocenter_file,
        AOD=True, release=DREL, header=HEADER, columns=COLUMNS)
    # return the SLR-derived geocenter solutions
    return DEG1.to_dict()
