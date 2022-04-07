#!/usr/bin/env python
u"""
read_tellus_geocenter.py
Written by Tyler Sutterley (04/2022)

Reads monthly geocenter spherical harmonic data files from GRACE Tellus
    Technical Notes (TN-13) calculated using GRACE/GRACE-FO measurements and
    Ocean Models of Degree 1

Datasets distributed by NASA PO.DAAC
https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1

Swenson, S., D. Chambers, and J. Wahr, "Estimating geocenter variations
    from a combination of GRACE and ocean model output", J. Geophys. Res.,
    113(B08410), 2008. doi:10.1029/2007JB005338

Sun, Y., R. Riva, and P. Ditmar, "Observed changes in the Earth's dynamic
    oblateness from GRACE data and geophysical models", J. Geodesy.,
    90(1), 81-89, 2016. doi:10.1007/s00190-015-0852-y

Sun, Y., R. Riva, and P. Ditmar, "Optimizing estimates of annual
    variations and trends in geocenter motion and J2 from a combination
    of GRACE data and geophysical models", J. Geophys. Res. Solid Earth,
    121, 2016. doi:10.1002/2016JB013073

CALLING SEQUENCE:
    geocenter = read_tellus_geocenter(geocenter_file)

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
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: function deprecated. merged with gravity_toolkit.geocenter
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use file not found exceptions
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: using utilities from time module
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

import warnings
import gravity_toolkit.geocenter

#-- PURPOSE: read geocenter data from PO.DAAC
def read_tellus_geocenter(geocenter_file, HEADER=True, JPL=False):
    """
    Reads monthly geocenter files computed by JPL Tellus using
    GRACE/GRACE-FO measurements and Ocean Models of degree 1

    Parameters
    ----------
    geocenter_file: str
        degree 1 file
    HEADER: bool, default True
        file contains header text to be skipped
    JPL: bool, default True
        Use JPL TN-13 geocenter files with self-attraction and loading

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
    warnings.warn("Deprecated. Please use gravity_toolkit.geocenter instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    DEG1 = gravity_toolkit.geocenter().from_tellus(geocenter_file,
        header=HEADER, JPL=JPL)
    #-- return the JPL GRACE Tellus geocenter solutions
    return DEG1.to_dict()
