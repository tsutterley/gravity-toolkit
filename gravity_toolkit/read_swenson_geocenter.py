#!/usr/bin/env python
u"""
read_swenson_geocenter.py
Written by Tyler Sutterley (04/2022)

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
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: function deprecated. merged with gravity_toolkit.geocenter
    Updated 09/2021: use functions for converting to and from GRACE months
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
import warnings
import gravity_toolkit.geocenter

#-- PURPOSE: read geocenter data from Sean Swenson
def read_swenson_geocenter(geocenter_file, HEADER=True):
    """
    Reads monthly geocenter files computed by Sean Swenson using
    GRACE/GRACE-FO measurements and Ocean Models of degree 1

    Parameters
    ----------
    geocenter_file: str
        degree 1 file
    HEADER: bool, default True
        file contains header text to be skipped

    Returns
    -------
    C10: float
        cosine d1/o0 spherical harmonic coefficients
    C11: float
        cosine d1/o1 spherical harmonic coefficients
    S11: float
        sine d1/o1 spherical harmonic coefficients
    month: int
        GRACE/GRACE-FO month
    time: float
        date of each month in year-decimal
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.geocenter instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    DEG1 = gravity_toolkit.geocenter().from_swenson(geocenter_file,
        header=HEADER)
    #-- return the geocenter solutions from Sean Swenson
    return DEG1.to_dict()
