#!/usr/bin/env python
u"""
read_gravis_geocenter.py
Written by Tyler Sutterley (04/2022)

Reads monthly geocenter spherical harmonic data files from
    GFZ GravIS calculated using GRACE/GRACE-FO measurements
    and Ocean Models of degree 1

Dataset distributed by GFZ
    ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
        GRAVIS-2B_GFZOP_GEOCENTER_0002.dat

CALLING SEQUENCE:
    geocenter = read_gravis_geocenter(geocenter_file)

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

REFERENCES:
    Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential
        GSM Coefficients GFZ RL06 (Level-2B Product)."
        V. 0002. GFZ Data Services, (2019).
        https://doi.org/10.5880/GFZ.GRAVIS_06_L2B

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: function deprecated. merged with gravity_toolkit.geocenter
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Written 05/2021
"""
import warnings
import gravity_toolkit.geocenter

#-- PURPOSE: read geocenter data from GFZ GravIS SLR/GRACE solutions
def read_gravis_geocenter(geocenter_file, HEADER=True):
    """
    Reads monthly geocenter spherical harmonic data files from
        GFZ GravIS calculated using GRACE/GRACE-FO measurements
        and Ocean Models of degree 1

    Parameters
    ----------
    geocenter_file: str
        degree 1 file
    HEADER: bool, default True
        File contains header text to be skipped

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
    DEG1 = gravity_toolkit.geocenter().from_gravis(geocenter_file,header=HEADER)
    #-- return the GFZ GravIS geocenter solutions
    return DEG1.to_dict()
