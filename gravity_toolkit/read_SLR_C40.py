#!/usr/bin/env python
u"""
read_SLR_C40.py
Written by Tyler Sutterley (09/2022)

Reads monthly degree 4 zonal spherical harmonic data files from SLR

Dataset distributed by CSR
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr
        gsfc_slr_5x5c61s61.txt

CALLING SEQUENCE:
    SLR_C40 = read_SLR_C40(SLR_file)

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
    Written 09/2022
"""
import warnings
import gravity_toolkit.SLR

# PURPOSE: read Degree 4 zonal data from Satellite Laser Ranging (SLR)
def read_SLR_C40(*args, **kwargs):
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
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.SLR instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.SLR.C40(*args,**kwargs)
