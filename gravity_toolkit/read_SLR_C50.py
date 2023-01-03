#!/usr/bin/env python
u"""
read_SLR_C50.py
Written by Yara Mohajerani and Tyler Sutterley (04/2022)

Reads monthly degree 5 zonal spherical harmonic data files from SLR

Dataset distributed by CSR
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr
        gsfc_slr_5x5c61s61.txt

CALLING SEQUENCE:
    SLR_C50 = read_SLR_C50(SLR_file)

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
import warnings
import gravity_toolkit.SLR

# PURPOSE: read Degree 5 zonal data from Satellite Laser Ranging (SLR)
def read_SLR_C50(*args, **kwargs):
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
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.SLR instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.SLR.C50(*args,**kwargs)
