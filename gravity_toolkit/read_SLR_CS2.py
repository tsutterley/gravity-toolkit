#!/usr/bin/env python
u"""
read_SLR_CS2.py
Written by Hugo Lecomte and Tyler Sutterley (04/2022)

Reads monthly degree 2,m (figure axis and azimuthal dependence)
    spherical harmonic data files from satellite laser ranging (SLR)

Dataset distributed by CSR
    http://download.csr.utexas.edu/pub/slr/degree_2/
        C21_S21_RL06.txt or C22_S22_RL06.txt
Dataset distributed by GFZ
    ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
        GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr

CALLING SEQUENCE:
    SLR_2m = read_SLR_CS2(SLR_file)

INPUTS:
    SLR_file:
        CSR 2,1: C21_S21_RL06.txt
        CSR 2,2: C22_S22_RL06.txt
        GFZ: GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat
        GSFC: GSFC_C21_S21.txt

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    DATE: mid-point of monthly solution for calculating 28-day arc averages

OUTPUTS:
    C2m: SLR degree 2 order m cosine stokes coefficients
    S2m: SLR degree 2 order m sine stokes coefficients
    eC2m: SLR degree 2 order m cosine stokes coefficient error
    eS2m: SLR degree 2 order m sine stokes coefficient error
    month: GRACE/GRACE-FO month of measurement (Apr. 2002 = 004)
    time: date of SLR measurement

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_SLR_harmonics.py: low-degree spherical harmonic coefficients from SLR

REFERENCES:
    Cheng et al., " Variations of the Earth's figure axis from satellite
        laser ranging and GRACE", Journal of Geophysical Research,
        116, B01409, (2011). https://doi.org/10.1029/2010JB000850
    Dahle et al., "The GFZ GRACE RL06 Monthly Gravity Field Time Series:
        Processing Details, and Quality Assessment", Remote Sensing,
        11(18), 2116, (2019). https://doi.org/10.3390/rs11182116
    Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential
        GSM Coefficients GFZ RL06 (Level-2B Product)."
        V. 0002. GFZ Data Services, (2019).
        https://doi.org/10.5880/GFZ.GRAVIS_06_L2B
    Chen el al., "Assessment of degree-2 order-1 gravitational changes
        from GRACE and GRACE Follow-on, Earth rotation, satellite laser
        ranging, and models", Journal of Geodesy, 95(38), (2021).
        https://doi.org/10.1007/s00190-021-01492-x

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 11/2021: reader for new weekly 5x5+6,1 fields from NASA GSFC
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 08/2021: output empty spherical harmonic errors for GSFC
    Updated 06/2021: added GSFC 7-day SLR figure axis solutions
    Updated 05/2021: added GFZ GravIS GRACE/SLR low degree solutions
    Updated 04/2021: use adjust_months function to fix special months cases
    Written 11/2020
"""
import warnings
import gravity_toolkit.SLR

# PURPOSE: read Degree 2,m data from Satellite Laser Ranging (SLR)
def read_SLR_CS2(*args, **kwargs):
    """
    Reads CS2,m spherical harmonic coefficients from SLR measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    ORDER: int, default 1
        Spherical harmonic order to extract from low-degree fields
    DATE: float or NoneType, default None
        Mid-point of monthly solution for calculating 28-day arc averages
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    C2m: float
        SLR degree 2 order m cosine stokes coefficients
    S2m: float
        SLR degree 2 order m sine stokes coefficients
    eC2m: float
        SLR degree 2 order m cosine stokes coefficient error
    eS2m: float
        SLR degree 2 order m sine stokes coefficient error
    month: int
        GRACE/GRACE-FO month of measurement
    time: float
        date of SLR measurement
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.SLR instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.SLR.CS2(*args,**kwargs)
