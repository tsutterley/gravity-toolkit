#!/usr/bin/env python
u"""
read_SLR_C30.py
Written by Yara Mohajerani and Tyler Sutterley (04/2022)

Reads monthly degree 3 zonal spherical harmonic data files from SLR

Dataset distributed by NASA PO.DAAC
    https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/gracefo/docs
        TN-14_C30_C30_GSFC_SLR.txt
    ftp://ftp.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GFZ
    ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
        GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat

CALLING SEQUENCE:
    SLR_C30 = read_SLR_C30(SLR_file)

INPUTS:
    SLR_file:
        CSR: CSR_Monthly_5x5_Gravity_Harmonics.txt
        GFZ: GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat
        GSFC: TN-14_C30_C30_GSFC_SLR.txt
        LARES: C30_LARES_filtered.txt

OUTPUTS:
    data: SLR degree 3 order 0 cosine stokes coefficients (C30)
    error: SLR degree 3 order 0 cosine stokes coefficient error (eC30)
    month: GRACE/GRACE-FO month of measurement (April 2002 = 004)
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
    read_SLR_harmonics.py: low-degree spherical harmonic coefficients from SLR

REFERENCES:
    Loomis, Rachlin, and Luthcke, "Improved Earth Oblateness Rate Reveals
        Increased Ice Sheet Losses and Mass-Driven Sea Level Rise",
        Geophysical Research Letters, 46(12), 6910-6917, (2019).
        https://doi.org/10.1029/2019GL082929
    Loomis, Rachlin, Wiese, Landerer, and Luthcke, "Replacing GRACE/GRACE-FO
        C30 with satellite laser ranging: Impacts on Antarctic Ice Sheet
        mass change". Geophysical Research Letters, 47, (2020).
        https://doi.org/10.1029/2019GL085488
    Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential
        GSM Coefficients GFZ RL06 (Level-2B Product)."
        V. 0002. GFZ Data Services, (2019).
        https://doi.org/10.5880/GFZ.GRAVIS_06_L2B

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: added GFZ GravIS GRACE/SLR low degree solutions
        define int/float precision to prevent deprecation warning
    Updated 04/2021: renamed SLR monthly 5x5 function from CSR
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
import warnings
import gravity_toolkit.SLR

# PURPOSE: read Degree 3 zonal data from Satellite Laser Ranging (SLR)
def read_SLR_C30(*args, **kwargs):
    """
    Reads C30 spherical harmonic coefficients from SLR measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    C30_MEAN: float, default 9.5717395773300e-07
        Mean C30 to add to LARES C30 anomalies
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    data: float
        SLR degree 3 order 0 cosine stokes coefficients
    error: float
        SLR degree 3 order 0 cosine stokes coefficient error
    month: int
        GRACE/GRACE-FO month of measurement
    time: float
        date of SLR measurement
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.SLR instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.SLR.C30(*args,**kwargs)
