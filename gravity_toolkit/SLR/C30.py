#!/usr/bin/env python
u"""
C30.py
Written by Yara Mohajerani and Tyler Sutterley (01/2023)

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
    SLR_C30 = gravity_toolkit.SLR.C30(SLR_file)

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
    Updated 01/2023: refactored satellite laser ranging read functions
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
import os
import re
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.read_SLR_harmonics

# PURPOSE: read Degree 3 zonal data from Satellite Laser Ranging (SLR)
def C30(SLR_file, C30_MEAN=9.5717395773300e-07, HEADER=True):
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

    # check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')
    # output dictionary with input data
    dinput = {}

    if bool(re.search(r'TN-(14)',SLR_file,re.I)):

        # SLR C30 RL06 file from PO.DAAC produced by GSFC
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
            # only read lines where C30 data exists (don't read NaN lines)
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
                dinput['data'][t] = np.float64(line_contents[5])
                dinput['error'][t] = np.float64(line_contents[7])*1e-10
                # GRACE/GRACE-FO month of SLR solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                # add to t count
                t += 1
        # verify that there imported C30 solutions
        # (TN-14 data format has changed in the past)
        if (t == 0):
            raise Exception('No GSFC C30 data imported')
        # truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]
    elif bool(re.search(r'C30_LARES',SLR_file,re.I)):
        # read LARES filtered values
        LARES_input = np.loadtxt(SLR_file,skiprows=1)
        dinput['time'] = LARES_input[:,0].copy()
        # convert C30 from anomalies to absolute
        dinput['data'] = 1e-10*LARES_input[:,1] + C30_MEAN
        # filtered data does not have errors
        dinput['error'] = np.zeros_like(LARES_input[:,1])
        # calculate GRACE/GRACE-FO month
        dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    elif bool(re.search(r'GRAVIS-2B_GFZOP',SLR_file,re.I)):
        # Combined GRACE/SLR solution file produced by GFZ
        # Column  1: MJD of BEGINNING of solution data span
        # Column  2: Year and fraction of year of BEGINNING of solution span
        # Column  6: Replacement C(3,0)
        # Column  7: Replacement C(3,0) - mean C(3,0) (1.0E-10)
        # Column  8: C(3,0) formal standard deviation (1.0E-12)
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
            HEADER = not bool(re.match(r'PRODUCT:+',line))
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
            # check for empty lines
            if (count > 0):
                # reading decimal year for start of span
                dinput['time'][t] = np.float64(line_contents[1])
                # Spherical Harmonic data for line
                dinput['data'][t] = np.float64(line_contents[5])
                dinput['error'][t] = np.float64(line_contents[7])*1e-10
                # GRACE/GRACE-FO month of SLR solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                # add to t count
                t += 1
        # truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]
    else:
        # CSR 5x5 + 6,1 file from CSR and extract C3,0 coefficients
        Ylms = gravity_toolkit.read_SLR_harmonics(SLR_file, HEADER=True)
        # extract dates, C30 harmonics and errors
        dinput['time'] = Ylms['time'].copy()
        dinput['data'] = Ylms['clm'][3,0,:].copy()
        dinput['error'] = Ylms['error']['clm'][3,0,:].copy()
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

    # return the SLR-derived degree 3 zonal solutions
    return dinput
