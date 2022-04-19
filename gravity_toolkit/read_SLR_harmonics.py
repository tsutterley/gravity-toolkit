#!/usr/bin/env python
u"""
read_SLR_harmonics.py
Written by Tyler Sutterley (04/2022)

Reads in low-degree spherical harmonic coefficients calculated from
    Satellite Laser Ranging (SLR) measurements

Dataset distributed by UTCSR
    ftp://ftp.csr.utexas.edu/outgoing/cheng/slrgeo.5d561_187_naod
    http://download.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt
Dataset distributed by GSFC
    https://earth.gsfc.nasa.gov/geo/data/slr

OPTIONS:
    SCALE: scale factor for converting to fully-normalized spherical harmonics
    HEADER: file contains header text to be skipped (default: True)

OUTPUTS:
    clm: Cosine spherical harmonic coefficients
    slm: Sine spherical harmonic coefficients
    error/clm: Cosine spherical harmonic coefficient uncertainty
    error/slm: Sine spherical harmonic coefficients uncertainty
    MJD: output date as Modified Julian Day
    time: output date in year-decimal

REFERENCES:
    Cheng, Ries, and Tapley, "Variations of the Earth's Figure Axis from
        Satellite Laser Ranging and GRACE", Journal of Geophysical Research,
        116, B01409, (2011). https://doi.org/10.1029/2010JB000850
    Loomis, Rachlin, Wiese, Landerer, and Luthcke, "Replacing GRACE/GRACE‐FO
        C30 with satellite laser ranging: Impacts on Antarctic Ice Sheet
        mass change". Geophysical Research Letters, 47, (2020).
        https://doi.org/10.1029/2019GL085488

NOTES:
    Degree-1 variations are not included with the new 5x5 product
    Thanks to Hugo Lecomte (University of Strasbourg) for pointing out the
        change in the CSR SLR file format (04/2021)

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
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: added function for converting from 7-day arcs
    Updated 11/2021: renamed module. added reader for GSFC weekly 5x5 fields
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: renamed module. new SLR 5x5 format from CSR (see notes)
        add file read checks (for mean harmonics header and arc number)
        use file not found exceptions for test
    Updated 12/2020: using utilities from time module
    Updated 07/2020: added function docstrings
    Updated 07/2019: following new format with mean field in header and no C6,0
    Updated 10/2018: using future division for python3 Compatibility
    Updated 10/2017: include the 6,0 and 6,1 coefficients in output Ylms
    Written 10/2017
"""
from __future__ import print_function, division

import os
import re
import numpy as np
import gravity_toolkit.time

#-- PURPOSE: wrapper function for calling individual readers
def read_SLR_harmonics(SLR_file, **kwargs):
    """
    Wrapper function for reading spherical harmonic coefficients
    from Satellite Laser Ranging (SLR) measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file
    **kwargs: dict
        keyword arguments for input readers
    """
    if bool(re.search(r'gsfc_slr_5x5c61s61',SLR_file,re.I)):
        return read_GSFC_weekly_6x1(SLR_file, **kwargs)
    elif bool(re.search(r'CSR_Monthly_5x5_Gravity_Harmonics',SLR_file,re.I)):
        return read_CSR_monthly_6x1(SLR_file, **kwargs)
    else:
        raise Exception('Unknown SLR file format {0}'.format(SLR_file))

#-- PURPOSE: read monthly degree harmonic data from Satellite Laser Ranging (SLR)
def read_CSR_monthly_6x1(SLR_file, SCALE=1e-10, HEADER=True):
    """
    Reads in monthly low degree and order spherical harmonic coefficients
    from Satellite Laser Ranging (SLR) measurements

    Parameters
    ----------
    SLR_file: str
        Satellite Laser Ranging file from CSR
    SCALE: float, default 1e-10
        Scale factor for converting to fully-normalized spherical harmonics
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    clm: float
        Cosine spherical harmonic coefficients
    slm: float
        Sine spherical harmonic coefficients
    error/clm: float
        Cosine spherical harmonic coefficient uncertainty
    error/slm: float
        Sine spherical harmonic coefficients uncertainty
    MJD: float
        output date as Modified Julian Day
    time: float
        output date in year-decimal

    References
    ----------
    .. [Cheng2010] M. Cheng, J. C. Ries, and B. D. Tapley,
        "Variations of the Earth's figure axis from satellite laser ranging
        and GRACE", *Journal of Geophysical Research*, 116(B01409), (2010).
        `doi: 10.1029/2010JB000850 <https://doi.org/10.1029/2010JB000850>`_
    """
    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')

    #-- read the file and get contents
    with open(os.path.expanduser(SLR_file), mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()
    file_lines = len(file_contents)

    #-- spherical harmonic degree range (5x5 with 6,1)
    #-- new 5x5 fields no longer include geocenter components
    LMIN = 2
    LMAX = 6
    n_harm = (LMAX**2 + 3*LMAX - LMIN**2 - LMIN)//2 - 5

    #-- counts the number of lines in the header
    count = 0
    indice = None
    #-- Reading over header text
    while HEADER:
        #-- file line at count
        line = file_contents[count]
        #-- find end within line to set HEADER flag to False when found
        HEADER = not bool(re.match(r'end\sof\sheader',line))
        if bool(re.match(80*r'=',line)):
            indice = count + 1
        #-- add 1 to counter
        count += 1

    #-- verify that mean field header indice was found
    if not indice:
        raise Exception('Mean field header not found')

    #-- number of dates within the file
    n_dates = (file_lines - count)//(n_harm + 1)

    #-- read mean fields from the header
    mean_Ylms = {}
    mean_Ylm_error = {}
    mean_Ylms['clm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylms['slm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylm_error['clm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylm_error['slm'] = np.zeros((LMAX+1,LMAX+1))
    for i in range(n_harm):
        #-- split the line into individual components
        line = file_contents[indice+i].split()
        #-- degree and order for the line
        l1 = np.int64(line[0])
        m1 = np.int64(line[1])
        #-- fill mean field Ylms
        mean_Ylms['clm'][l1,m1] = np.float64(line[2].replace('D','E'))
        mean_Ylms['slm'][l1,m1] = np.float64(line[3].replace('D','E'))
        mean_Ylm_error['clm'][l1,m1] = np.float64(line[4].replace('D','E'))
        mean_Ylm_error['slm'][l1,m1] = np.float64(line[5].replace('D','E'))

    #-- output spherical harmonic fields
    Ylms = {}
    Ylms['error'] = {}
    Ylms['MJD'] = np.zeros((n_dates))
    Ylms['time'] = np.zeros((n_dates))
    Ylms['clm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylms['slm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylms['error']['clm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylms['error']['slm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    #-- input spherical harmonic anomalies and errors
    Ylm_anomalies = {}
    Ylm_anomaly_error = {}
    Ylm_anomalies['clm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylm_anomalies['slm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylm_anomaly_error['clm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylm_anomaly_error['slm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    #-- for each date
    for d in range(n_dates):
        #-- split the date line into individual components
        line_contents = file_contents[count].split()
        #-- verify arc number from iteration and file
        IARC = int(line_contents[0])
        assert (IARC == (d+1))
        #-- modified Julian date of the middle of the month
        Ylms['MJD'][d] = np.mean(np.array(line_contents[5:7],dtype=np.float64))
        #-- date of the mid-point of the arc given in years
        YY,MM = np.array(line_contents[3:5])
        Ylms['time'][d] = gravity_toolkit.time.convert_calendar_decimal(YY,MM)
        #-- add 1 to counter
        count += 1

        #-- read the anomaly field
        for i in range(n_harm):
            #-- split the line into individual components
            line = file_contents[count].split()
            #-- degree and order for the line
            l1 = np.int64(line[0])
            m1 = np.int64(line[1])
            #-- fill anomaly field Ylms and rescale to output
            Ylm_anomalies['clm'][l1,m1,d] = np.float64(line[2])*SCALE
            Ylm_anomalies['slm'][l1,m1,d] = np.float64(line[3])*SCALE
            Ylm_anomaly_error['clm'][l1,m1,d] = np.float64(line[6])*SCALE
            Ylm_anomaly_error['slm'][l1,m1,d] = np.float64(line[7])*SCALE
            #-- add 1 to counter
            count += 1

        #-- calculate full coefficients and full errors
        Ylms['clm'][:,:,d] = Ylm_anomalies['clm'][:,:,d] + mean_Ylms['clm'][:,:]
        Ylms['slm'][:,:,d] = Ylm_anomalies['slm'][:,:,d] + mean_Ylms['slm'][:,:]
        Ylms['error']['clm'][:,:,d]=np.sqrt(Ylm_anomaly_error['clm'][:,:,d]**2 +
            mean_Ylm_error['clm'][:,:]**2)
        Ylms['error']['slm'][:,:,d]=np.sqrt(Ylm_anomaly_error['slm'][:,:,d]**2 +
            mean_Ylm_error['slm'][:,:]**2)

    #-- return spherical harmonic fields and date information
    return Ylms

#-- PURPOSE: read weekly degree harmonic data from Satellite Laser Ranging (SLR)
def read_GSFC_weekly_6x1(SLR_file, SCALE=1.0, HEADER=True):
    """
    Reads weekly 5x5 spherical harmonic coefficients with 1 coefficient from
        degree 6 calculated from satellite laser ranging measurements

    Parameters
    ----------
    SLR_file: str
        Satellite laser ranging file from GSFC
    SCALE: float, default 1.0
        Scale factor for converting to fully-normalized spherical harmonics
    HEADER: bool, default True
        File contains header text to be skipped

    Returns
    -------
    clm: float
        Cosine spherical harmonic coefficients
    slm: float
        Sine spherical harmonic coefficients
    MJD: float
        output date as Modified Julian Day
    time: float
        output date in year-decimal

    References
    ----------
    .. [Loomis2020] B. D. Loomis, K. E. Rachlin, D. N. Wiese, F. W. Landerer,
        and S. B. Luthcke, "Replacing GRACE/GRACE‐FO *C*\ :sub:`30` with
        satellite laser ranging: Impacts on Antarctic Ice Sheet mass change".
        *Geophysical Research Letters*, 47, (2020).
        `doi: 10.1029/2019GL085488 <https://doi.org/10.1029/2019GL085488>`_
    """
    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')

    #-- read the file and get contents
    with open(os.path.expanduser(SLR_file), mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()
    file_lines = len(file_contents)

    #-- spherical harmonic degree range (5x5 with 6,1)
    LMIN = 2
    LMAX = 6
    n_harm = (LMAX**2 + 3*LMAX - LMIN**2 - LMIN)//2 - 5

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text
    while HEADER:
        #-- file line at count
        line = file_contents[count]
        #-- find the final line within the header text
        #-- to set HEADER flag to False when found
        HEADER = not bool(re.search(r'Product:',line))
        #-- add 1 to counter
        count += 1

    #-- number of dates within the file
    n_dates = (file_lines - count)//(n_harm + 1)
    #-- output spherical harmonic fields
    Ylms = {}
    Ylms['MJD'] = np.zeros((n_dates))
    Ylms['time'] = np.zeros((n_dates))
    Ylms['clm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    Ylms['slm'] = np.zeros((LMAX+1,LMAX+1,n_dates))
    #-- for each date
    for d in range(n_dates):
        #-- split the date line into individual components
        line_contents = file_contents[count].split()
        #-- modified Julian date of the beginning of the week
        Ylms['MJD'][d] = np.float64(line_contents[0])
        #-- date of the mid-point of the arc given in years
        Ylms['time'][d] = np.float64(line_contents[1])
        #-- add 1 to counter
        count += 1

        #-- read the spherical harmonic field
        for i in range(n_harm):
            #-- split the line into individual components
            line_contents = file_contents[count].split()
            #-- degree and order for the line
            l1 = np.int64(line_contents[0])
            m1 = np.int64(line_contents[1])
            #-- Spherical Harmonic data rescaled to output
            Ylms['clm'][l1,m1,d] = np.float64(line_contents[2])*SCALE
            Ylms['slm'][l1,m1,d] = np.float64(line_contents[3])*SCALE
            #-- add 1 to counter
            count += 1

    #-- return spherical harmonic fields and date information
    return Ylms

#-- PURPOSE: interpolate harmonics from 7-day to monthly
def convert_weekly(t_in, d_in, DATE=[], NEIGHBORS=28):
    """
    Interpolate harmonics from 7-day to 28-day

    Parameters
    ----------
    t_in: float
        Weekly time
    d_in: float
        Weekly harmonics
    DATE: list, default []
        Output monthly time for central averages
    NEIGHBORS: int, default 28
        Number of days to use in average

    Returns
    -------
    time: float
        output date in year-decimal
    month: int
        GRACE/GRACE-FO month
    data: float
        monthly spherical harmonic coefficients
    """
    #-- duplicate time and harmonics
    tdec = np.repeat(t_in, 7)
    data = np.repeat(d_in, 7)
    #-- calculate daily dates to use in centered moving average
    tdec += (np.mod(np.arange(len(tdec)),7) - 3.5)/365.25
    #-- calculate moving-average solution from 7-day arcs
    dinput = {}
    dinput['time'] = np.zeros_like(DATE)
    dinput['data'] = np.zeros_like(DATE,dtype='f8')
    #-- for each output monthly date
    for i,D in enumerate(DATE):
        #-- find all dates within NEIGHBORS days of mid-point
        isort = np.argsort((tdec - D)**2)[:NEIGHBORS]
        #-- calculate monthly mean of date and data
        dinput['time'][i] = np.mean(tdec[isort])
        dinput['data'][i] = np.mean(data[isort])
    #-- GRACE/GRACE-FO month
    dinput['month'] = gravity_toolkit.time.calendar_to_grace(dinput['time'])
    dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])
    #-- return the moving averages
    return dinput
