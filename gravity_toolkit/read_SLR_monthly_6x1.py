#!/usr/bin/env python
u"""
read_SLR_monthly_6x1.py
Written by Tyler Sutterley (05/2021)

Reads in monthly 5x5 spherical harmonic coefficients with 1
    coefficient from degree 6 all calculated from SLR measurements

Dataset distributed by UTCSR
    ftp://ftp.csr.utexas.edu/outgoing/cheng/slrgeo.5d561_187_naod
    http://download.csr.utexas.edu/pub/slr/degree_5/
        CSR_Monthly_5x5_Gravity_Harmonics.txt

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)

OUTPUTS:
    clm: Cosine spherical harmonic coefficients
    slm: Sine spherical harmonic coefficients
    error/clm: Cosine spherical harmonic coefficient uncertainty
    error/slm: Sine spherical harmonic coefficients uncertainty
    MJD: output date as Modified Julian Day
    time: output date in year-decimal

REFERENCES:
    Cheng, M., J. C.  Ries, and B. D. Tapley, 'Variations of the Earth's Figure
    Axis from Satellite Laser Ranging and GRACE', J. Geophys. Res., 116, B01409,
    2011, DOI:10.1029/2010JB000850.

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

#-- PURPOSE: read low degree harmonic data from Satellite Laser Ranging (SLR)
def read_SLR_monthly_6x1(input_file, HEADER=True):
    """
    Reads in monthly low degree and order spherical harmonic coefficients
    from Satellite Laser Ranging (SLR) measurements

    Arguments
    ---------
    input_file: input satellite laser ranging file from CSR

    Keyword arguments
    -----------------
    HEADER: file contains header text to be skipped

    Returns
    -------
    clm: Cosine spherical harmonic coefficients
    slm: Sine spherical harmonic coefficients
    error/clm: Cosine spherical harmonic coefficient uncertainty
    error/slm: Sine spherical harmonic coefficients uncertainty
    MJD: output date as Modified Julian Day
    time: output date in year-decimal
    """
    #-- check that SLR file exists
    if not os.access(os.path.expanduser(input_file), os.F_OK):
        raise FileNotFoundError('SLR file not found in file system')

    #-- read the file and get contents
    with open(os.path.expanduser(input_file),'r') as f:
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
            #-- fill anomaly field Ylms (variations and sigmas scaled by 1.0e10)
            Ylm_anomalies['clm'][l1,m1,d] = np.float64(line[2])*1e-10
            Ylm_anomalies['slm'][l1,m1,d] = np.float64(line[3])*1e-10
            Ylm_anomaly_error['clm'][l1,m1,d] = np.float64(line[6])*1e-10
            Ylm_anomaly_error['slm'][l1,m1,d] = np.float64(line[7])*1e-10
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
