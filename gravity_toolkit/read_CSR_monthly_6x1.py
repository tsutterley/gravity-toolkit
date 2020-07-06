#!/usr/bin/env python
u"""
read_CSR_monthly_6x1.py
Written by Tyler Sutterley (07/2020)

Reads in monthly 5x5 spherical harmonic coefficients with 1
    coefficient from degree 6 all calculated from SLR measurements

Dataset distributed by UTCSR
    ftp://ftp.csr.utexas.edu/outgoing/cheng/slrgeo.5d561_187_naod

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)

OUTPUTS:
    clm: Cosine spherical harmonic coefficients
    slm: Sine spherical harmonic coefficients
    error/clm: Cosine spherical harmonic coefficient uncertainty
    error/slm: Sine spherical harmonic coefficients uncertainty
    MJD: output date as Modified Julian Day
    time: output date in year-decimal

REFERENCE:
    Cheng, M., J. C.  Ries, and B. D. Tapley, 'Variations of the Earth's Figure
    Axis from Satellite Laser Ranging and GRACE', J. Geophys. Res., 116, B01409,
    2011, DOI:10.1029/2010JB000850.

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    convert_calendar_decimal.py: converts from calendar dates to decimal years

UPDATE HISTORY:
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
from gravity_toolkit.convert_calendar_decimal import convert_calendar_decimal

#-- PURPOSE: read low degree harmonic data from Satellite Laser Ranging (SLR)
def read_CSR_monthly_6x1(input_file, HEADER=True):
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

    #-- read the file and get contents
    with open(os.path.expanduser(input_file),'r') as f:
        file_contents = f.read().splitlines()
    file_lines = len(file_contents)

    #-- spherical harmonic degree range (full 5x5 with 6,1)
    LMIN = 1
    LMAX = 6
    n_harm = (LMAX**2 + 3*LMAX - LMIN**2 - LMIN)//2 - 5

    #-- counts the number of lines in the header
    count = 0
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

    #-- number of dates within the file
    n_dates = (file_lines - count)//(n_harm + 1)

    #-- read mean fields from the header
    mean_Ylms = {}
    mean_Ylm_error = {}
    mean_Ylms['clm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylms['slm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylm_error['clm'] = np.zeros((LMAX+1,LMAX+1))
    mean_Ylm_error['slm'] = np.zeros((LMAX+1,LMAX+1))
    for i in range(n_harm+1):
        #-- split the line into individual components
        line = file_contents[indice+i].split()
        #-- degree and order for the line
        l1 = np.int(line[0])
        m1 = np.int(line[1])
        #-- fill mean field Ylms
        mean_Ylms['clm'][l1,m1] = np.float(line[2].replace('D','E'))
        mean_Ylms['slm'][l1,m1] = np.float(line[3].replace('D','E'))
        mean_Ylm_error['clm'][l1,m1] = np.float(line[4].replace('D','E'))
        mean_Ylm_error['slm'][l1,m1] = np.float(line[5].replace('D','E'))

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
        #-- modified Julian date of the middle of the month
        Ylms['MJD'][d] = np.mean(np.array(line_contents[5:7],dtype=np.float))
        #-- date of the mid-point of the arc given in years
        YY,MM = np.array(line_contents[3:5])
        Ylms['time'][d] = convert_calendar_decimal(YY,MM)
        #-- add 1 to counter
        count += 1

        #-- read the anomaly field
        for i in range(n_harm):
            #-- split the line into individual components
            line = file_contents[count].split()
            #-- degree and order for the line
            l1 = np.int(line[0])
            m1 = np.int(line[1])
            #-- fill anomaly field Ylms (variations and sigmas scaled by 1.0e10)
            Ylm_anomalies['clm'][l1,m1,d] = np.float(line[2])*1e-10
            Ylm_anomalies['slm'][l1,m1,d] = np.float(line[3])*1e-10
            Ylm_anomaly_error['clm'][l1,m1,d] = np.float(line[6])*1e-10
            Ylm_anomaly_error['slm'][l1,m1,d] = np.float(line[7])*1e-10
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
