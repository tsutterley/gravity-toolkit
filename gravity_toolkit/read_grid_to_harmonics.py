#!/usr/bin/env python
u"""
read_grid_to_harmonics.py
Written by Hugo Lecomte (12/2020)

Reads netCDF file with grid data and extracts spherical harmonic from those data
Correct data for drift in pole tide following Wahr et al. (2015)
Parses date of GRACE/GRACE-FO data from filename

Design for JPL MASCON netCDF data available on
https://podaac-tools.jpl.nasa.gov/drive/files
In the folder /allData/tellus/retired/L3/mascon/RL06/JPL/v02

INPUTS:
    input_file: GRACE/GRACE-FO Level-3 netCDF grid data file
    LMAX: Maximum degree of spherical harmonics (degree of truncation)

OPTIONS:
    MMAX: Maximum order of spherical harmonics (order of truncation)
        default is the maximum spherical harmonic degree
    POLE_TIDE: correct GSM data for pole tide drift following Wahr et al. (2015)

OUTPUTS:
    time: mid-month date in year-decimal
    start: start date of range as Julian day
    end: end date of range as Julian day
    clm: cosine spherical harmonics of input data (LMAX,MMAX)
    slm: sine spherical harmonics of input data (LMAX,MMAX)
    eclm: cosine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)
    eslm: sine spherical harmonic uncalibrated standard deviations (LMAX,MMAX)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Written 12/2020
"""
import os
import re
import io
import numpy as np
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.hdf5_read import hdf5_read
from gravity_toolkit.utilities import get_data_path
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.gen_stokes import gen_stokes

#-- PURPOSE: read Level-3 GRACE and GRACE-FO netCDF files
def read_grid_to_harmonics(input_file, VARNAME, LMAX, MMAX=None, LONNAME='lon',
                         LATNAME='lat', TIMENAME='time', UNITS=1, POLE_TIDE=False):
    """
    Reads netCDF or HDF5 file with grid data and extracts spherical harmonic from those data
    Correct data prior to Release 6 for pole tide drift
    Parses date of GRACE/GRACE-FO data from filename

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-3 netCDF grid data file
    VARNAME: z variable name in the file
    LMAX: Maximum degree of spherical harmonics (degree of truncation)

    Keyword arguments
    -----------------
    MMAX: Maximum order of spherical harmonics
    LONNAME: longitude variable name in the file
    LATNAME: latitude variable name in the file
    TIMENAME: time variable name in the file
    UNITS: input data units
        1: cm of water thickness
        2: Gtons of mass
        3: kg/m^2
    POLE_TIDE: correct for pole tide drift following Wahr et al. (2015)

    Returns
    -------
    clm: GRACE/GRACE-FO cosine spherical harmonics
    slm: GRACE/GRACE-FO sine spherical harmonics
    time: time of each GRACE/GRACE-FO measurement (mid-month)
    month: GRACE/GRACE-FO months of input datasets
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    title: string denoting low degree zonals replacement, geocenter usage and corrections
    directory: directory of exact GRACE/GRACE-FO product
    """

    #-- parse filename
    pfx,center,time,realm,release,v_id,sfx = parse_file(input_file)

    #-- read file content
    if input_file[-3:] == '.nc':
        file_contents = ncdf_read(input_file, DATE=True, VARNAME=VARNAME, LONNAME=LONNAME,
                                  LATNAME=LATNAME, TIMENAME=TIMENAME, ATTRIBUTES=True,
                                  TITLE=True, COMPRESSION=sfx)
    elif input_file[-4:] == '.hdf' or input_file[-3:] == '.h5' or input_file[-5:] == '.hdf5':
        file_contents = hdf5_read(input_file, DATE=True, VARNAME=VARNAME, LONNAME=LONNAME,
                                  LATNAME=LATNAME, TIMENAME=TIMENAME, ATTRIBUTES=True,
                                  TITLE=True, COMPRESSION=sfx)

    #-- load love numbers
    hl, kl, ll = read_love_numbers(get_data_path(['data', 'love_numbers']), REFERENCE='CF')

    #-- set maximum spherical harmonic order
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX

    #-- number of dates in data
    n_time = file_contents['time'].shape[0]
    #-- Spherical harmonic coefficient matrices to be filled from data file
    grace_clm = np.zeros((LMAX + 1, MMAX + 1, n_time))
    grace_slm = np.zeros((LMAX + 1, MMAX + 1, n_time))
    #-- Time matrix to fill
    tdec = np.zeros((n_time))
    month = np.zeros((n_time))
    #-- output dimensions
    lout = np.arange(LMAX + 1)
    mout = np.arange(MMAX + 1)

    #-- for each date, conversion to spherical harmonics
    for i in range(n_time):
        harmo = gen_stokes(file_contents['data'][i, :, :],
                   file_contents['lon'][:], file_contents['lat'][:],
                           LMAX=LMAX, MMAX=MMAX, UNITS=UNITS, LOVE=(hl, kl, ll))

        grace_clm[:, :, i] = harmo['clm']
        grace_slm[:, :, i] = harmo['slm']

    #-- extract GRACE date information from input file name
    start_yr = np.float(time[:4])

    #-- variables initialization for date conversion
    current_year = start_yr
    current_month = 1
    cmp_past_dpm = 0
    cmp_past_dpy = 0
    if (start_yr % 4) == 0:#-- Leap Year (% = modulus)
        dpy = 366.0
        dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:#-- Standard Year
        dpy = 365.0
        dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    #-- for each date, conversion to month and decimal year
    for i in range(n_time):
        #-- Month iteration
        while file_contents['time'][i] - cmp_past_dpm > dpm[(current_month - 1)%12]:
            current_month += 1
            cmp_past_dpm += dpm[(current_month - 1)%12]

        #-- Year iteration
        while file_contents['time'][i] - cmp_past_dpy > dpy:
            current_year += 1
            cmp_past_dpy += dpy
            if (current_year % 4) == 0:  #-- Leap Year (% = modulus)
                dpy = 366.0
                dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            else:  #-- Standard Year
                dpy = 365.0
                dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        tdec[i] = current_year + (file_contents['time'][i] - cmp_past_dpy)/dpy
        month[i] = current_month

        #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
        #-- Accelerometer shutoffs make this relation between month number
        #-- and date more complicated as days from other months are used
        #-- May15 (month 161) is centered in Apr15 (160)
        if (month[i] == 160) and (month[i] == month[i - 1]):
            month[i] = month[i - 1] + 1

    #-- extract GRACE and GRACE-FO file informations
    title = file_contents['attributes']

    #-- Correct Pole Tide following Wahr et al. (2015) 10.1002/2015JB011986
    if POLE_TIDE:
        for i in range(n_time):
            #-- time since 2000.0
            dt = tdec[i] - 2000.0

            #-- JPL Pole Tide Correction
            #-- values for IERS mean pole [2010]
            if tdec[i] < 2010.0:
                a = np.array([0.055974,1.8243e-3,1.8413e-4,7.024e-6])
                b = np.array([-0.346346,-1.7896e-3,1.0729e-4,0.908e-6])
            elif tdec[i] >= 2010.0:
                a = np.array([0.023513,7.6141e-3,0.0,0.0])
                b = np.array([-0.358891,0.6287e-3,0.0,0.0])
            #-- calculate m1 and m2 values
            m1 = np.copy(a[0])
            m2 = np.copy(b[0])
            for x in range(1,4):
                m1 += a[x]*dt**x
                m2 += b[x]*dt**x
            #-- pole tide values for JPL
            #-- JPL remove the IERS mean pole from m1 and m2
            #-- before computing their harmonic solutions
            C21_PT = -1.551e-9*(m1 - 0.62e-3*dt) - 0.012e-9*(m2 + 3.48e-3*dt)
            S21_PT = 0.021e-9*(m1 - 0.62e-3*dt) - 1.505e-9*(m2 + 3.48e-3*dt)
            #-- correct GRACE spherical harmonics for pole tide
            #-- note: -= means grace_xlm = grace_xlm - PT
            grace_clm[2, 1, i] -= C21_PT
            grace_clm[2, 1, i] -= S21_PT

    #-- return the GRACE data, GRACE date (mid-month in decimal), and the
    #-- start and end days as Julian dates
    return {'clm': grace_clm, 'slm': grace_slm, 'time': tdec, 'month': month,
            'l': lout, 'm': mout, 'title': title, 'directory': os.path.split(input_file)[0]}

#-- PURPOSE: extract parameters from filename
def parse_file(input_file):
    """
    Extract parameters from filename

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    """
    #-- compile numerical expression operator for parameters from files
    #-- JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    regex_pattern = r'(.*?)\.(.*?)\.(.*?)\.(.*?)\.(.*?)\.(.*?)\.(\w{2,})'
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- extract parameters from input filename
    if isinstance(input_file, io.IOBase):
        return rx.findall(input_file.filename).pop()
    else:
        return rx.findall(os.path.basename(input_file)).pop()
