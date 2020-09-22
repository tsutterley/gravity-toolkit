#!/usr/bin/env python
u"""
ncdf_stokes.py
Written by Tyler Sutterley (07/2020)

Writes spherical harmonic coefficients to netCDF4 files

CALLING SEQUENCE:
    ncdf_stokes(clm1,slm1,linp,minp,tinp,month,FILENAME=output_netcdf4_file)

INPUTS:
    clm1: cosine spherical harmonic coefficients
    slm1: sine spherical harmonic coefficients
    linp: spherical harmonic degree (l)
    minp: spherical harmonic order (m)
    tinp: date of measurement
    month: GRACE/GRACE-FO month

OPTIONS:
    FILENAME: output netCDF4 filename
    UNITS: spherical harmonic units
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    MONTHS_NAME: name of months variable within netCDF4 file
    MONTHS_UNITS: months variable units
    MONTHS_LONGNAME: months variable description
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing netCDF4 file
    VERBOSE: will print to screen the netCDF4 structure parameters
    DATE: harmonics have date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 03/2020: only include title if not None
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2019: don't include time (HH:MM:SS) in creation date
    Updated 07/2019: added creation date as a global attribute
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 12/2018: using python dictionaries to improve readability
    Updated 10/2018: using future division for python3 Compatibility
    Updated 02/2017: added MONTHS_UNITS, MONTHS_LONGNAME, MONTHS_NAME parameters
        aligned TIME_LONGNAME and TIME_UNITS with attributes
    Updated 07/2016: using netCDF4-python
    Updated 06/2016: using __future__ print function
    Updated 03/2016: direct calculation of number of harmonics n_harm
    Updated 07/2015: forgot to add case for DATE=False
    Updated 06/2015: can output single netcdf with multiple dates
    Updated 05/2015: minor change for MMAX != LMAX
    Updated 05/2014: new parameters for time attributes
    Updated 02/2014: minor update to if statements
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 05/2013 made UNITS an option in case converting the units to
        mass harmonics or other harmonic variant
    Updated 03/2013: added units to clm and slm as 'Geodesy Normalization'
        switched I/O to column arrays for smaller file sizes and compatibility
            between languages
        made date an option for datasets that have no date (e.g. GIA)
    Updated 01/2013 to add time and GRACE/GRACE-FO month number
    Written 07/2012
"""
from __future__ import print_function, division

import time
import netCDF4
import numpy as np

def ncdf_stokes(clm1, slm1, linp, minp, tinp, month, FILENAME=None,
    UNITS='Geodesy_Normalization', TIME_UNITS=None, TIME_LONGNAME=None,
    MONTHS_NAME='month', MONTHS_UNITS='number', MONTHS_LONGNAME='GRACE_month',
    TITLE=None, DATE=True, CLOBBER=True, VERBOSE=False):
    """
    Writes spherical harmonic coefficients to netCDF4 files

    Arguments
    ---------
    clm1: cosine spherical harmonic coefficients
    slm1: sine spherical harmonic coefficients
    linp: spherical harmonic degree (l)
    minp: spherical harmonic order (m)
    tinp: date of measurement
    month: GRACE/GRACE-FO month

    Keyword arguments
    -----------------
    FILENAME: netCDF4 filename
    UNITS: spherical harmonic units
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    MONTHS_NAME: name of months variable within netCDF4 file
    MONTHS_UNITS: months variable units
    MONTHS_LONGNAME: months variable description
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing netCDF4 file
    VERBOSE: will print to screen the netCDF4 structure parameters
    DATE: harmonics have date information
    """

    #-- setting NetCDF clobber attribute
    clobber = 'w' if CLOBBER else 'a'
    #-- opening netCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, clobber, format="NETCDF4")

    #-- Maximum spherical harmonic degree (LMAX) and order (MMAX)
    LMAX = np.max(linp)
    MMAX = np.max(minp)
    #-- Calculating the number of cos and sin harmonics up to LMAX
    #-- taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm = (LMAX**2 + 3*LMAX - (LMAX-MMAX)**2 - (LMAX-MMAX))//2 + 1

    #-- Restructuring output matrix to array format
    #-- will reduce matrix size and insure compatibility between platforms
    if DATE:
        if (np.ndim(tinp) == 0):
            n_time = 1
            clm = np.zeros((n_harm))
            slm = np.zeros((n_harm))
        else:
            n_time = len(tinp)
            clm = np.zeros((n_harm,n_time))
            slm = np.zeros((n_harm,n_time))
    else:
        n_time = 0
        clm = np.zeros((n_harm))
        slm = np.zeros((n_harm))

    #-- restructured degree and order
    lout = np.zeros((n_harm,), dtype=np.int32)
    mout = np.zeros((n_harm,), dtype=np.int32)
    #-- create counter variable lm
    lm = 0
    for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
        for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
            lout[lm] = np.int(l)
            mout[lm] = np.int(m)
            if (DATE and (n_time > 1)):
                clm[lm,:] = clm1[l,m,:]
                slm[lm,:] = slm1[l,m,:]
            else:
                clm[lm] = clm1[l,m]
                slm[lm] = slm1[l,m]
            #-- add 1 to lm counter variable
            lm += 1

    #-- Defining the netCDF dimensions
    fileID.createDimension('lm', n_harm)
    if DATE:
        fileID.createDimension('time', n_time)

    #-- defining the netCDF variables
    nc = {}
    #-- degree and order
    nc['l'] = fileID.createVariable('l', 'i', ('lm',))
    nc['m'] = fileID.createVariable('m', 'i', ('lm',))
    #-- spherical harmonics
    if (DATE and (n_time > 1)):
        nc['clm'] = fileID.createVariable('clm', 'd', ('lm','time',))
        nc['slm'] = fileID.createVariable('slm', 'd', ('lm','time',))
    else:
        nc['clm'] = fileID.createVariable('clm', 'd', ('lm',))
        nc['slm'] = fileID.createVariable('slm', 'd', ('lm',))
    if DATE:
        #-- time (in decimal form)
        nc['time'] = fileID.createVariable('time', 'd', ('time',))
        #-- GRACE/GRACE-FO month (or integer date)
        nc['month'] = fileID.createVariable(MONTHS_NAME, 'i', ('time',))

    #-- filling netCDF variables
    nc['l'][:] = lout.copy()
    nc['m'][:] = mout.copy()
    nc['clm'][:] = clm.copy()
    nc['slm'][:] = slm.copy()
    if DATE:
        nc['time'][:] = tinp
        nc['month'][:] = month

    #-- Defining attributes for degree and order
    nc['l'].long_name = 'spherical_harmonic_degree'#-- SH degree long name
    nc['l'].units = 'Wavenumber'#-- SH degree units
    nc['m'].long_name = 'spherical_harmonic_order'#-- SH order long name
    nc['m'].units = 'Wavenumber'#-- SH order units
    #-- Defining attributes for harmonics
    nc['clm'].long_name = 'cosine_spherical_harmonics'
    nc['clm'].units = UNITS
    nc['slm'].long_name = 'sine_spherical_harmonics'
    nc['slm'].units = UNITS
    if DATE:
        #-- Defining attributes for date and month
        nc['time'].long_name = TIME_LONGNAME
        nc['time'].units = TIME_UNITS
        nc['month'].long_name = MONTHS_LONGNAME
        nc['month'].units = MONTHS_UNITS
    #-- global variable of netCDF file
    if TITLE:
        fileID.TITLE = TITLE
    #-- date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output netCDF structure information
    if VERBOSE:
        print(FILENAME)
        print(list(fileID.variables.keys()))

    #-- Closing the netCDF file
    fileID.close()
