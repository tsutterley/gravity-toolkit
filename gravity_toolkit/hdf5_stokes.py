#!/usr/bin/env python
u"""
hdf5_stokes.py
Written by Tyler Sutterley (03/2020)

Writes spherical harmonic coefficients to HDF5 files

CALLING SEQUENCE:
    hdf5_stokes(clm, slm, linp, minp, tinp, month, FILENAME=output_HDF5_file)

INPUTS:
    clm: Cosine Stokes Coefficient
    slm: Sine Stokes Coefficient
    linp: degree (l)
    minp: order (m)
    tinp: date of measurement
    month: GRACE/GRACE-FO month

OPTIONS:
    FILENAME: output filename HDF5
    UNITS: spherical harmonic units
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    MONTHS_NAME: name of months variable within HDF5 file
    MONTHS_UNITS: months variable units
    MONTHS_LONGNAME: months variable description
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    VERBOSE: will print to screen the HDF5 structure parameters
    DATE: harmonics have date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (http://h5py.org)

UPDATE HISTORY:
    Updated 03/2020: only include title if not None
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2019: don't include time (HH:MM:SS) in creation date
    Updated 07/2019: added creation date as a global attribute
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 12/2018: using python dictionaries to improve readability
    Updated 10/2018: using future division for python3 Compatibility
    Updated 02/2017: added MONTHS_UNITS, MONTHS_LONGNAME, MONTHS_NAME parameters
        aligned TIME_LONGNAME and TIME_UNITS with attributes
        can output a HDF5 file with multiple dates similar to the netcdf program
    Updated 06/2016: using __future__ print function
    Updated 03/2016: direct calculation of number of harmonics n_harm
    Updated 05/2015: minor change for MMAX != LMAX
    Updated 11/2014: got back to writing this
        in working condition with updated attributes as in netcdf equivalent
    Updated 12/2013: converted ncdf code to HDF5 code (alternative data type)
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 05/2013 made UNITS an option in case converting the units to
        mass harmonics or other harmonic variant
    Updated 03/2013: added units to clm and slm as 'Geodesy Normalization'
        switched I/O to column arrays for smaller file sizes and compatibility
            between languages
        made date an option for datasets that have no date
    Updated 01/2013 to add time and GRACE/GRACE-FO month number
    Written 07/2012
"""
from __future__ import print_function, division

import time
import h5py
import numpy as np

def hdf5_stokes(clm1, slm1, linp, minp, tinp, month, FILENAME=None,
    UNITS='Geodesy_Normalization', TIME_UNITS=None, TIME_LONGNAME=None,
    MONTHS_NAME='month', MONTHS_UNITS='number', MONTHS_LONGNAME='GRACE_month',
    TITLE=None, DATE=True, CLOBBER=True, VERBOSE=False):

    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- opening HDF5 file for writing
    fileID = h5py.File(FILENAME, clobber)

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

    #-- Defining the HDF5 dataset variables
    h5 = {}
    h5['l'] = fileID.create_dataset('l', (n_harm,), \
        data=lout, dtype=np.int, compression='gzip')
    h5['m'] = fileID.create_dataset('m', (n_harm,), \
        data=mout, dtype=np.int, compression='gzip')
    if DATE:
        h5['time'] = fileID.create_dataset('time', (n_time,), \
            data=tinp, dtype=np.float, compression='gzip')
        h5['month'] = fileID.create_dataset(MONTHS_NAME, (n_time,), \
            data=month, dtype=np.int, compression='gzip')
    #-- if more than 1 date in file
    if (n_time > 1):
        h5['clm'] = fileID.create_dataset('clm', (n_harm,n_time,), \
            data=clm, dtype=np.float, compression='gzip')
        h5['slm'] = fileID.create_dataset('slm', (n_harm,n_time,), \
            data=slm, dtype=np.float, compression='gzip')
    else:
        h5['clm'] = fileID.create_dataset('clm', (n_harm,), \
            data=clm, dtype=np.float, compression='gzip')
        h5['slm'] = fileID.create_dataset('slm', (n_harm,), \
            data=slm, dtype=np.float, compression='gzip')

    #-- filling HDF5 dataset attributes
    #-- Defining attributes for degree and order
    h5['l'].attrs['long_name'] = 'spherical_harmonic_degree'#-- degree long name
    h5['l'].attrs['units'] = 'Wavenumber'#-- SH degree units
    h5['m'].attrs['long_name'] = 'spherical_harmonic_order'#-- order long name
    h5['m'].attrs['units'] = 'Wavenumber'#-- SH order units
    #-- Defining attributes for dataset
    h5['clm'].attrs['long_name'] = 'cosine_spherical_harmonics'
    h5['clm'].attrs['units'] = UNITS
    h5['slm'].attrs['long_name'] = 'sine_spherical_harmonics'
    h5['slm'].attrs['units'] = UNITS
    if DATE:
        #-- Defining attributes for date and month (or integer date)
        h5['time'].attrs['long_name'] = TIME_LONGNAME
        h5['time'].attrs['units'] = TIME_UNITS
        h5['month'].attrs['long_name'] = MONTHS_LONGNAME
        h5['month'].attrs['units'] = MONTHS_UNITS
    #-- description of file
    if TITLE:
        fileID.attrs['description'] = TITLE
    #-- date created
    fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output HDF5 structure information
    if VERBOSE:
        print(FILENAME)
        print(list(fileID.keys()))

    #-- Closing the HDF5 file
    fileID.close()
