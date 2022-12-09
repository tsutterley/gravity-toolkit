#!/usr/bin/env python
u"""
hdf5_stokes.py
Written by Tyler Sutterley (11/2021)

Writes spherical harmonic coefficients to HDF5 files

CALLING SEQUENCE:
    hdf5_stokes(clm1,slm1,linp,minp,tinp,month,FILENAME=output_HDF5_file)

INPUTS:
    clm1: cosine spherical harmonic coefficients
    slm1: sine spherical harmonic coefficients
    linp: spherical harmonic degree (l)
    minp: spherical harmonic order (m)
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
    TITLE: description attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    DATE: harmonics have date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org)

UPDATE HISTORY:
    Updated 11/2021: use remapped dictionary for filling HDF5 variables
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 12/2020: added REFERENCE option to set file attribute
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
import logging
import numpy as np
import warnings

# attempt imports
try:
    import h5py
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

def hdf5_stokes(clm1, slm1, linp, minp, tinp, month, **kwargs):
    """
    Writes spherical harmonic coefficients to HDF5 files

    Parameters
    ----------
    clm1: cosine spherical harmonic coefficients
    slm1: sine spherical harmonic coefficients
    linp: spherical harmonic degree (l)
    minp: spherical harmonic order (m)
    tinp: date of measurement
    month: GRACE/GRACE-FO month

    FILENAME: HDF5 filename
    UNITS: spherical harmonic units
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    MONTHS_NAME: name of months variable within HDF5 file
    MONTHS_UNITS: months variable units
    MONTHS_LONGNAME: months variable description
    TITLE: description attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    DATE: harmonics have date information
    """
    # set default keyword arguments
    kwargs.setdefault('FILENAME',None)
    kwargs.setdefault('UNITS','Geodesy_Normalization')
    kwargs.setdefault('TIME_UNITS',None)
    kwargs.setdefault('TIME_LONGNAME',None)
    kwargs.setdefault('MONTHS_NAME','month')
    kwargs.setdefault('MONTHS_UNITS','number')
    kwargs.setdefault('MONTHS_LONGNAME','GRACE_month')
    kwargs.setdefault('TITLE',None)
    kwargs.setdefault('REFERENCE',None)
    kwargs.setdefault('DATE',True)
    kwargs.setdefault('CLOBBER',True)
    # set deprecation warning
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use harmonics.to_HDF5",
        DeprecationWarning)

    # setting HDF5 clobber attribute
    clobber = 'w' if kwargs['CLOBBER'] else 'w-'
    # opening HDF5 file for writing
    fileID = h5py.File(kwargs['FILENAME'], clobber)

    # Maximum spherical harmonic degree (LMAX) and order (MMAX)
    LMAX = np.max(linp)
    MMAX = np.max(minp)
    # Calculating the number of cos and sin harmonics up to LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm = (LMAX**2 + 3*LMAX - (LMAX-MMAX)**2 - (LMAX-MMAX))//2 + 1

    # dictionary with output variables
    output = {}
    # restructured degree and order
    output['l'] = np.zeros((n_harm,), dtype=np.int32)
    output['m'] = np.zeros((n_harm,), dtype=np.int32)
    # Restructuring output matrix to array format
    # will reduce matrix size and insure compatibility between platforms
    if kwargs['DATE']:
        n_time = len(np.atleast_1d(tinp))
        output['time'] = np.copy(tinp)
        output['month'] = np.copy(month)
        if (n_time == 1):
            output['clm'] = np.zeros((n_harm))
            output['slm'] = np.zeros((n_harm))
        else:
            output['clm'] = np.zeros((n_harm,n_time))
            output['slm'] = np.zeros((n_harm,n_time))
    else:
        n_time = 0
        output['clm'] = np.zeros((n_harm))
        output['slm'] = np.zeros((n_harm))

    # create counter variable lm
    lm = 0
    for m in range(0,MMAX+1):# MMAX+1 to include MMAX
        for l in range(m,LMAX+1):# LMAX+1 to include LMAX
            output['l'][lm] = np.int64(l)
            output['m'][lm] = np.int64(m)
            if (kwargs['DATE'] and (n_time > 1)):
                output['clm'][lm,:] = clm1[l,m,:]
                output['slm'][lm,:] = slm1[l,m,:]
            else:
                output['clm'][lm] = clm1[l,m]
                output['slm'][lm] = slm1[l,m]
            # add 1 to lm counter variable
            lm += 1

    # Defining the HDF5 dataset variables
    h5 = {}
    for key,val in output.items():
        h5[key] = fileID.create_dataset(key, val.shape,
            data=val, dtype=val.dtype, compression='gzip')

    # filling HDF5 dataset attributes
    # Defining attributes for degree and order
    h5['l'].attrs['long_name'] = 'spherical_harmonic_degree'# degree long name
    h5['l'].attrs['units'] = 'Wavenumber'# SH degree units
    h5['m'].attrs['long_name'] = 'spherical_harmonic_order'# order long name
    h5['m'].attrs['units'] = 'Wavenumber'# SH order units
    # Defining attributes for dataset
    h5['clm'].attrs['long_name'] = 'cosine_spherical_harmonics'
    h5['clm'].attrs['units'] = kwargs['UNITS']
    h5['slm'].attrs['long_name'] = 'sine_spherical_harmonics'
    h5['slm'].attrs['units'] = kwargs['UNITS']
    if kwargs['DATE']:
        # Defining attributes for date and month (or integer date)
        h5['time'].attrs['long_name'] = kwargs['TIME_LONGNAME']
        h5['time'].attrs['units'] = kwargs['TIME_UNITS']
        h5['month'].attrs['long_name'] = kwargs['MONTHS_LONGNAME']
        h5['month'].attrs['units'] = kwargs['MONTHS_UNITS']
    # description of file
    if kwargs['TITLE']:
        fileID.attrs['description'] = kwargs['TITLE']
    # reference of file
    if kwargs['REFERENCE']:
        fileID.attrs['reference'] = kwargs['REFERENCE']
    # date created
    fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())

    # Output HDF5 structure information
    logging.info(kwargs['FILENAME'])
    logging.info(list(fileID.keys()))

    # Closing the HDF5 file
    fileID.close()
