#!/usr/bin/env python
u"""
hdf5_write.py
Written by Tyler Sutterley (11/2021)

Writes spatial data to HDF5 files

CALLING SEQUENCE:
    hdf5_write(data, lon, lat, tim, FILENAME=output_HDF5_file)

INPUTS:
    data: z data
    lon: longitude array
    lat: latitude array
    tim: time array

OPTIONS:
    FILENAME: output filename HDF5
    VARNAME: z variable name in HDF5 file
    LONNAME: longitude variable name in HDF5 file
    LATNAME: latitude variable name in HDF5 file
    UNITS: z variable units
    LONGNAME: z variable description
    FILL_VALUE: missing value for z variable
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    TITLE: description attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    DATE: data has date information

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
    Updated 04/2020: added option DATE if including time data
    Updated 03/2020: only include title if not None
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019 for public release
    Updated 08/2019: don't include time (HH:MM:SS) in creation date
    Updated 07/2019: added creation date as a global attribute
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 08/2018: use n_time variable for output HDF5 dimensions
    Updated 03/2018: added option TIMENAME to specify the variable name of time
    Updated 02/2017: TIME_LONGNAME and TIME_UNITS with attributes,
        updated TIME_LONGNAME to Date_in_Decimal_Years
    Updated 06/2016: using __future__ print function
    Updated 05/2016: output data types same as input data types
    Updated 11/2014: got back to writing this
        in working condition with updated attributes as in netcdf equivalent
    Updated 12/2013: converted ncdf code to HDF5 code (alternative data type)
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 01/2013: adding time as a variable
    Updated 10/2012: changed from variable names x and y to lon and lat.
    Written 07/2012
"""
from __future__ import print_function

import time
import h5py
import logging
import numpy as np
import warnings

def hdf5_write(data, lon, lat, tim, **kwargs):
    """
    Writes spatial data to HDF5 files

    Parameters
    ----------
    data: z data
    lon: longitude array
    lat: latitude array
    tim: time array

    FILENAME: HDF5 filename
    VARNAME: z variable name in HDF5 file
    LONNAME: longitude variable name in HDF5 file
    LATNAME: latitude variable name in HDF5 file
    UNITS: z variable units
    LONGNAME: z variable description
    FILL_VALUE: missing value for z variable
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    TITLE: description attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    DATE: data has date information
    """
    #-- set default keyword arguments
    kwargs.setdefault('FILENAME',None)
    kwargs.setdefault('VARNAME','z')
    kwargs.setdefault('LONNAME','lon')
    kwargs.setdefault('LATNAME','lat')
    kwargs.setdefault('TIMENAME','time')
    kwargs.setdefault('UNITS',None)
    kwargs.setdefault('LONGNAME',None)
    kwargs.setdefault('FILL_VALUE',None)
    kwargs.setdefault('TIME_UNITS',None)
    kwargs.setdefault('TIME_LONGNAME',None)
    kwargs.setdefault('TITLE',None)
    kwargs.setdefault('REFERENCE',None)
    kwargs.setdefault('DATE',True)
    kwargs.setdefault('CLOBBER',True)
    #-- set deprecation warning
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use spatial.to_HDF5",
        DeprecationWarning)

    #-- setting HDF5 clobber attribute
    clobber = 'w' if kwargs['CLOBBER'] else 'w-'

    #-- opening HDF5 file for writing
    fileID = h5py.File(kwargs['FILENAME'], clobber)

    #-- create output dictionary with key mapping
    output = {}
    output[kwargs['LONNAME']] = np.copy(lon)
    output[kwargs['LATNAME']] = np.copy(lat)
    dimensions = [kwargs['LATNAME'],kwargs['LONNAME']]
    #-- extend with date variables
    if kwargs['DATE']:
        output[kwargs['TIMENAME']] = np.array(tim,dtype='f')
        output[kwargs['VARNAME']] = np.atleast_3d(data)
        dimensions.append(kwargs['TIMENAME'])
    else:
        output[kwargs['VARNAME']] = np.copy(data)

    #-- Defining the HDF5 dataset variables
    h5 = {}
    for key,val in output.items():
        h5[key] = fileID.create_dataset(key, val.shape, data=val,
            dtype=val.dtype, compression='gzip')
    #-- add dimensions
    for i,dim in enumerate(dimensions):
        h5[kwargs['VARNAME']].dims[i].label = dim
        h5[kwargs['VARNAME']].dims[i].attach_scale(h5[dim])

    #-- filling HDF5 dataset attributes
    #-- Defining attributes for longitude and latitude
    h5[kwargs['LONNAME']].attrs['long_name'] = 'longitude'
    h5[kwargs['LONNAME']].attrs['units'] = 'degrees_east'
    h5[kwargs['LATNAME']].attrs['long_name'] = 'latitude'
    h5[kwargs['LATNAME']].attrs['units'] = 'degrees_north'
    #-- Defining attributes for dataset
    h5[kwargs['VARNAME']].attrs['long_name'] = kwargs['LONGNAME']
    h5[kwargs['VARNAME']].attrs['units'] = kwargs['UNITS']
    #-- Dataset contains missing values
    if (kwargs['FILL_VALUE'] is not None):
        h5[kwargs['VARNAME']].attrs['_FillValue'] = kwargs['FILL_VALUE']
    #-- Defining attributes for date
    if kwargs['DATE']:
        h5[kwargs['TIMENAME']].attrs['long_name'] = kwargs['TIME_LONGNAME']
        h5[kwargs['TIMENAME']].attrs['units'] = kwargs['TIME_UNITS']
    #-- description of file
    if kwargs['TITLE']:
        fileID.attrs['description'] = kwargs['TITLE']
    #-- reference of file
    if kwargs['REFERENCE']:
        fileID.attrs['reference'] = kwargs['REFERENCE']
    #-- date created
    fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output HDF5 structure information
    logging.info(kwargs['FILENAME'])
    logging.info(list(fileID.keys()))

    #-- Closing the HDF5 file
    fileID.close()
