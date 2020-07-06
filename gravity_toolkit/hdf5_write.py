#!/usr/bin/env python
u"""
hdf5_write.py
Written by Tyler Sutterley (07/2020)

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
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    VERBOSE: will print to screen the HDF5 structure parameters
    DATE: data has date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org)

UPDATE HISTORY:
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
import numpy as np

def hdf5_write(data, lon, lat, tim, FILENAME=None, VARNAME='z', LONNAME='lon',
    LATNAME='lat', TIMENAME='time', UNITS=None, LONGNAME=None, FILL_VALUE=None,
    TIME_UNITS=None, TIME_LONGNAME=None, TITLE=None, DATE=True, CLOBBER=True,
    VERBOSE=False):
    """
    Writes spatial data to HDF5 files

    Arguments
    ---------
    data: z data
    lon: longitude array
    lat: latitude array
    tim: time array

    Keyword arguments
    -----------------
    FILENAME: HDF5 filename
    VARNAME: z variable name in HDF5 file
    LONNAME: longitude variable name in HDF5 file
    LATNAME: latitude variable name in HDF5 file
    UNITS: z variable units
    LONGNAME: z variable description
    FILL_VALUE: missing value for z variable
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    TITLE: title attribute of dataset
    CLOBBER: will overwrite an existing HDF5 file
    VERBOSE: will print to screen the HDF5 structure parameters
    DATE: data has date information
    """

    #-- setting HDF5 clobber attribute
    if CLOBBER in ('Y','y'):
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- Dimensions of time parameters
    n_time = 1 if (np.ndim(tim) == 0) else len(tim)

    #-- opening HDF5 file for writing
    fileID = h5py.File(FILENAME, clobber)
    #-- Defining the HDF5 dataset variables
    h5 = {}
    h5[LONNAME] = fileID.create_dataset(LONNAME, lon.shape, data=lon,
        dtype=lon.dtype, compression='gzip')
    h5[LATNAME] = fileID.create_dataset(LATNAME, lat.shape, data=lat,
        dtype=lat.dtype, compression='gzip')
    h5[VARNAME] = fileID.create_dataset(VARNAME, data.shape, data=data,
        dtype=data.dtype, fillvalue=FILL_VALUE, compression='gzip')
    if DATE:
        h5[TIMENAME] = fileID.create_dataset(TIMENAME, (n_time,), data=tim,
            dtype=np.float, compression='gzip')
    #-- add dimensions
    h5[VARNAME].dims[0].label=LATNAME
    h5[VARNAME].dims[0].attach_scale(h5[LATNAME])
    h5[VARNAME].dims[1].label=LONNAME
    #-- if more than 1 date in file
    if (n_time > 1):
        h5[VARNAME].dims[2].label=TIMENAME
        h5[VARNAME].dims[2].attach_scale(h5[TIMENAME])

    #-- filling HDF5 dataset attributes
    #-- Defining attributes for longitude and latitude
    h5[LONNAME].attrs['long_name'] = 'longitude'
    h5[LONNAME].attrs['units'] = 'degrees_east'
    h5[LATNAME].attrs['long_name'] = 'latitude'
    h5[LATNAME].attrs['units'] = 'degrees_north'
    #-- Defining attributes for dataset
    h5[VARNAME].attrs['long_name'] = LONGNAME
    h5[VARNAME].attrs['units'] = UNITS
    #-- Dataset contains missing values
    if (FILL_VALUE is not None):
        h5[VARNAME].attrs['_FillValue'] = FILL_VALUE
    #-- Defining attributes for date
    if DATE:
        h5[TIMENAME].attrs['long_name'] = TIME_LONGNAME
        h5[TIMENAME].attrs['units'] = TIME_UNITS
    #-- description of file
    if TITLE:
        fileID.attrs['description'] = TITLE
    #-- date created
    fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output HDF5 structure information
    if VERBOSE in ('Y','y'):
        print(FILENAME)
        print(list(fileID.keys()))

    #-- Closing the HDF5 file
    fileID.close()
