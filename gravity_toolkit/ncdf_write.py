#!/usr/bin/env python
u"""
ncdf_write.py
Written by Tyler Sutterley (11/2021)

Writes spatial data to COARDS-compliant netCDF4 files

CALLING SEQUENCE:
    ncdf_write(data, lon, lat, tim, FILENAME=output_netcdf4_file)

INPUTS:
    data: z data
    lon: longitude array
    lat: latitude array
    tim: time array

OPTIONS:
    FILENAME: output netCDF4 filename
    VARNAME: z variable name in netCDF4 file
    LONNAME: longitude variable name in netCDF4 file
    LATNAME: latitude variable name in netCDF4 file
    UNITS: z variable units
    LONGNAME: z variable description
    FILL_VALUE: missing value for z variable
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    TITLE: title attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing netCDF4 file
    DATE: data has date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
    Updated 11/2021: use remapped dictionary for filling netCDF4 variables
    Updated 10/2021: using python logging for handling verbose output
    Updated 12/2020: added REFERENCE option to set file attribute
    Updated 07/2020: added function docstrings
    Updated 04/2020: added option DATE if including time data
    Updated 03/2020: only include title if not None
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019 for public release
    Updated 08/2019: don't include time (HH:MM:SS) in creation date
    Updated 07/2019: added creation date as a global attribute
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 03/2018: added option TIMENAME to specify the variable name of time
    Updated 02/2017: TIME_LONGNAME and TIME_UNITS with attributes,
        updated TIME_LONGNAME to Date_in_Decimal_Years
    Updated 07/2016: using netCDF4-python with zlib compression
    Updated 06/2016: using __future__ print function
    Updated 05/2016: output data types same as input data types
    Updated 11/2014: new parameters for variable names and attributes
    Updated 05/2014: new parameters for time attributes, and missing values
    Updated 02/2014: minor update to if statements
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 01/2013: adding time as a variable
    Updated 10/2012: changed from variable names x and y to lon and lat.
    Written 07/2012
"""
from __future__ import print_function

import time
import logging
import netCDF4
import numpy as np
import warnings

def ncdf_write(data, lon, lat, tim, **kwargs):
    """
    Writes spatial data to COARDS-compliant netCDF4 files

    Parameters
    ----------
    data: z data
    lon: longitude array
    lat: latitude array
    tim: time array

    FILENAME: netCDF4 filename
    VARNAME: z variable name in netCDF4 file
    LONNAME: longitude variable name in netCDF4 file
    LATNAME: latitude variable name in netCDF4 file
    UNITS: z variable units
    LONGNAME: z variable description
    FILL_VALUE: missing value for z variable
    TIME_UNITS: time variable units
    TIME_LONGNAME: time variable description
    TITLE: title attribute of dataset
    REFERENCE: reference attribute of dataset
    CLOBBER: will overwrite an existing netCDF4 file
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
    warnings.warn("Deprecated. Please use spatial.to_netCDF4",
        DeprecationWarning)

    #-- setting NetCDF clobber attribute
    clobber = 'w' if kwargs['CLOBBER'] else 'a'
    #-- opening NetCDF file for writing
    #-- Create the NetCDF file
    fileID = netCDF4.Dataset(kwargs['FILENAME'], clobber, format="NETCDF4")

    #-- create output dictionary with key mapping
    output = {}
    output[kwargs['LONNAME']] = np.copy(lon)
    output[kwargs['LATNAME']] = np.copy(lat)
    dimensions = [kwargs['LATNAME'],kwargs['LONNAME']]
    #-- extend with date variables
    if kwargs['DATE']:
        output[kwargs['TIMENAME']] = np.atleast_1d(tim).astype('f')
        output[kwargs['VARNAME']] = np.atleast_3d(data)
        dimensions.append(kwargs['TIMENAME'])
    else:
        output[kwargs['VARNAME']] = np.copy(data)

    #-- defining the NetCDF dimensions and variables
    nc = {}
    #-- NetCDF dimensions
    for i,dim in enumerate(dimensions):
        fileID.createDimension(dim, len(output[dim]))
        nc[dim] = fileID.createVariable(dim, output[dim].dtype, (dim,))
    #-- NetCDF spatial data
    for key in [kwargs['VARNAME']]:
        nc[key] = fileID.createVariable(key, output[key].dtype,
            tuple(dimensions), fill_value=kwargs['FILL_VALUE'],
            zlib=True)
    #-- filling NetCDF variables
    for key,val in output.items():
        nc[key][:] = val.copy()

    #-- Defining attributes for longitude and latitude
    nc[kwargs['LONNAME']].long_name = 'longitude'
    nc[kwargs['LONNAME']].units = 'degrees_east'
    nc[kwargs['LATNAME']].long_name = 'latitude'
    nc[kwargs['LATNAME']].units = 'degrees_north'
    #-- Defining attributes for dataset
    nc[kwargs['VARNAME']].long_name = kwargs['LONGNAME']
    nc[kwargs['VARNAME']].units = kwargs['UNITS']
    #-- Defining attributes for date if applicable
    if kwargs['DATE']:
        nc[kwargs['TIMENAME']].long_name = kwargs['TIME_LONGNAME']
        nc[kwargs['TIMENAME']].units = kwargs['TIME_UNITS']
    #-- global variables of NetCDF file
    if kwargs['TITLE']:
        fileID.title = kwargs['TITLE']
    if kwargs['REFERENCE']:
        fileID.reference = kwargs['REFERENCE']
    #-- date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output NetCDF structure information
    logging.info(kwargs['FILENAME'])
    logging.info(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()
