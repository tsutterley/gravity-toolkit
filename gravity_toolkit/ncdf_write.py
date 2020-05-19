#!/usr/bin/env python
u"""
ncdf_write.py
Written by Tyler Sutterley (03/2020)

Writes spatial data to COARDS-compliant NetCDF4 files

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
    CLOBBER: will overwrite an existing netCDF4 file
    VERBOSE: will print to screen the netCDF4 structure parameters
    DATE: data has date information

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
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
import netCDF4
import numpy as np

def ncdf_write(data, lon, lat, tim, FILENAME=None, VARNAME='z', LONNAME='lon',
    LATNAME='lat', TIMENAME='time', UNITS=None, LONGNAME=None, FILL_VALUE=None,
    TIME_UNITS=None, TIME_LONGNAME=None, TITLE=None, DATE=True, CLOBBER=True,
    VERBOSE=False):

    #-- setting NetCDF clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'a'

    #-- opening NetCDF file for writing
    #-- Create the NetCDF file
    fileID = netCDF4.Dataset(FILENAME, clobber, format="NETCDF4")

    #-- Defining the NetCDF dimensions
    n_time = 1 if (np.ndim(tim) == 0) else len(tim)
    fileID.createDimension(LONNAME, len(lon))
    fileID.createDimension(LATNAME, len(lat))
    fileID.createDimension(TIMENAME, n_time)

    #-- defining the NetCDF variables
    nc = {}
    #-- lat and lon
    nc[LONNAME] = fileID.createVariable(LONNAME, lon.dtype, (LONNAME,))
    nc[LATNAME] = fileID.createVariable(LATNAME, lat.dtype, (LATNAME,))
    #-- spatial data
    if (n_time > 1):
        nc[VARNAME] = fileID.createVariable(VARNAME, data.dtype,
            (LATNAME,LONNAME,TIMENAME,), fill_value=FILL_VALUE, zlib=True)
    else:
        nc[VARNAME] = fileID.createVariable(VARNAME, data.dtype,
            (LATNAME,LONNAME,), fill_value=FILL_VALUE, zlib=True)
    #-- time
    if DATE:
        nc[TIMENAME] = fileID.createVariable(TIMENAME, 'f8', (TIMENAME,))

    #-- filling NetCDF variables
    nc[LONNAME][:] = lon
    nc[LATNAME][:] = lat
    nc[VARNAME][:,:] = data
    if DATE:
        nc[TIMENAME][:] = tim

    #-- Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    #-- Defining attributes for dataset
    nc[VARNAME].long_name = LONGNAME
    nc[VARNAME].units = UNITS
    #-- Defining attributes for date if applicable
    if DATE:
        nc[TIMENAME].long_name = TIME_LONGNAME
        nc[TIMENAME].units = TIME_UNITS
    #-- global variable of NetCDF file
    if TITLE:
        fileID.TITLE = TITLE
    #-- date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    #-- Output NetCDF structure information
    if VERBOSE:
        print(FILENAME)
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()
