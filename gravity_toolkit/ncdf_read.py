#!/usr/bin/env python
u"""
ncdf_read.py
Written by Tyler Sutterley (08/2020)

Reads spatial data from COARDS-compliant netCDF4 files

CALLING SEQUENCE:
    dinput = ncdf_read(filename, DATE=False, VERBOSE=False)

INPUTS:
    filename: netCDF4 file to be opened and read

OUTPUTS:
    data: z value of dataset
    lon: longitudinal array
    lat: latitudinal array
    time: time value of dataset (if specified by DATE)
    attributes: netCDF4 attributes (for variables and title)

OPTIONS:
    DATE: netCDF4 file has date information
    VERBOSE: will print to screen the netCDF4 structure parameters
    VARNAME: z variable name in netCDF4 file
    LONNAME: longitude variable name in netCDF4 file
    LATNAME: latitude variable name in netCDF4 file
    TIMENAME: time variable name in netCDF4 file
    ATTRIBUTES: netCDF4 variables contain attribute parameters
    TITLE: netCDF4 file contains title attribute parameter
    COMPRESSION: netCDF4 file is compressed using gzip or zip

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
    Updated 08/2020: flake8 compatible regular expression strings
        add options to read from gzip or zip compressed files
    Updated 07/2020: added function docstrings
    Updated 06/2020: output data as lat/lon following spatial module
        attempt to read fill value attribute and set to None if not present
    Updated 10/2019: changing Y/N flags to True/False
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 06/2018: extract fill_value and title without variable attributes
    Updated 07-09/2016: using netCDF4-python
    Updated 06/2016: using __future__ print, output filename if VERBOSE
    Updated 05/2016: will only transpose if data is 2 dimensional (not 3)
        added parameter to read the TITLE variable
    Updated 07/2015: updated read title for different cases with regex
    Updated 05/2015: added parameter TIMENAME for time variable name
    Updated 04/2015: fix attribute outputs (forgot to copy to new dictionary)
    Updated 02/2015: added copy for variable outputs
        fixes new error flag from mmap=True
    Updated 11/2014: new parameters for variable names and attributes
        all variables in a single python dictionary
    Updated 05/2014: new parameter for missing value
        new outputs: all attributes, fill value
        added try for TITLE attribute
        converting time to numpy array
    Updated 02/2014: minor update to if statements
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 01/2013: adding time variable
    Written 07/2012
"""
from __future__ import print_function

import os
import re
import gzip
import netCDF4
import zipfile
import numpy as np

def ncdf_read(filename, DATE=False, VERBOSE=False, VARNAME='z', LONNAME='lon',
    LATNAME='lat', TIMENAME='time', ATTRIBUTES=True, TITLE=True,
    COMPRESSION=None):
    """
    Reads spatial data from COARDS-compliant netCDF4 files

    Arguments
    ---------
    filename: netCDF4 file to be opened and read

    Keyword arguments
    -----------------
    DATE: netCDF4 file has date information
    VERBOSE: will print to screen the netCDF4 structure parameters
    VARNAME: z variable name in netCDF4 file
    LONNAME: longitude variable name in netCDF4 file
    LATNAME: latitude variable name in netCDF4 file
    TIMENAME: time variable name in netCDF4 file
    ATTRIBUTES: netCDF4 variables contain attribute parameters
    TITLE: netCDF4 file contains a description attribute
    COMPRESSION: netCDF4 file is compressed using gzip or zip

    Returns
    -------
    data: z value of dataset
    lon: longitudinal array
    lat: latitudinal array
    time: time value of dataset
    attributes: netCDF4 attributes
    """

    #-- Open the NetCDF4 file for reading
    if (COMPRESSION == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.expanduser(filename),'r') as f:
            fileID = netCDF4.Dataset(os.path.basename(filename),memory=f.read())
    elif (COMPRESSION == 'zip'):
        #-- read zipped file and extract file into in-memory file object
        fileBasename,fileExtension = os.path.splitext(filename)
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            #-- read bytes from zipfile as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(os.path.basename(filename),
                memory=z.read(fileBasename))
    else:
        #-- read netCDF4 dataset
        fileID = netCDF4.Dataset(os.path.expanduser(filename), 'r')
    #-- create python dictionary for output variables
    dinput = {}

    #-- Output NetCDF file information
    if VERBOSE:
        print(fileID.filepath())
        print(list(fileID.variables.keys()))

    #-- netcdf variable names
    NAMES = {}
    NAMES['lon'] = LONNAME
    NAMES['lat'] = LATNAME
    NAMES['data'] = VARNAME
    if DATE:
        NAMES['time'] = TIMENAME
    #-- for each variable
    for key,nckey in NAMES.items():
        #-- Getting the data from each NetCDF variable
        dinput[key] = fileID.variables[nckey][:]

    #-- switching data array to lat/lon if lon/lat
    sz = dinput['data'].shape
    if (dinput['data'].ndim == 2) and (len(dinput['lon']) == sz[0]):
        dinput['data'] = dinput['data'].T

    #-- getting attributes of included variables
    dinput['attributes'] = {}
    if ATTRIBUTES:
        #-- for each variable
        #-- get attributes for the included variables
        for key in NAMES.keys():
            dinput['attributes'][key] = [fileID.variables[NAMES[key]].units, \
                fileID.variables[NAMES[key]].long_name]
    #-- missing data fill value
    try:
        dinput['attributes']['_FillValue'] = fileID[VARNAME].attrs['_FillValue']
    except:
        dinput['attributes']['_FillValue'] = None
    #-- Global attribute (title of dataset)
    if TITLE:
        title, = [st for st in dir(fileID) if re.match(r'TITLE',st,re.I)]
        dinput['attributes']['title'] = getattr(fileID, title)

    #-- Closing the NetCDF file
    fileID.close()
    #-- return the output variable
    return dinput
