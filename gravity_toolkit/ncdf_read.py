#!/usr/bin/env python
u"""
ncdf_read.py
Written by Tyler Sutterley (11/2021)

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
    COMPRESSION: netCDF4 file is compressed or streaming as bytes
        gzip
        zip
        bytes

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
    Updated 11/2021: try to get more global attributes. use kwargs
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: prevent warnings with python3 compatible regex strings
    Updated 12/2020: try/except for getting variable unit attributes
        attempt to get a standard set of attributes from each variable
        add fallback for finding netCDF4 file within from zip files
        added bytes option for COMPRESSION if streaming from memory
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
import uuid
import gzip
import logging
import netCDF4
import zipfile
import numpy as np
import warnings

def ncdf_read(filename, **kwargs):
    """
    Reads spatial data from COARDS-compliant netCDF4 files

    Parameters
    ----------
    filename: netCDF4 file to be opened and read

    DATE: netCDF4 file has date information
    VERBOSE: will print to screen the netCDF4 structure parameters
    VARNAME: z variable name in netCDF4 file
    LONNAME: longitude variable name in netCDF4 file
    LATNAME: latitude variable name in netCDF4 file
    TIMENAME: time variable name in netCDF4 file
    COMPRESSION: netCDF4 file is compressed or streaming as bytes
        gzip
        zip
        bytes

    Returns
    -------
    data: z value of dataset
    lon: longitudinal array
    lat: latitudinal array
    time: time value of dataset
    attributes: netCDF4 attributes
    """
    #-- set default keyword arguments
    kwargs.setdefault('DATE',False)
    kwargs.setdefault('VARNAME','z')
    kwargs.setdefault('LONNAME','lon')
    kwargs.setdefault('LATNAME','lat')
    kwargs.setdefault('TIMENAME','time')
    kwargs.setdefault('COMPRESSION',None)
    #-- set deprecation warning
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use spatial.from_netCDF4",
        DeprecationWarning)

    #-- Open the NetCDF4 file for reading
    if (kwargs['COMPRESSION'] == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.expanduser(filename),'r') as f:
            fileID = netCDF4.Dataset(os.path.basename(filename),memory=f.read())
    elif (kwargs['COMPRESSION'] == 'zip'):
        #-- read zipped file and extract file into in-memory file object
        fileBasename,_ = os.path.splitext(os.path.basename(filename))
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            #-- first try finding a netCDF4 file with same base filename
            #-- if none found simply try searching for a netCDF4 file
            try:
                f,=[f for f in z.namelist() if re.match(fileBasename,f,re.I)]
            except:
                f,=[f for f in z.namelist() if re.search(r'\.nc(4)?$',f)]
            #-- read bytes from zipfile as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=z.read(f))
    elif (kwargs['COMPRESSION'] == 'bytes'):
        #-- read as in-memory (diskless) netCDF4 dataset
        fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=filename.read())
    else:
        #-- read netCDF4 dataset
        fileID = netCDF4.Dataset(os.path.expanduser(filename), 'r')
    #-- create python dictionary for output variables
    dinput = {}
    dinput['attributes'] = {}

    #-- Output NetCDF file information
    logging.info(fileID.filepath())
    logging.info(list(fileID.variables.keys()))

    #-- mapping between output keys and netCDF4 variable names
    keys = ['lon','lat','data']
    nckeys = [kwargs['LONNAME'],kwargs['LATNAME'],kwargs['VARNAME']]
    if kwargs['DATE']:
        keys.append('time')
        nckeys.append(kwargs['TIMENAME'])
    #-- list of variable attributes
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','_FillValue','missing_value']
    #-- for each variable
    for key,nckey in zip(keys,nckeys):
        #-- Getting the data from each NetCDF variable
        dinput[key] = np.squeeze(fileID.variables[nckey][:].data)
        #-- Getting attributes of included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                dinput['attributes'][key][attr] = \
                    fileID.variables[nckey].getncattr(attr)
            except (KeyError,ValueError,AttributeError):
                pass

    #-- switching data array to lat/lon if lon/lat
    sz = dinput['data'].shape
    if (dinput['data'].ndim == 2) and (len(dinput['lon']) == sz[0]):
        dinput['data'] = dinput['data'].T

    #-- Global attributes
    for att_name in ['title','description','reference']:
        try:
            ncattr, = [s for s in fileID.ncattrs()
                if re.match(att_name,s,re.I)]
            dinput['attributes'][att_name] = fileID.getncattr(ncattr)
        except (ValueError, KeyError, AttributeError):
            pass
    #-- Closing the NetCDF file
    fileID.close()
    #-- return the output variable
    return dinput
