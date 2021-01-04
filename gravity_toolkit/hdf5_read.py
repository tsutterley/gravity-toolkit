#!/usr/bin/env python
u"""
hdf5_read.py
Written by Tyler Sutterley (12/2020)

Reads spatial data from HDF5 files

CALLING SEQUENCE:
    dinput = hdf5_read(filename, DATE=False, VERBOSE=False)

INPUTS:
    filename: HDF5 file to be opened and read

OUTPUTS:
    data: z value of dataset
    lon: longitudinal array
    lat: latitudinal array
    time: time value of dataset (if specified by DATE)
    attributes: HDF5 attributes (for variables and title)

OPTIONS:
    DATE: HDF5 file has date information
    VERBOSE: will print to screen the HDF5 structure parameters
    VARNAME: z variable name in HDF5 file
    LONNAME: longitude variable name in HDF5 file
    LATNAME: latitude variable name in HDF5 file
    TIMENAME: time variable name in HDF5 file
    COMPRESSION: HDF5 file is compressed or streaming as bytes
        gzip
        zip
        bytes

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org)

UPDATE HISTORY:
    Updated 12/2020: try/except for getting variable unit attributes
        attempt to get a standard set of attributes from each variable
        add fallback for finding HDF5 file within from zip files
        added bytes option for COMPRESSION if streaming from memory
    Updated 08/2020: add options to read from gzip or zip compressed files
    Updated 07/2020: added function docstrings
    Updated 06/2020: output data as lat/lon following spatial module
        attempt to read fill value attribute and set to None if not present
    Updated 10/2019: changing Y/N flags to True/False
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 06/2018: extract fill_value and title without variable attributes
    Updated 06/2016: using __future__ print function
    Updated 05/2016: will only transpose if data is 2 dimensional (not 3)
        added parameter to read the TITLE variable
    Updated 05/2015: added parameter TIMENAME for time variable name
    Updated 11/2014: got back to writing this
        in working condition with updated attributes as in netcdf equivalent
    Updated 12/2013: converted ncdf code to HDF5 code (alternative data type)
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 05/2013: converted to Python
    Updated 03/2013: converted to Octave
    Updated 01/2013: adding time variable
    Written 07/2012 for GMT and for archiving datasets
        Motivation for archival: netCDF files are much smaller than ascii
        files and more portable/transferable than IDL .sav files
        (possible to connect with geostatistics packages in R?)
"""
from __future__ import print_function

import os
import io
import re
import gzip
import h5py
import zipfile
import numpy as np

def hdf5_read(filename, DATE=False, VERBOSE=False, VARNAME='z', LONNAME='lon',
    LATNAME='lat', TIMENAME='time', COMPRESSION=None):
    """
    Reads spatial data from HDF5 files

    Arguments
    ---------
    filename: HDF5 file to be opened and read

    Keyword arguments
    -----------------
    DATE: HDF5 file has date information
    VERBOSE: will print to screen the HDF5 structure parameters
    VARNAME: z variable name in HDF5 file
    LONNAME: longitude variable name in HDF5 file
    LATNAME: latitude variable name in HDF5 file
    TIMENAME: time variable name in HDF5 file
    COMPRESSION: HDF5 file is compressed or streaming as bytes
        gzip
        zip
        bytes

    Returns
    -------
    data: z value of dataset
    lon: longitudinal array
    lat: latitudinal array
    time: time value of dataset
    attributes: HDF5 attributes
    """

    #-- Open the HDF5 file for reading
    if (COMPRESSION == 'gzip'):
        #-- read gzip compressed file and extract into in-memory file object
        with gzip.open(os.path.expanduser(filename),'r') as f:
            fid = io.BytesIO(f.read())
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (COMPRESSION == 'zip'):
        #-- read zipped file and extract file into in-memory file object
        fileBasename,_ = os.path.splitext(os.path.basename(filename))
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            #-- first try finding a HDF5 file with same base filename
            #-- if none found simply try searching for a HDF5 file
            try:
                zname,=[f for f in z.namelist() if re.match(fileBasename,f,re.I)]
            except:
                zname,=[f for f in z.namelist() if re.search('\.H(DF)?5$',f,re.I)]
            #-- read bytes from zipfile into in-memory BytesIO object
            fid = io.BytesIO(z.read(zname))
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (COMPRESSION == 'bytes'):
        #-- read as in-memory (diskless) HDF5 dataset
        fileID = h5py.File(filename, 'r')
    else:
        #-- read HDF5 dataset
        fileID = h5py.File(os.path.expanduser(filename), 'r')
    #-- allocate python dictionary for output variables
    dinput = {}
    dinput['attributes'] = {}

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- mapping between output keys and HDF5 variable names
    keys = ['lon','lat','data']
    h5keys = [LONNAME,LATNAME,VARNAME]
    if DATE:
        keys.append('time')
        h5keys.append(TIMENAME)

    #-- list of variable attributes
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','_FillValue','missing_value']
    #-- for each variable
    for key,h5key in zip(keys,h5keys):
        #-- Getting the data from each HDF5 variable
        dinput[key] = fileID[h5key][:].copy()
        #-- Getting attributes of included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            try:
                dinput['attributes'][key][attr] = fileID[h5key].attrs[attr]
            except (KeyError, AttributeError):
                pass

    #-- switching data array to lat/lon if lon/lat
    sz = dinput['data'].shape
    if (dinput['data'].ndim == 2) and (len(dinput['lon']) == sz[0]):
        dinput['data'] = dinput['data'].T

    #-- Global attribute description
    try:
        dinput['attributes']['title'] = fileID.attrs['description']
    except (KeyError, AttributeError):
        pass

    #-- Closing the HDF5 file
    fileID.close()
    return dinput
