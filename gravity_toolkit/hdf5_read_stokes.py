#!/usr/bin/env python
u"""
hdf5_read_stokes.py
Written by Tyler Sutterley (08/2020)

Reads spherical harmonic data from HDF5 files

CALLING SEQUENCE:
    Ylms = hdf5_read_stokes(filename, DATE=True, VERBOSE=False)

INPUTS:
    filename: HDF5 file to be opened and read

OUTPUTS:
    clm: cosine spherical harmonic coefficients
    slm: sine spherical harmonic coefficients
    l: degree (l)
    m: order (m)
    time: time of measurement (if specified by DATE)
    month: GRACE/GRACE-FO month (if specified by DATE)
    attributes: HDF5 attributes for:
        spherical harmonics (clm,slm), variables (l,m,time,month), and title

OPTIONS:
    DATE: HDF5 file has date information
    VERBOSE: will print to screen the HDF5 structure parameters
    ATTRIBUTES: HDF5 variables contain attribute parameters
    COMPRESSION: HDF5 file is compressed using gzip or zip

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org)

UPDATE HISTORY:
    Updated 08/2020: add options to read from gzip or zip compressed files
    Updated 07/2020: added function docstrings
    Updated 03/2020: added ATTRIBUTES option to check if file has attributes
    Updated 10/2019: changing Y/N flags to True/False.  check if time is array
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 06/2016: using __future__ print function
    Updated 02/2016: capitalized LMAX and MMAX variables to match other programs
    Updated 05/2015: minor change for MMAX != LMAX
    Updated 11/2014: got back to writing this
        in working condition with updated attributes as in netcdf equivalent
    Updated 12/2013: converted ncdf code to HDF5 code (alternative data type)
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 03/2013: switched I/O to column arrays instead of matrix
    Written 07/2012
"""
from __future__ import print_function

import os
import re
import io
import gzip
import h5py
import zipfile
import numpy as np

def hdf5_read_stokes(filename, DATE=True, VERBOSE=False, ATTRIBUTES=True,
    COMPRESSION=None):
    """
    Reads spherical harmonic data from HDF5 files

    Arguments
    ---------
    filename: HDF5 file to be opened and read

    Keyword arguments
    -----------------
    DATE: HDF5 file has date information
    VERBOSE: will print to screen the HDF5 structure parameters
    ATTRIBUTES: HDF5 variables contain attribute parameters
    COMPRESSION: HDF5 file is compressed using gzip or zip

    Returns
    -------
    clm: cosine spherical harmonic coefficients
    slm: sine spherical harmonic coefficients
    l: degree
    m: order
    time: time of measurement
    month: GRACE/GRACE-FO month
    attributes: HDF5 attributes for variables and file
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
        fileBasename,fileExtension = os.path.splitext(filename)
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            #-- read bytes from zipfile into in-memory BytesIO object
            fid = io.BytesIO(z.read(fileBasename))
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    else:
        #-- read HDF5 dataset
        fileID = h5py.File(os.path.expanduser(filename), 'r')
    #-- allocate python dictionary for output variables
    dinput = {}

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- Getting the data from each HDF5 variable
    #-- converting HDF5 objects into numpy arrays
    ll = np.array(fileID['l'][:])
    mm = np.array(fileID['m'][:])
    #-- Spherical harmonic files have date information
    if DATE:
        dinput['time'] = fileID['time'][:].copy()
        dinput['month'] = fileID['month'][:].copy()
        n_time = len(dinput['time'])
    else:
        n_time = 0

    #-- Restructuring input array back into matrix format
    LMAX = np.max(ll)
    MMAX = np.max(mm)
    #-- LMAX+1 to include LMAX (LMAX+1 elements)
    dinput['l'] = np.arange(0,LMAX+1)
    dinput['m'] = np.arange(0,MMAX+1)
    #-- convert input clm/slm to numpy arrays
    CLM = np.array(fileID['clm'][:])
    SLM = np.array(fileID['slm'][:])
    #-- size of the input grids
    n_harm, = fileID['l'].shape
    #-- import spherical harmonic data
    if (DATE and (n_time > 1)):
        #-- contains multiple dates
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm],:] = CLM[lm,:]
            dinput['slm'][ll[lm],mm[lm],:] = SLM[lm,:]
    else:
        #-- contains either no dates or a single date
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm]] = CLM[lm]
            dinput['slm'][ll[lm],mm[lm]] = SLM[lm]

    #-- Getting attributes of clm/slm and included variables
    if ATTRIBUTES:
        dinput['attributes'] = {}
        dinput['attributes']['l'] = [fileID['l'].attrs['units'], \
            fileID['l'].attrs['long_name']]
        dinput['attributes']['m'] = [fileID['m'].attrs['units'], \
            fileID['m'].attrs['long_name']]
        dinput['attributes']['clm'] = [fileID['clm'].attrs['units'], \
            fileID['clm'].attrs['long_name']]
        dinput['attributes']['slm'] = [fileID['slm'].attrs['units'], \
            fileID['slm'].attrs['long_name']]
        #-- time attributes
        if DATE:
            dinput['attributes']['time'] = [fileID['time'].attrs['units'], \
                fileID['time'].attrs['long_name']]
            dinput['attributes']['month'] = [fileID['month'].attrs['units'], \
                fileID['month'].attrs['long_name']]
        #-- Global attribute description
        dinput['attributes']['title'] = fileID.attrs['description']

    #-- Closing the HDF5 file
    fileID.close()

    #-- return the output variable
    return dinput
