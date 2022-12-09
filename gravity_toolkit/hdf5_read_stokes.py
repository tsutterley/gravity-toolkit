#!/usr/bin/env python
u"""
hdf5_read_stokes.py
Written by Tyler Sutterley (11/2021)

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
        spherical harmonics (clm,slm)
        variables (l,m,time,month)
        file description

OPTIONS:
    DATE: HDF5 file has date information
    COMPRESSION: HDF5 file is compressed or streaming as bytes
        gzip
        zip
        bytes

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org)

UPDATE HISTORY:
    Updated 11/2021: try to get more global attributes. use kwargs
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: prevent warnings with python3 compatible regex strings
    Updated 12/2020: try/except for getting variable unit attributes
        add fallback for finding HDF5 file within from zip files
        added bytes option for COMPRESSION if streaming from memory
    Updated 09/2020: use try/except for reading attributes
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
import io
import re
import gzip
import logging
import zipfile
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

def hdf5_read_stokes(filename, **kwargs):
    """
    Reads spherical harmonic data from HDF5 files

    Parameters
    ----------
    filename: HDF5 file to be opened and read

    DATE: HDF5 file has date information
    VERBOSE: will print to screen the HDF5 structure parameters
    COMPRESSION: HDF5 file is compressed or streaming as bytes
        gzip
        zip
        bytes

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
    # set default keyword arguments
    kwargs.setdefault('DATE',True)
    kwargs.setdefault('COMPRESSION',None)
    # set deprecation warning
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use harmonics.from_HDF5",
        DeprecationWarning)

    # Open the HDF5 file for reading
    if (kwargs['COMPRESSION'] == 'gzip'):
        # read gzip compressed file and extract into in-memory file object
        with gzip.open(os.path.expanduser(filename),'r') as f:
            fid = io.BytesIO(f.read())
        # set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        # rewind to start of file
        fid.seek(0)
        # read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (kwargs['COMPRESSION'] == 'zip'):
        # read zipped file and extract file into in-memory file object
        fileBasename,_ = os.path.splitext(os.path.basename(filename))
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            # first try finding a HDF5 file with same base filename
            # if none found simply try searching for a HDF5 file
            try:
                f,=[f for f in z.namelist() if re.match(fileBasename,f,re.I)]
            except:
                f,=[f for f in z.namelist() if re.search(r'\.H(DF)?5$',f,re.I)]
            # read bytes from zipfile into in-memory BytesIO object
            fid = io.BytesIO(z.read(f))
        # set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        # rewind to start of file
        fid.seek(0)
        # read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (kwargs['COMPRESSION'] == 'bytes'):
        # read as in-memory (diskless) HDF5 dataset
        fileID = h5py.File(filename, 'r')
    else:
        # read HDF5 dataset
        fileID = h5py.File(os.path.expanduser(filename), 'r')
    # allocate python dictionary for output variables
    dinput = {}
    dinput['attributes'] = {}

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # output variable keys
    h5keys = ['l','m','clm','slm']
    # Getting the data from each HDF5 variable
    # converting HDF5 objects into numpy arrays
    ll = np.array(fileID['l'][:])
    mm = np.array(fileID['m'][:])
    # Spherical harmonic files have date information
    if kwargs['DATE']:
        h5keys.extend(['time','month'])
        dinput['time'] = fileID['time'][:].copy()
        dinput['month'] = fileID['month'][:].copy()
        n_time = len(dinput['time'])
    else:
        n_time = 0

    # Restructuring input array back into matrix format
    LMAX = np.max(ll)
    MMAX = np.max(mm)

    # LMAX+1 to include LMAX (LMAX+1 elements)
    dinput['l'] = np.arange(0,LMAX+1)
    dinput['m'] = np.arange(0,MMAX+1)
    # convert input clm/slm to numpy arrays
    CLM = np.array(fileID['clm'][:])
    SLM = np.array(fileID['slm'][:])
    # size of the input grids
    n_harm, = fileID['l'].shape
    # import spherical harmonic data
    if (kwargs['DATE'] and (n_time > 1)):
        # contains multiple dates
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm],:] = CLM[lm,:]
            dinput['slm'][ll[lm],mm[lm],:] = SLM[lm,:]
    else:
        # contains either no dates or a single date
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm]] = CLM[lm]
            dinput['slm'][ll[lm],mm[lm]] = SLM[lm]

    # Getting attributes of clm/slm and included variables
    # get attributes for the included variables
    for key in h5keys:
        try:
            dinput['attributes'][key] = [
                fileID[key].attrs['units'],
                fileID[key].attrs['long_name']
                ]
        except (KeyError, AttributeError):
            pass
    # Global attributes
    for att_name in ['title','description','reference']:
        try:
            dinput['attributes'][att_name] = fileID.attrs[att_name]
        except (ValueError, KeyError, AttributeError):
            pass

    # Closing the HDF5 file
    fileID.close()

    # return the output variable
    return dinput
