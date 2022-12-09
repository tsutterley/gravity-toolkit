#!/usr/bin/env python
u"""
ncdf_read_stokes.py
Written by Tyler Sutterley (11/2021)

Reads spherical harmonic data from netCDF4 files

CALLING SEQUENCE:
    Ylms = ncdf_read_stokes(filename, DATE=True, VERBOSE=False)

INPUTS:
    filename: netCDF4 file to be opened and read

OUTPUTS:
    clm: Cosine Stokes Coefficient
    slm: Sine Stokes Coefficient
    l: degree (l)
    m: order (m)
    time: time of measurement (if specified by DATE)
    month: GRACE/GRACE-FO month (if specified by DATE)
    attributes: netCDF4 attributes for:
        spherical harmonics (clm,slm)
        variables (l,m,time,month)
        file title

OPTIONS:
    DATE: netCDF4 file has date information
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
        add fallback for finding netCDF4 file within from zip files
        added bytes option for COMPRESSION if streaming from memory
    Updated 09/2020: use try/except for reading attributes
    Updated 08/2020: flake8 compatible regular expression strings
        add options to read from gzip or zip compressed files
    Updated 07/2020: added function docstrings
    Updated 03/2020: added ATTRIBUTES option to check if file has attributes
    Updated 10/2019: changing Y/N flags to True/False
    Updated 03/2019: print variables keys in list for Python3 compatibility
    Updated 09/2016: slicing of clm and slm on numpy arrays not netcdf variables
    Updated 07/2016: using netCDF4-python
    Updated 06/2016: using __future__ print, output filename if VERBOSE
    Updated 02/2016: capitalized LMAX and MMAX variables to match other programs
    Updated 07/2015: updated read title for different cases with regex
    Updated 06/2015: can input a single netcdf with multiple dates
    Updated 05/2015: minor change for MMAX != LMAX
    Updated 02/2015: simplified attributes with for loop
    Updated 01/2015: added copy for variable outputs
        fixes new error flag from mmap=True
    Updated 11/2014: all variables in a single python dictionary
    Updated 05/2014: converted time and month to numpy arrays
    Updated 05/2014: output all attributes under single variable
        added try for TITLE attribute
    Updated 02/2014: minor update to if statements
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 07/2013: switched from Scientific Python to Scipy
    Updated 03/2013: switched I/O to column arrays instead of matrix
    Written 07/2012
"""
from __future__ import print_function

import os
import re
import uuid
import gzip
import logging
import zipfile
import numpy as np
import warnings

# attempt imports
try:
    import netCDF4
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("netCDF4 not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

def ncdf_read_stokes(filename, **kwargs):
    """
    Reads spherical harmonic data from netCDF4 files

    Parameters
    ----------
    filename: netCDF4 file to be opened and read

    DATE: netCDF4 file has date information
    COMPRESSION: netCDF4 file is compressed or streaming as bytes
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
    attributes: netCDF4 attributes for variables and file
    """
    # set default keyword arguments
    kwargs.setdefault('DATE',True)
    kwargs.setdefault('COMPRESSION',None)
    # set deprecation warning
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use harmonics.from_netCDF4",
        DeprecationWarning)

    # Open the NetCDF4 file for reading
    if (kwargs['COMPRESSION'] == 'gzip'):
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.expanduser(filename),'r') as f:
            fileID = netCDF4.Dataset(os.path.basename(filename),memory=f.read())
    elif (kwargs['COMPRESSION'] == 'zip'):
        # read zipped file and extract file into in-memory file object
        fileBasename,_ = os.path.splitext(os.path.basename(filename))
        with zipfile.ZipFile(os.path.expanduser(filename)) as z:
            # first try finding a netCDF4 file with same base filename
            # if none found simply try searching for a netCDF4 file
            try:
                f,=[f for f in z.namelist() if re.match(fileBasename,f,re.I)]
            except:
                f,=[f for f in z.namelist() if re.search(r'\.nc(4)?$',f)]
            # read bytes from zipfile as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=z.read(f))
    elif (kwargs['COMPRESSION'] == 'bytes'):
        # read as in-memory (diskless) netCDF4 dataset
        fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=filename.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(os.path.expanduser(filename), 'r')
    # create python dictionary for output variables
    dinput = {}
    dinput['attributes'] = {}

    # Output NetCDF file information
    logging.info(fileID.filepath())
    logging.info(list(fileID.variables.keys()))

    # Getting the data from each NetCDF variable
    # converting NetCDF objects into numpy arrays
    nckeys = ['l','m','clm','slm']
    ll = fileID.variables['l'][:].copy()
    mm = fileID.variables['m'][:].copy()
    clm = fileID.variables['clm'][:].copy()
    slm = fileID.variables['slm'][:].copy()
    # read date variables if specified
    if kwargs['DATE']:
        nckeys.extend(['time','month'])
        dinput['time'] = fileID.variables['time'][:].copy()
        dinput['month'] = fileID.variables['month'][:].copy()
        n_time = len(dinput['time'])
    else:
        n_time = 0

    # Restructuring input array back into matrix format
    LMAX = np.max(ll)
    MMAX = np.max(mm)
    # output spherical harmonic degree and order
    # LMAX+1 to include LMAX (LMAX+1 elements)
    dinput['l'] = np.arange(0,LMAX+1)
    dinput['m'] = np.arange(0,MMAX+1)
    # number of harmonics
    n_harm, = fileID.variables['l'].shape
    # import spherical harmonic data
    if (kwargs['DATE'] and (n_time > 1)):
        # contains multiple dates
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm],:] = clm[lm,:]
            dinput['slm'][ll[lm],mm[lm],:] = slm[lm,:]
    else:
        # contains either no dates or a single date
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm]] = clm[lm]
            dinput['slm'][ll[lm],mm[lm]] = slm[lm]

    # Getting attributes of clm/slm and included variables
    # get attributes for the included variables
    for key in nckeys:
        try:
            dinput['attributes'][key] = [
                fileID.variables[key].units,
                fileID.variables[key].long_name
                ]
        except (KeyError,ValueError,AttributeError):
            pass
    # Global attributes
    for att_name in ['title','description','reference']:
        try:
            ncattr, = [s for s in fileID.ncattrs()
                if re.match(att_name,s,re.I)]
            dinput['attributes'][att_name] = fileID.getncattr(ncattr)
        except (ValueError, KeyError, AttributeError):
            pass

    # Closing the NetCDF file
    fileID.close()

    # return output variable
    return dinput
