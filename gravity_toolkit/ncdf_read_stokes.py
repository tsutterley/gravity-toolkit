#!/usr/bin/env python
u"""
ncdf_read_stokes.py
Written by Tyler Sutterley (07/2020)

Reads spherical harmonic data from netCDF4 files

CALLING SEQUENCE:
    file_inp = ncdf_read_stokes(filename, DATE=True, VERBOSE=False)

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
        spherical harmonics (clm,slm), variables (l,m,time,month), and title

OPTIONS:
    DATE: netCDF4 file has date information
    ATTRIBUTES: netCDF4 variables contain attribute parameters
    VERBOSE: will print to screen the netCDF4 structure parameters

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)

UPDATE HISTORY:
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

import netCDF4
import numpy as np
import re

def ncdf_read_stokes(filename, DATE=True, ATTRIBUTES=True, VERBOSE=False):
    """
    Reads spherical harmonic data from netCDF4 files

    Arguments
    ---------
    filename: netCDF4 file to be opened and read

    Keyword arguments
    -----------------
    DATE: netCDF4 file has date information
    ATTRIBUTES: netCDF4 variables contain attribute parameters
    VERBOSE: will print to screen the netCDF4 structure parameters

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

    #-- Open the NetCDF file for reading
    fileID = netCDF4.Dataset(filename, 'r')
    #-- create python dictionary for output variables
    dinput = {}
    #-- create python dictionary for variable attributes
    attributes = {}

    #-- Output NetCDF file information
    if VERBOSE:
        print(fileID.filepath())
        print(list(fileID.variables.keys()))

    #-- Getting the data from each NetCDF variable
    #-- converting NetCDF objects into numpy arrays
    ll = fileID.variables['l'][:].copy()
    mm = fileID.variables['m'][:].copy()
    clm = fileID.variables['clm'][:].copy()
    slm = fileID.variables['slm'][:].copy()
    #-- save date variables if specified
    if DATE:
        dinput['time'] = fileID.variables['time'][:].copy()
        dinput['month'] = fileID.variables['month'][:].copy()
        n_time = len(dinput['time'])
    else:
        n_time = 0

    #-- Restructuring input array back into matrix format
    LMAX = np.max(ll)
    MMAX = np.max(mm)
    #-- output spherical harmonic degree and order
    #-- LMAX+1 to include LMAX (LMAX+1 elements)
    dinput['l'] = np.arange(0,LMAX+1)
    dinput['m'] = np.arange(0,MMAX+1)
    #-- number of harmonics
    n_harm, = fileID.variables['l'].shape
    #-- import spherical harmonic data
    if (DATE and (n_time > 1)):
        #-- contains multiple dates
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1,n_time))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm],:] = clm[lm,:]
            dinput['slm'][ll[lm],mm[lm],:] = slm[lm,:]
    else:
        #-- contains either no dates or a single date
        dinput['clm'] = np.zeros((LMAX+1,MMAX+1))
        dinput['slm'] = np.zeros((LMAX+1,MMAX+1))
        for lm in range(n_harm):
            dinput['clm'][ll[lm],mm[lm]] = clm[lm]
            dinput['slm'][ll[lm],mm[lm]] = slm[lm]

    #-- Getting attributes of clm/slm and included variables
    if ATTRIBUTES:
        #-- for each variable
        #-- get attributes for the included variables
        for key in dinput.keys():
            attributes[key] = [fileID.variables[key].units, \
                fileID.variables[key].long_name]
        #-- put attributes in output python dictionary
        dinput['attributes'] = attributes
        #-- Global attribute (title of dataset)
        rx = re.compile('TITLE',re.IGNORECASE)
        title, = [st for st in dir(fileID) if rx.match(st)]
        dinput['attributes']['title'] = getattr(fileID, title)

    #-- Closing the NetCDF file
    fileID.close()

    #-- return output variable
    return dinput
