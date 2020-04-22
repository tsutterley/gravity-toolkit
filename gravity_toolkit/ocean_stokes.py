#!/usr/bin/env python
u"""
ocean_stokes.py
Written by Tyler Sutterley (04/2020)

Reads a land-sea mask and converts to a series of spherical harmonics

INPUTS:
    LANDMASK: Mask file to use as input following Sutterley et al. (2020)
        1x1 degree mask distributed from UCAR as part of NCL
            https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc
        1,0.5 and 0.25 degree masks distributed from ORNL as part of ISLSCP
            https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html
    LMAX: maximum spherical harmonic degree of the output harmonics

OPTIONS:
    MMAX: maximum spherical harmonic order of the output harmonics
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

OUTPUTS:
    clm: Cosine spherical harmonic coefficients (geodesy normalization)
    slm: Sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES
    ncdf_read.py: reads a spatial grid from a netcdf file
    gen_stokes.py: converts a spatial field into a series of spherical harmonics

UPDATE HISTORY:
    Updated 04/2020 for public release
    Updated 04/2020: added option LOVE for passing load love numbers
    Updated 09/2017: added option MASK for different input masks
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 03/2015
"""
import os
import numpy as np
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.gen_stokes import gen_stokes

def ocean_stokes(LANDMASK, LMAX, MMAX=None, LOVE=None):
    #-- maximum spherical harmonic order
    MMAX = np.copy(LMAX) if MMAX is None else MMAX
    #-- Read Land-Sea Mask of specified input file
    #-- 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    #-- Open the land-sea NetCDF file for reading
    landsea = ncdf_read(os.path.expanduser(LANDMASK),
        VARNAME='LSMASK', DATE=False, ATTRIBUTES=False)
    #-- create land function
    land_function = np.zeros_like(landsea['data'])
    #-- combine land and island levels for land function
    indx,indy = np.nonzero((landsea['data'] >= 1) & (landsea['data'] <= 3))
    land_function[indx,indy] = 1.0
    #-- ocean function reciprocal of land function
    ocean_function = 1.0 - land_function
    #-- convert to spherical harmonics (1 cm w.e.)
    ocean_Ylms = gen_stokes(ocean_function, landsea['lon'], landsea['lat'],
        UNITS=1, LMIN=0, LMAX=LMAX, MMAX=MMAX, LOVE=LOVE)
    #-- return the spherical harmonic coefficients
    return ocean_Ylms
