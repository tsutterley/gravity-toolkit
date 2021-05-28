#!/usr/bin/env python
u"""
gen_harmonics.py
Written by Tyler Sutterley (05/2021)
Converts data from the spatial domain to spherical harmonic coefficients

Differs from the gen_stokes() function as it does not
    compute the solid earth elastic terms

CALLING SEQUENCE:
    Ylms = gen_harmonics(data, lon, lat, LMIN=0, LMAX=60, DDEG=0.5)

INPUTS:
    data: data magnitude
    lon: longitude array
    lat: latitude array

OUTPUTS:
    Ylms: harmonics object
        clm: 4-pi normalized cosine spherical harmonic coefficients
        slm: 4-pi normalized sine spherical harmonic coefficients
        l: spherical harmonic degree to LMAX
        m: spherical harmonic order to MMAX

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    PLM: input fully normalized associated Legendre polynomials

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5

REFERENCE:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)

UPDATE HISTORY:
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Updated 04/2020: include degrees and orders in output dictionary
    Updated 10/2017: updated comments and cleaned up code
    Updated 08/2017: Using Holmes and Featherstone relation for Plms
    Written 09/2016
"""
import numpy as np
import gravity_toolkit.harmonics
from gravity_toolkit.plm_holmes import plm_holmes

def gen_harmonics(data, lon, lat, LMAX=60, MMAX=None, PLM=0):
    """
    Converts data from the spatial domain to spherical harmonic coefficients

    Arguments
    ---------
    data: data magnitude
    lon: longitude array
    lat: latitude array

    Keyword arguments
    -----------------
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    PLM: input Legendre polynomials

    Returns
    -------
    clm: cosine spherical harmonic coefficients
    slm: sine spherical harmonic coefficients
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX
    """

    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- dimensions of the longitude and latitude arrays
    nlon = np.int64(len(lon))
    nlat = np.int64(len(lat))

    #-- grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    #-- longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    #-- colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    #-- convert latitude and longitude to float if integers
    lon = lon.astype(np.float64)
    lat = lat.astype(np.float64)
    #-- reformatting longitudes to range 0:360 (if previously -180:180)
    if np.count_nonzero(lon < 0):
        lon[lon < 0] += 360.0
    #-- calculate longitude and colatitude arrays in radians
    phi = np.reshape(lon,(1,nlon))*np.pi/180.0#-- reshape to 1xnlon
    th = (90.0 - np.squeeze(lat))*np.pi/180.0#-- remove singleton dimensions

    #-- reforming data to lonXlat if input latXlon
    sz = np.shape(data)
    dinput = np.transpose(data) if (sz[0] == nlat) else np.copy(data)

    #-- Calculating cos/sin of phi arrays (output [m,phi])
    #-- LMAX+1 as there are LMAX+1 elements between 0 and LMAX
    m = np.arange(MMAX+1)[:, np.newaxis]
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))

    #-- Multiplying sin(th) with differentials of theta and phi
    #-- to calculate the integration factor at each latitude
    int_fact = np.sin(th)*dphi*dth
    coeff = 1.0/(4.0*np.pi)

    #-- Calculates plms using Holmes and Featherstone (2002) recursion relation
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    #-- added option to precompute plms to improve computational speed
    if (np.ndim(PLM) == 0):
        plmout,dplm = plm_holmes(LMAX,np.cos(th))
    else:
        plmout = PLM

    #-- Multiply plms by integration factors [sin(theta)*dtheta*dphi]
    #-- truncate plms to maximum spherical harmonic order if MMAX < LMAX
    m = np.arange(MMAX+1)
    for j in range(0,nlat):
        plm[:,m,j] = plmout[:,m,j]*int_fact[j]

    #-- Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    #-- Multiplying gridded data with sin/cos of m#phis (output [m,theta])
    #-- This will sum through all phis in the dot product
    dcos = np.dot(ccos,dinput)
    dsin = np.dot(ssin,dinput)
    for l in range(0,LMAX+1):
        mm = np.min([MMAX,l])#-- truncate to MMAX if specified (if l > MMAX)
        m = np.arange(0,mm+1)#-- mm+1 elements between 0 and mm
        #-- Summing product of plms and data over all latitudes
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        #-- convert to output normalization (4-pi normalized harmonics)
        Ylms.clm[l,m] = coeff*yclm[l,m]
        Ylms.slm[l,m] = coeff*yslm[l,m]

    #-- return the output spherical harmonics object
    return Ylms
