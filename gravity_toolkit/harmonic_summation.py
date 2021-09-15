#!/usr/bin/env python
u"""
harmonic_summation.py
Written by Tyler Sutterley (07/2020)

Returns the spatial field for a series of spherical harmonics

CALLING SEQUENCE:
    spatial = harmonic_summation(clm1, slm1, lon, lat, LMIN=0, LMAX=60)

INPUTS:
    clm1: cosine spherical harmonic coefficients in output units
    slm1: sine spherical harmonic coefficients in output units
    lon: longitude array for output spatial field
    lat: latitude array for output spatial field

OPTIONS:
    LMIN: Lower bound of Spherical Harmonic Degrees
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    PLM: Fully-normalized associated Legendre polynomials

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: Computes fully-normalized associated Legendre polynomials

UPDATE HISTORY:
    Updated 09/2020: added influence of Earth oblateness in option
    Updated 07/2020: added function docstrings
    Updated 05/2015: added parameter MMAX for MMAX != LMAX.
    Written 05/2013
"""
import numpy as np
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.units import units

def harmonic_summation(clm1,slm1,lon,lat,LMIN=0,LMAX=0,MMAX=None,PLM=None,ellps=False):
    """
    Converts data from spherical harmonic coefficients to a spatial field

    Arguments
    ---------
    clm1: cosine spherical harmonic coefficients in output units
    slm1: sine spherical harmonic coefficients in output units
    lon: longitude array
    lat: latitude array

    Keyword arguments
    -----------------
    LMIN: Lower bound of Spherical Harmonic Degrees
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    PLM: Fully-normalized associated Legendre polynomials
    ellps: Consideration of the Earth oblateness

    Returns
    -------
    spatial: spatial field
    """

    #-- if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- Longitude in radians
    phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]
    #-- Colatitude in radians
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    thmax = len(th)

     #-- create a unit object to get constants
    unit = units(lmax=LMAX)

    #--  Calculate fourier coefficients from legendre coefficients
    d_cos = np.zeros((MMAX+1,thmax))#-- [m,th]
    d_sin = np.zeros((MMAX+1,thmax))#-- [m,th]
    if PLM is None:
        #-- if plms are not pre-computed: calculate Legendre polynomials
        PLM,dPLM = plm_holmes(LMAX,np.cos(th))

    #-- Truncating harmonics to degree and order LMAX
    #-- removing coefficients below LMIN and above MMAX
    mm = np.arange(0,MMAX+1)
    clm = np.zeros((LMAX+1,MMAX+1))
    slm = np.zeros((LMAX+1,MMAX+1))
    clm[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    slm[LMIN:LMAX+1,mm] = slm1[LMIN:LMAX+1,mm]
    for k in range(0,thmax):
        #-- summation over all spherical harmonic degrees
        if ellps:
            d_cos[:,k] = np.sum(PLM[:,mm,k]*clm[:,mm],axis=0)*(np.sqrt(1 - (2*unit.flat - unit.flat**2)*np.sin(th[k])**2)/(1 - unit.flat))**(unit.l + 2)
            d_sin[:,k] = np.sum(PLM[:,mm,k]*slm[:,mm],axis=0)*(np.sqrt(1 - (2*unit.flat - unit.flat**2)*np.sin(th[k])**2)/(1 - unit.flat))**(unit.l + 2)
        else:
            d_cos[:, k] = np.sum(PLM[:, mm, k] * clm[:, mm], axis=0)
            d_sin[:, k] = np.sum(PLM[:, mm, k] * slm[:, mm], axis=0)

    #-- Final signal recovery from fourier coefficients
    m = np.arange(0,MMAX+1)[:,np.newaxis]
    # -- Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m, phi))
    ssin = np.sin(np.dot(m, phi))
    #-- summation of cosine and sine harmonics
    s = np.dot(np.transpose(ccos),d_cos) + np.dot(np.transpose(ssin),d_sin)

    #-- return output data
    return s
