#!/usr/bin/env python
u"""
combine_harmonics.py
Written by Tyler Sutterley (05/2015)

Returns the spatial field for a series of spherical harmonics

CALLING SEQUENCE:
    spatial = combine_harmonics(clm1, slm1, lon, lat, LMIN=0, LMAX=60)

INPUTS:
    clm: cosine spherical harmonic coefficients in output units
    slm: sine spherical harmonic coefficients in output units
    lon: longitude
    lat: latitude

OPTIONS:
    LMIN: Lower bound of Spherical Harmonic Degrees
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    PLM: plm coefficients (if computed outside the function)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: Computes fully-normalized associated Legendre polynomials

UPDATE HISTORY:
    Updated 05/2015: added parameter MMAX for MMAX != LMAX.
    Written 05/2013
"""
import numpy as np
from gravity_toolkit.plm_holmes import plm_holmes

def combine_harmonics(clm1,slm1,lon,lat,LMIN=0,LMAX=0,MMAX=None,PLM=None):

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
        d_cos[:,k] = np.sum(PLM[:,mm,k]*clm[:,mm],axis=0)
        d_sin[:,k] = np.sum(PLM[:,mm,k]*slm[:,mm],axis=0)

    #-- Final signal recovery from fourier coefficients
    m = np.arange(0,MMAX+1)[:,np.newaxis]
    #-- Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))
    #-- summation of cosine and sine harmonics
    s = np.dot(np.transpose(ccos),d_cos) + np.dot(np.transpose(ssin),d_sin)

    #-- return output data
    return s
