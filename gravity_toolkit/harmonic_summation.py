#!/usr/bin/env python
u"""
harmonic_summation.py
Written by Tyler Sutterley (01/2023)

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
    associated_legendre.py: Computes fully-normalized associated
        Legendre polynomials

UPDATE HISTORY:
    Updated 01/2023: refactored associated legendre polynomials
        added fft-based transform function
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: added function docstrings
    Updated 05/2015: added parameter MMAX for MMAX != LMAX.
    Written 05/2013
"""
import numpy as np
from gravity_toolkit.associated_legendre import plm_holmes

def harmonic_summation(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None, PLM=None):
    """
    Converts data from spherical harmonic coefficients to a spatial field

    Parameters
    ----------
    clm1: float
        cosine spherical harmonic coefficients in output units
    slm1: float
        sine spherical harmonic coefficients in output units
    lon: float
        longitude array
    lat: float
        latitude array
    LMIN: int, default 0
        Lower bound of Spherical Harmonic Degrees
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    PLM: float or NoneType, default None
        Fully-normalized associated Legendre polynomials

    Returns
    -------
    spatial: float
        spatial field
    """

    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Longitude in radians
    phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]
    # colatitude in radians
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    thmax = len(th)

    # Calculate fourier coefficients from legendre coefficients
    d_cos = np.zeros((MMAX+1,thmax))# [m,th]
    d_sin = np.zeros((MMAX+1,thmax))# [m,th]
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    mm = np.arange(0,MMAX+1)
    clm = np.zeros((LMAX+1,MMAX+1))
    slm = np.zeros((LMAX+1,MMAX+1))
    clm[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    slm[LMIN:LMAX+1,mm] = slm1[LMIN:LMAX+1,mm]
    for k in range(0,thmax):
        # summation over all spherical harmonic degrees
        d_cos[:,k] = np.sum(PLM[:,mm,k]*clm[:,mm],axis=0)
        d_sin[:,k] = np.sum(PLM[:,mm,k]*slm[:,mm],axis=0)

    # Final signal recovery from fourier coefficients
    m = np.arange(0,MMAX+1)[:,np.newaxis]
    # Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))
    # summation of cosine and sine harmonics
    s = np.dot(np.transpose(ccos),d_cos) + np.dot(np.transpose(ssin),d_sin)

    # return output data
    return s

def harmonic_transform(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None, PLM=None):
    """
    Converts data from spherical harmonic coefficients to a spatial field
    using Fast-Fourier Transforms

    Parameters
    ----------
    clm1: float
        cosine spherical harmonic coefficients in output units
    slm1: float
        sine spherical harmonic coefficients in output units
    lon: float
        longitude array
    lat: float
        latitude array
    LMIN: int, default 0
        Lower bound of Spherical Harmonic Degrees
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    PLM: float or NoneType, default None
        Fully-normalized associated Legendre polynomials

    Returns
    -------
    spatial: float
        spatial field
    """
    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # verify that longitudes cover the complete sphere
    assert np.isclose(np.min(lon), 0.0)
    assert np.isclose(np.max(lon), 360.0)
    # number of longitudinal points
    phimax = len(np.squeeze(lon))
    # colatitude in radians
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    thmax = len(th)

    # combined Ylms and Fourier coefficients (complex)
    Ylms = np.zeros((LMAX+1,MMAX+1),dtype=np.complex128)
    delta_M = np.zeros((MMAX+1,thmax),dtype=np.complex128)# [m,th]
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # Real (cosine) and imaginary (sine) components
    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    Ylms[LMIN:LMAX+1,:MMAX+1] = clm1[LMIN:LMAX+1,0:MMAX+1] - \
        slm1[LMIN:LMAX+1,0:MMAX+1]*1j
    # calculate Ylms summation for each theta band
    for k in range(0,thmax):
        # summation over all spherical harmonic degrees
        delta_M[:,k] = np.sum(PLM[:,:,k]*Ylms[:,:],axis=0)/2.0

    # output spatial field from FFT transformation
    s = np.zeros((phimax,thmax))
    # calculate fft for each theta band (over phis with axis=0)
    s[:-1,:] = 2.0*(phimax-1)*np.fft.ifft(delta_M,n=phimax-1,axis=0).real
    # complete sphere (values at 360 == values at 0)
    s[-1,:] = s[0,:]

    # return output data
    return s
