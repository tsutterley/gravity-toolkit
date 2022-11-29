#!/usr/bin/env python
u"""
harmonic_gradients.py
Original IDL code calc_grad.pro written by Sean Swenson
Adapted by Tyler Sutterley (10/2022)

Calculates the zonal and meridional gradients of a scalar field
from a series of spherical harmonics

CALLING SEQUENCE:
    gradient = harmonic_gradients(clm, slm, lon, lat)

INPUTS:
    clm1: cosine spherical harmonic coefficients in output units
    slm1: sine spherical harmonic coefficients in output units
    lon: longitude array for output spatial field
    lat: latitude array for output spatial field

OPTIONS:
    LMIN: Lower bound of Spherical Harmonic Degrees
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    fourier_legendre.py: Computes the Fourier coefficients of the associated
        Legendre functions

UPDATE HISTORY:
    Updated 10/2022: cleaned up program for public release
    Updated 07/2020: added function docstrings
    Updated 06/2019: using Python3 compatible division
    Updated 05/2015: code updates
    Written 05/2013
"""
from __future__ import division
import numpy as np
from gravity_toolkit.fourier_legendre import legendre_gradient

def harmonic_gradients(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None):
    """
    Calculates the gradient of a scalar field from a series of
    spherical harmonics

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

    Returns
    -------
    gradients: float
        zonal and meridional gradient fields
    """

    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Longitude in radians
    phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]
    # Colatitude in radians
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    thmax = len(np.squeeze(lat))
    phimax = len(np.squeeze(lon))

    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    mm = np.arange(0,MMAX+1)
    clm = np.zeros((LMAX+1,MMAX+1))
    slm = np.zeros((LMAX+1,MMAX+1))
    clm[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    slm[LMIN:LMAX+1,mm] = slm1[LMIN:LMAX+1,mm]
    # spherical harmonic degree and order
    ll = np.arange(0,LMAX+1)[np.newaxis, :]# lmax+1 to include lmax
    mm = np.arange(0,MMAX+1)[:, np.newaxis]# mmax+1 to include mmax

    # generate Vlm coefficients (vlm and wlm)
    vlm, wlm = legendre_gradient(LMAX, MMAX)

    dlm = np.zeros((LMAX+1,LMAX+1,2))
    #  minus sign is because lat and theta change with opposite sign
    for l in range(0,LMAX+1):
        dlm[l,:,0] = -clm[l,:]*np.sqrt((l+1.0)*l)
        dlm[l,:,1] = -slm[l,:]*np.sqrt((l+1.0)*l)

    m_even = np.arange(0,MMAX+2,2)
    m_odd = np.arange(1,MMAX,2)

    #  Calculate fourier coefficients from legendre coefficients
    d_cos = np.zeros((LMAX+1,thmax,2))
    d_sin = np.zeros((LMAX+1,thmax,2))
    cnk = np.cos(np.dot(th[:,np.newaxis],ll))
    snk = np.sin(np.dot(th[:,np.newaxis],ll))

    wtmp = np.zeros((len(m_even),LMAX+1,2))
    vtmp = np.zeros((len(m_even),LMAX+1,2))
    # m = even terms (vlm,wlm sine series)
    for n in range(0,LMAX+1):
        wtmp[:,n,0] = np.sum(wlm[:,m_even,n]*dlm[:,m_even,0],axis=0)
        wtmp[:,n,1] = np.sum(wlm[:,m_even,n]*dlm[:,m_even,1],axis=0)
        vtmp[:,n,0] = np.sum(vlm[:,m_even,n]*dlm[:,m_even,0],axis=0)
        vtmp[:,n,1] = np.sum(vlm[:,m_even,n]*dlm[:,m_even,1],axis=0)

    d_cos[m_even,:,0] = np.dot(wtmp[:,:,1],np.transpose(snk))
    d_sin[m_even,:,0] = np.dot(-wtmp[:,:,0],np.transpose(snk))
    d_cos[m_even,:,1] = np.dot(vtmp[:,:,1],np.transpose(snk))
    d_sin[m_even,:,1] = np.dot(-vtmp[:,:,0],np.transpose(snk))

    # m = odd terms (vlm,wlm cosine series)
    wtmp = np.zeros((len(m_odd),LMAX+1,2))
    vtmp = np.zeros((len(m_odd),LMAX+1,2))
    for n in range(0,LMAX+1):
        wtmp[:,n,0] = np.sum(wlm[:,m_odd,n]*dlm[:,m_odd,0],axis=0)
        wtmp[:,n,1] = np.sum(wlm[:,m_odd,n]*dlm[:,m_odd,1],axis=0)
        vtmp[:,n,0] = np.sum(vlm[:,m_odd,n]*dlm[:,m_odd,0],axis=0)
        vtmp[:,n,1] = np.sum(vlm[:,m_odd,n]*dlm[:,m_odd,1],axis=0)

    d_cos[m_odd,:,0] = np.dot(wtmp[:,:,1],np.transpose(cnk))
    d_sin[m_odd,:,0] = np.dot(-wtmp[:,:,0],np.transpose(cnk))
    d_cos[m_odd,:,1] = np.dot(vtmp[:,:,1],np.transpose(cnk))
    d_sin[m_odd,:,1] = np.dot(-vtmp[:,:,0],np.transpose(cnk))

    # Calculating cos(m*phi) and sin(m*phi)
    ccos = np.cos(np.dot(mm,phi))
    ssin = np.sin(np.dot(mm,phi))
    # Final signal recovery from fourier coefficients
    gradients = np.zeros((phimax,thmax,2))
    gradients[:,:,0] = np.dot(np.transpose(ccos), d_cos[:,:,0]) + \
        np.dot(np.transpose(ssin), d_sin[:,:,0])
    gradients[:,:,1] = np.dot(np.transpose(ccos), d_cos[:,:,1]) + \
        np.dot(np.transpose(ssin), d_sin[:,:,1])
    # return the gradient fields
    return gradients
