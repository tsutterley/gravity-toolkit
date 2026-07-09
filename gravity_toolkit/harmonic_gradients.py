#!/usr/bin/env python
u"""
harmonic_gradients.py
Original IDL code calc_grad.pro written by Sean Swenson
Adapted by Tyler Sutterley (07/2026)

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
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
    Updated 03/2023: improve typing for variables in docstrings
        added geostrophic currents program from Wahr et al. (2002)
    Updated 10/2022: cleaned up program for public release
    Updated 07/2020: added function docstrings
    Updated 06/2019: using Python3 compatible division
    Updated 05/2015: code updates
    Written 05/2013
"""
from __future__ import division
import numpy as np
from gravity_toolkit.fourier_legendre import legendre_gradient
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.units import units

def harmonic_gradients(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None):
    """
    Calculates the gradient of a scalar field from a series of
    spherical harmonics :cite:p:`Driscoll:1994bp`

    Parameters
    ----------
    clm1: np.ndarray
        cosine spherical harmonic coefficients in output units
    slm1: np.ndarray
        sine spherical harmonic coefficients in output units
    lon: np.ndarray
        longitude array
    lat: np.ndarray
        latitude array
    LMIN: int, default 0
        Lower bound of Spherical Harmonic Degrees
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders

    Returns
    -------
    gradients: np.ndarray
        zonal and meridional gradient fields
    """

    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Longitude in radians
    phi = np.radians(np.squeeze(lon))
    # Colatitude in radians
    th = np.radians(90.0 - np.squeeze(lat))
    thmax = len(th)

    # spherical harmonic degree and order
    ll = np.arange(0,LMAX+1)# lmax+1 to include lmax
    mm = np.arange(0,MMAX+1)# mmax+1 to include mmax
    # real (cosine) and imaginary (sine) components
    Ylm = np.zeros((LMAX+1, MMAX+1), dtype=np.complex128)
    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    Ylm.real[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm].copy()
    Ylm.imag[LMIN:LMAX+1,mm] = -slm1[LMIN:LMAX+1,mm].copy()
    dlm = np.einsum("l...,lm...->lm", np.sqrt((ll+1.0)*ll), -1j*Ylm)

    # generate Vlm coefficients (vlm and wlm)
    Vlmk, Wlmk = legendre_gradient(LMAX, MMAX)
    # even and odd spherical harmonic orders
    m_even = np.arange(0,MMAX+2,2)
    m_odd = np.arange(1,MMAX,2)

    # Euler's formula for theta * k and m * phi
    k_th = np.exp(1j * np.einsum("h...,k...->kh...", th, ll))
    m_phi = np.exp(1j * np.einsum("m...,p...->mp...", mm, phi))
    # Calculate fourier coefficients from legendre coefficients
    d = np.zeros((LMAX+1,thmax,2), dtype=np.complex128)
    wtmp = np.einsum("lmk...,lm...->mk", Wlmk, dlm)
    vtmp = np.einsum("lmk...,lm...->mk", Vlmk, dlm)
    d[m_even,:,0] = np.einsum("mk...,kh...->mh", wtmp[m_even,:], k_th.imag)
    d[m_even,:,1] = np.einsum("mk...,kh...->mh", vtmp[m_even,:], k_th.imag)
    d[m_odd,:,0] = np.einsum("mk...,kh...->mh", wtmp[m_odd,:], k_th.real)
    d[m_odd,:,1] = np.einsum("mk...,kh...->mh", vtmp[m_odd,:], k_th.real)
    # calculate the zonal and meridional gradients of the scalar field
    gradients = np.einsum("mp...,mhd...->phd...", m_phi, d)
    # return the gradient fields and drop imaginary component
    return gradients.real

def geostrophic_currents(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None, RAD=0,
    DENSITY=1.035, LOVE=None, PLM=None):
    r"""
    Converts data from spherical harmonic coefficients to spatial
    fields of approximate ocean geostrophic currents following
    :cite:p:`Wahr:1998hy,Wahr:2002ie`

    Parameters
    ----------
    clm1: np.ndarray
        cosine spherical harmonic coefficients
    slm1: np.ndarray
        sine spherical harmonic coefficients
    lon: np.ndarray
        longitude array
    lat: np.ndarray
        latitude array
    LMIN: int, default 0
        Lower bound of Spherical Harmonic Degrees
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    RAD: float, default 0.0
        Gaussian smoothing radius (km)
    LMAX: int, default 0
        Upper bound of Spherical Harmonic Degrees
    DENSITY: float, default 1.035
        Average density of seawater at depth in g/cm\ :sup:`3`
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    PLM: np.ndarray or NoneType, default None
        Fully-normalized associated Legendre polynomials

    Returns
    -------
    currents: np.ndarray
        zonal and meridional current fields [cm/s]
    """

    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Longitude in radians
    phi = np.radians(np.squeeze(lon))
    phmax = len(lon)
    # colatitude in radians
    th = np.radians(90.0 - np.squeeze(lat))
    thmax = len(th)

    # Gaussian Smoothing
    if (RAD != 0):
        wl = 2.0*np.pi*gauss_weights(RAD, LMAX)
    else:
        # else = 1
        wl = np.ones((LMAX+1))

    # Setting units factor for output
    # extract arrays of kl, hl, and ll Love Numbers
    factors = units(lmax=LMAX).harmonic(*LOVE)
    coeff = factors.g_wmo*factors.rho_e/(6.0*factors.omega*DENSITY)

    # if plms are not pre-computed: calculate Legendre polynomials
    if PLM is None:
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # smooth harmonics and convert to output units
    clm = np.zeros((LMAX+1, MMAX+1, 2))
    slm = np.zeros((LMAX+1, MMAX+1, 2))
    # zonal flow harmonics (equation 3)
    # differentiating Legendre polynomials with respect to longitude
    for l in range(1, LMAX):
        # truncate to degree and order
        mm = np.arange(0, np.min([l,MMAX])+1)
        temp1 = (l - 1.0)/(1.0 + LOVE.kl[l-1]) * \
            np.sqrt((l**2 - mm**2)*(2.0*l - 1.0)/(2.0*l + 1))
        temp2 = (l + 2.0)/(1.0 + LOVE.kl[l+1]) * \
            np.sqrt(((l+1)**2 - mm**2)*(2.0*l + 3.0)/(2.0*l + 1))
        clm[l,mm,0] = coeff*wl[l]*(temp1*clm1[l-1,mm] - temp2*clm1[l+1,mm])
        slm[l,mm,0] = coeff*wl[l]*(temp1*slm1[l-1,mm] - temp2*slm1[l+1,mm])
    # meridional flow harmonics (equation 4)
    # differentiating Legendre polynomials with respect to colatitude
    for l in range(0, LMAX+1):
        # truncate to degree and order
        mm = np.arange(0, np.min([l,MMAX])+1)
        temp = mm*(2.0*l + 1.0)/(1.0 + LOVE.kl[l])
        clm[l,mm,1] = -coeff*wl[l]*temp*slm1[l,mm]
        slm[l,mm,1] = coeff*wl[l]*temp*clm1[l,mm]

    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    mm = np.arange(0, MMAX+1)
    # real (cosine) and imaginary (sine) components
    Ylm = clm[LMIN:LMAX+1,:MMAX+1,:] - 1j * slm[LMIN:LMAX+1,:MMAX+1,:]
    # convolve legendre polynomials and truncate to degree and order
    iint = 1.0/(np.cos(th)*np.sin(th))
    plm = np.einsum("h...,lmh...->lmh...", iint, PLM[LMIN:LMAX+1,:MMAX+1,:])
    # summation over all spherical harmonic degrees
    pconv = np.einsum("lmh...,lmd...->mhd...", plm, Ylm)

    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    m_phi = np.exp(1j * np.einsum("m...,p...->mp...", mm, phi))
    # output geostrophic current fields
    currents = np.empty((phmax,thmax,2))
    # summation of cosine and sine harmonics
    currents[:] = np.einsum("mp...,mhd...->phd...", m_phi, pconv)

    # return the current fields and drop imaginary component
    return currents.real
