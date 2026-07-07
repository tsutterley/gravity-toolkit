#!/usr/bin/env python
u"""
harmonic_summation.py
Written by Tyler Sutterley (07/2026)

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
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    units.py: class for converting spherical harmonic data to specific units

UPDATE HISTORY:
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
    Updated 04/2023: allow love numbers to be None for custom units case
    Updated 03/2023: allow units inputs to be strings for named types
        improve typing for variables in docstrings
        minor refactor in line ordering for readability
    Updated 02/2023: set custom units as top option in if/else statements
    Updated 01/2023: refactored associated legendre polynomials
        added wrapper function for smoothing and converting to output units
        added fft-based transform function
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: added function docstrings
    Updated 05/2015: added parameter MMAX for MMAX != LMAX.
    Written 05/2013
"""
import numpy as np
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.units import units

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

    # longitude in radians
    phi = np.radians(np.squeeze(lon))
    # colatitude in radians
    th = np.radians(90.0 - np.squeeze(lat))

    # if plms are not pre-computed: calculate Legendre polynomials
    if PLM is None:
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # spherical harmonic order
    mm = np.arange(0,MMAX+1)# mmax+1 to include mmax
    # real (cosine) and imaginary (sine) components
    Ylm = np.zeros((LMAX+1, MMAX+1), dtype=np.complex128)
    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    Ylm.real[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    Ylm.imag[LMIN:LMAX+1,mm] = -slm1[LMIN:LMAX+1,mm]
    # Calculate fourier coefficients from legendre coefficients
    # summation over all spherical harmonic degrees
    pconv = np.einsum("lmh...,lm...->mh...", PLM, Ylm)
    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    m_phi = np.exp(1j * np.einsum("m...,p...->mp...", mm, phi))
    # summation of cosine and sine harmonics
    spatial = np.einsum("mp...,mh...->ph...", m_phi, pconv)
    # return output data and drop imaginary component
    return spatial.real

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
    th = np.radians(90.0 - np.squeeze(lat))
    thmax = len(th)

    # if plms are not pre-computed: calculate Legendre polynomials
    if PLM is None:
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # spherical harmonic order
    mm = np.arange(0,MMAX+1)# mmax+1 to include mmax
    # real (cosine) and imaginary (sine) components
    Ylm = np.zeros((LMAX+1, MMAX+1), dtype=np.complex128)
    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    Ylm.real[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    Ylm.imag[LMIN:LMAX+1,mm] = -slm1[LMIN:LMAX+1,mm]
    # calculate Ylms summation for each theta band
    delta_M = np.einsum("lmh...,lm...->mh...", PLM, Ylm / 2.0)

    # output spatial field from FFT transformation
    s = np.zeros((phimax,thmax))
    # calculate fft for each theta band (over phis with axis=0)
    s[:-1,:] = 2.0*(phimax-1)*np.fft.ifft(delta_M,n=phimax-1,axis=0).real
    # complete sphere (values at 360 == values at 0)
    s[-1,:] = s[0,:]

    # return output data
    return s

def stokes_summation(clm1, slm1, lon, lat,
    LMIN=0, LMAX=60, MMAX=None, RAD=0, UNITS=0, LOVE=None, PLM=None):
    r"""
    Converts data from spherical harmonic coefficients to a spatial field
    :cite:p:`Wahr:1998hy`

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
    UNITS: int, default 0
        Output data units

            - ``1``: cm water equivalent thickness (cm w.e., g/cm\ :sup:`2`)
            - ``2``: mm geoid height
            - ``3``: mm elastic crustal deformation :cite:p:`Davis:2004il`
            - ``4``: microGal gravitational perturbation
            - ``5``: mbar equivalent surface pressure
            - ``6``: cm viscoelastic crustal uplift (GIA) :cite:p:`Wahr:2000ek`
            - list: custom degree-dependent unit conversion factor
    LMAX: int, default 0
        Upper bound of Spherical Harmonic Degrees
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    PLM: np.ndarray or NoneType, default None
        Fully-normalized associated Legendre polynomials

    Returns
    -------
    spatial: np.ndarray
        spatial field
    """
    # if LMAX is not specified, will use the size of the input harmonics
    if (LMAX == 0):
        LMAX = np.shape(clm1)[0]-1
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Gaussian Smoothing
    if (RAD != 0):
        wl = 2.0*np.pi*gauss_weights(RAD, LMAX)
    else:
        # else = 1
        wl = np.ones((LMAX+1))

    # Setting units factor for output
    # dfactor is the degree dependent coefficients
    factors = units(lmax=LMAX)
    if isinstance(UNITS, (list,np.ndarray)):
        # custom units
        dfactor = np.copy(UNITS)
    elif isinstance(UNITS, str):
        # named units
        dfactor = factors.harmonic(*LOVE).get(UNITS)
    elif isinstance(UNITS, int):
        # use named unit codes
        dfactor = factors.harmonic(*LOVE).get(units.bycode(UNITS))
    else:
        raise ValueError(f'Unknown units {UNITS}')

    # spherical harmonic order
    mm = np.arange(0,MMAX+1)# mmax+1 to include mmax
    # smooth harmonics and convert to output units
    clm = np.einsum("l,l,lm->lm", wl, dfactor, clm1[:, mm])
    slm = np.einsum("l,l,lm->lm", wl, dfactor, slm1[:, mm])

    # return the spatial field
    return harmonic_summation(clm, slm, lon, lat,
        LMIN=LMIN, LMAX=LMAX, MMAX=MMAX, PLM=PLM)
