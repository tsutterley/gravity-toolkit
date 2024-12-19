#!/usr/bin/env python
u"""
harmonic_summation.py
Written by Tyler Sutterley (03/2023)

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

    # Longitude in radians
    phi = (np.squeeze(lon)*np.pi/180.0)[np.newaxis,:]
    # colatitude in radians
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    thmax = len(th)

    # if plms are not pre-computed: calculate Legendre polynomials
    if PLM is None:
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # Truncating harmonics to degree and order LMAX
    # removing coefficients below LMIN and above MMAX
    mm = np.arange(0, MMAX+1)
    clm = np.zeros((LMAX+1, MMAX+1))
    slm = np.zeros((LMAX+1, MMAX+1))
    clm[LMIN:LMAX+1,mm] = clm1[LMIN:LMAX+1,mm]
    slm[LMIN:LMAX+1,mm] = slm1[LMIN:LMAX+1,mm]
    # Calculate fourier coefficients from legendre coefficients
    d_cos = np.zeros((MMAX+1,thmax))# [m,th]
    d_sin = np.zeros((MMAX+1,thmax))# [m,th]
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

    # if plms are not pre-computed: calculate Legendre polynomials
    if PLM is None:
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # combined Ylms and Fourier coefficients (complex)
    Ylms = np.zeros((LMAX+1, MMAX+1),dtype=np.complex128)
    delta_M = np.zeros((MMAX+1,thmax),dtype=np.complex128)# [m,th]
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

    # truncate to degree and order
    mm = np.arange(0, MMAX+1)
    # smooth harmonics and convert to output units
    clm = np.zeros((LMAX+1, MMAX+1))
    slm = np.zeros((LMAX+1, MMAX+1))
    for l in range(0, LMAX+1):# LMAX+1 to include LMAX
        clm[l,:] = wl[l]*dfactor[l]*clm1[l,mm]
        slm[l,:] = wl[l]*dfactor[l]*slm1[l,mm]

    # return the spatial field
    return harmonic_summation(clm, slm, lon, lat,
        LMIN=LMIN, LMAX=LMAX, MMAX=MMAX, PLM=PLM)
