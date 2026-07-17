#!/usr/bin/env python
"""
gen_harmonics.py
Written by Tyler Sutterley (07/2026)
Converts data from the spatial domain to spherical harmonic coefficients
Does not compute the solid Earth elastic response or convert units

CALLING SEQUENCE:
    Ylms = gen_harmonics(data, lon, lat, LMIN=0, LMAX=60)

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
        or fourier coefficients of Legendre polynomials
    METHOD: conversion method for calculating harmonics
        integration: for global grids
        fourier: for regional or global grids

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    fourier_legendre.py: Computes the Fourier coefficients of the associated
        Legendre functions
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCES:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)

UPDATE HISTORY:
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
        added custom weighting function for gridded data
    Updated 03/2023: improve typing for variables in docstrings
    Updated 01/2023: refactored associated legendre polynomials
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: merged integration and fourier harmonics programs
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Updated 04/2020: include degrees and orders in output dictionary
    Updated 10/2017: updated comments and cleaned up code
    Updated 08/2017: Using Holmes and Featherstone relation for Plms
    Updated 08/2015: changed from sys.exit to raise ValueError
    Updated 05/2015: updated output for MMAX != LMAX
    Written 05/2013
"""

import numpy as np
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.fourier_legendre import fourier_legendre


def gen_harmonics(data, lon, lat, **kwargs):
    """
    Converts data from the spatial domain to spherical harmonic coefficients

    Parameters
    ----------
    data: np.ndarray
        data magnitude
    lon: np.ndarray
        longitude array
    lat: np.ndarray
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    PLM: np.ndarray, default 0
        Fully normalized associated Legendre polynomials or
        Fourier coefficients of Legendre polynomials
    METHOD: str, default 'integration'
        Conversion method for calculating harmonics

            - ``'integration'``: for global grids
            - ``'fourier'``: for regional or global grids

    Returns
    -------
    clm: np.ndarray
        cosine spherical harmonic coefficients (4-pi normalized)
    slm: np.ndarray
        sine spherical harmonic coefficients (4-pi normalized)
    l: np.ndarray
        spherical harmonic degree to LMAX
    m: np.ndarray
        spherical harmonic order to MMAX
    """
    # set default keyword arguments
    kwargs.setdefault('LMAX', 60)
    kwargs.setdefault('MMAX', None)
    kwargs.setdefault('PLM', 0)
    kwargs.setdefault('METHOD', 'integration')
    # upper bound of spherical harmonic orders (default = LMAX)
    if kwargs['MMAX'] is None:
        kwargs['MMAX'] = np.copy(kwargs['LMAX'])
    # convert latitude and longitude to float if integers
    lon = lon.astype(np.float64)
    lat = lat.astype(np.float64)
    # reforming data to lonXlat if input latXlon
    sz = np.shape(data)
    dinput = np.transpose(data) if (sz[0] == len(lat)) else np.copy(data)
    # convert spatial field into spherical harmonics
    if kwargs['METHOD'].lower() == 'integration':
        Ylms = integration(dinput, lon, lat, **kwargs)
    elif kwargs['METHOD'].lower() == 'fourier':
        Ylms = fourier(dinput, lon, lat, **kwargs)
    # return the output spherical harmonics object
    return Ylms


def integration(
    data,
    lon,
    lat,
    LMAX=60,
    MMAX=None,
    WEIGHT=None,
    PLM=0,
    **kwargs,
):
    """
    Converts data from the spatial domain to spherical harmonic coefficients

    Parameters
    ----------
    data: float
        data magnitude
    lon: float
        longitude array
    lat: float
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    WEIGHT: np.ndarray or NoneType, default None
        Custom latitudinal weighting function for gridded data
    PLM: float, default 0
        input Legendre polynomials

    Returns
    -------
    clm: float
        cosine spherical harmonic coefficients
    slm: float
        sine spherical harmonic coefficients
    l: int
        spherical harmonic degree to LMAX
    m: int
        spherical harmonic order to MMAX
    """
    # calculate longitude and colatitude arrays in radians
    phi = np.radians(np.squeeze(lon))
    th = np.radians(90.0 - np.squeeze(lat))
    # reformatting longitudes to range 0:360 (if previously -180:180)
    phi = np.where(phi < 0, phi + 2.0 * np.pi, phi)
    # grid dimensions
    nlat = np.int64(len(th))

    # LMAX+1 as there are LMAX+1 elements between 0 and LMAX
    ll = np.arange(LMAX + 1)
    mm = np.arange(MMAX + 1)
    # Calculating cos/sin of phi arrays (output [m,phi])
    m_phi = np.exp(1j * np.einsum('m...,p...->mp...', mm, phi))

    # use an integration factor for gridded data or
    # calculate from sin(theta)*dtheta*dphi
    int_fact = np.zeros((nlat))
    if WEIGHT is not None:
        # Weighting function for integrating gridded data
        int_fact[:] = np.broadcast_to(np.atleast_1d(WEIGHT), nlat)
    else:
        # Multiplying sin(th) with differentials of theta and phi
        # to calculate the integration factor at each latitude
        dphi = np.abs(phi[1] - phi[0])
        dth = np.abs(th[1] - th[0])
        int_fact[:] = np.sin(th) * dphi * dth
    # normalizing coefficients
    coeff = 1.0 / (4.0 * np.pi)

    # Calculate polynomials using Holmes and Featherstone (2002) relation
    if np.ndim(PLM) == 0:
        PLM, dplm = plm_holmes(LMAX, np.cos(th))
    # Multiply plms by integration factors [sin(theta)*dtheta*dphi]
    # truncate plms to maximum spherical harmonic order if MMAX < LMAX
    plm = np.einsum(
        'lmh...,h...->lmh...', PLM[: LMAX + 1, : MMAX + 1, :], int_fact
    )
    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX + 1, MMAX + 1))
    Ylms.slm = np.zeros((LMAX + 1, MMAX + 1))
    # Multiplying gridded data with sin/cos of m#phis (output [m,theta])
    # This will sum through all phis in the dot product
    d = np.einsum('mp...,ph...->mh...', m_phi, data)
    # Summing product of plms and data over all latitudes
    ylm = np.einsum('lmh...,mh...->lm...', plm, d)
    # convert to output normalization (4-pi normalized harmonics)
    # truncate to MMAX if specified (if l > MMAX)
    Ylms.clm = coeff * ylm.real[: LMAX + 1, : MMAX + 1]
    Ylms.slm = coeff * ylm.imag[: LMAX + 1, : MMAX + 1]
    # return the output spherical harmonics object
    return Ylms


def fourier(
    data,
    lon,
    lat,
    LMAX=60,
    MMAX=None,
    PLM=0,
    **kwargs,
):
    """
    Computes the spherical harmonic coefficients of a spatial field

    Parameters
    ----------
    data: float
        data magnitude
    lon: float
        longitude array
    lat: float
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    PLM: float, default 0
        input Fourier coefficients of Legendre polynomials

    Returns
    -------
    clm: float
        cosine spherical harmonic coefficients
    slm: float
        sine spherical harmonic coefficients
    l: int
        spherical harmonic degree to LMAX
    m: int
        spherical harmonic order to MMAX
    """

    # dimensions of the longitude and latitude arrays
    nlon = np.int64(len(lon))
    nlat = np.int64(len(lat))
    # calculate longitude and colatitude arrays in radians
    phi = np.radians(np.squeeze(lon))
    th = np.radians(90.0 - np.squeeze(lat))
    # reformatting longitudes to range 0:360 (if previously -180:180)
    phi = np.where(phi < 0, phi + 2.0 * np.pi, phi)
    # grid step in radians
    dphi = np.abs(phi[1] - phi[0])
    dth = np.abs(th[1] - th[0])

    # MMAX+1 to include MMAX
    mm = np.arange(MMAX + 1)
    # Calculate cos and sin coefficients of signal
    m_phi = np.exp(1j * np.einsum('m...,p...->mp...', mm, phi))
    d = np.einsum('mp...,ph...->mh...', m_phi, data)
    # normalize coefficients
    d[0, :] *= 1.0 / nlon
    d[1:, :] *= 2.0 / nlon

    # Calculate cos and sin coefficients of theta component
    # Because the function is defined on (0,pi)
    # it can be expanded in just cosine terms.
    # this routine assumes that 0 and pi are not included
    f = np.zeros((MMAX + 1, MMAX + 1), dtype=np.complex128)
    m_even = slice(0, MMAX + 1, 2)
    m_odd = slice(1, MMAX, 2)

    if np.isclose([th[0], th[nlat - 1]], [0.0, np.pi]).all():
        # global case (includes poles)
        # non-endpoints
        k_th = np.exp(1j * np.einsum('h...,k...->kh...', th[1 : nlat - 1], mm))
        f[m_even, :] = 2.0 * np.einsum(
            'mh...,kh...->mk', d[m_even, 1 : nlat - 1], k_th.real
        )
        f[m_odd, :] = 2.0 * np.einsum(
            'mh...,kh...->mk', d[m_odd, 1 : nlat - 1], k_th.imag
        )
        # endpoints
        k_th = np.exp(1j * mm * th[0])
        f[m_even, :] += np.einsum('m...,k...->mk', d[m_even, 0], k_th)
        f[m_odd, :] += np.einsum('m...,k...->mk', d[m_odd, 0], k_th)
        k_th = np.exp(1j * mm * th[nlat - 1])
        f[m_even, :] += np.einsum('m...,k...->mk', d[m_even, nlat - 1], k_th)
        f[m_odd, :] += np.einsum('m...,k...->mk', d[m_odd, nlat - 1], k_th)
    elif not np.isclose([th[0], th[nlat - 1]], [0.0, np.pi]).any():
        k_th = np.exp(1j * np.einsum('h...,k...->kh...', th, mm))
        f[m_even, :] = 2.0 * np.einsum(
            'mh...,kh...->mk', d[m_even, :], k_th.real
        )
        f[m_odd, :] = 2.0 * np.einsum('mh...,kh...->mk', d[m_odd, :], k_th.imag)
    else:
        raise ValueError('Latitude coordinates incompatible')

    # Normalize theta fourier coefficients
    f[:, 0] *= 1.0 / (2.0 * nlat)
    f[:, 1 : MMAX + 1] *= 1.0 / nlat
    # Correct normalization for the incomplete coverage of the sphere
    f[:] *= nlon * dphi / (2.0 * np.pi) * nlat * dth / np.pi

    # Calculate cos and sin coefficients of Legendre functions
    # Expand m = even terms in a cosine series
    # Expand m = odd terms in a sine series
    # Both are stride 2
    if np.ndim(PLM) == 0:
        Almk = fourier_legendre(LMAX, MMAX)
    else:
        # use precomputed alms to improve computational speed
        Almk = PLM

    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX + 1, MMAX + 1))
    Ylms.slm = np.zeros((LMAX + 1, MMAX + 1))

    # calculate spherical harmonics for m == even terms
    # even l terms (l even, m even, k even)
    l_even = slice(0, LMAX + 1, 2)
    n_even = np.arange(m_even.start, m_even.stop, m_even.step)
    k_even = np.zeros((len(n_even), len(n_even)))
    for k in range(0, MMAX + 2, 2):
        k_even[:, k // 2] = 0.5 * (
            1.0 / (1.0 - n_even - k)
            + 1.0 / (1.0 + n_even - k)
            + 1.0 / (1.0 - n_even + k)
            + 1.0 / (1.0 + n_even + k)
        )
    # calculate summation over coefficients
    Aeven = np.einsum(
        'lmk...,kn...->lmn...', Almk[l_even, m_even, m_even], k_even
    )
    Yeven = np.einsum('lmn...,mn...->lm...', Aeven, f[m_even, m_even])
    Ylms.clm[l_even, m_even] = Yeven.real
    Ylms.slm[l_even, m_even] = Yeven.imag

    # odd l terms (l odd, m even, k odd)
    l_odd = slice(1, LMAX, 2)
    n_odd = np.arange(m_odd.start, m_odd.stop, m_odd.step)
    k_odd = np.zeros((len(n_odd), len(n_odd)))
    for k in range(1, MMAX + 1, 2):
        k_odd[:, (k - 1) // 2] = 0.5 * (
            1.0 / (1.0 - n_odd - k)
            + 1.0 / (1.0 + n_odd - k)
            + 1.0 / (1.0 - n_odd + k)
            + 1.0 / (1.0 + n_odd + k)
        )
    # calculate summation over coefficients
    Aodd = np.einsum('lmk...,kn...->lmn...', Almk[l_odd, m_even, m_odd], k_odd)
    Yodd = np.einsum('lmn...,mn...->lm...', Aodd, f[m_even, m_odd])
    Ylms.clm[l_odd, m_even] = Yodd.real
    Ylms.slm[l_odd, m_even] = Yodd.imag

    # calculate spherical harmonics for m == odd terms
    # even l terms (l even, m odd, k even)
    l_even = slice(2, LMAX + 1, 2)  # do not in include l=0
    n_even = np.arange(m_even.start, m_even.stop, m_even.step)
    k_even = np.zeros((len(n_even), len(n_even)))
    for k in range(0, MMAX + 2, 2):
        k_even[:, k // 2] = 0.5 * (
            -1.0 / (1.0 - n_even - k)
            + 1.0 / (1.0 + n_even - k)
            + 1.0 / (1.0 - n_even + k)
            - 1.0 / (1.0 + n_even + k)
        )
    Aeven = np.einsum(
        'lmk...,kn...->lmn...', Almk[l_even, m_odd, m_even], k_even
    )
    Yeven = np.einsum('lmn...,mn...->lm...', Aeven, f[m_odd, m_even])
    Ylms.clm[l_even, m_odd] = Yeven.real
    Ylms.slm[l_even, m_odd] = Yeven.imag

    # odd l terms (l odd, m odd, k odd)
    l_odd = slice(1, LMAX, 2)
    n_odd = np.arange(m_odd.start, m_odd.stop, m_odd.step)
    k_odd = np.zeros((len(n_odd), len(n_odd)))
    for k in range(1, MMAX + 1, 2):
        k_odd[:, (k - 1) // 2] = 0.5 * (
            -1.0 / (1.0 - n_odd - k)
            + 1.0 / (1.0 + n_odd - k)
            + 1.0 / (1.0 - n_odd + k)
            - 1.0 / (1.0 + n_odd + k)
        )
    # calculate summation over coefficients
    Aodd = np.einsum('lmk...,kn...->lmn...', Almk[l_odd, m_odd, m_odd], k_odd)
    Yodd = np.einsum('lmn...,mn...->lm...', Aodd, f[m_odd, m_odd])
    Ylms.clm[l_odd, m_odd] = Yodd.real
    Ylms.slm[l_odd, m_odd] = Yodd.imag

    # Divide by Plm normalization
    Ylms.clm[:, 0] /= 2.0
    Ylms.slm[:, 0] /= 2.0
    Ylms.clm[:, 1 : MMAX + 1] /= 4.0
    Ylms.slm[:, 1 : MMAX + 1] /= 4.0

    # return the output spherical harmonics object
    return Ylms
