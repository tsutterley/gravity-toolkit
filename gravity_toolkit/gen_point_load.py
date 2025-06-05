#!/usr/bin/env python
u"""
gen_point_load.py
Written by Tyler Sutterley (04/2023)
Calculates gravitational spherical harmonic coefficients for point masses

CALLING SEQUENCE:
    Ylms = gen_point_load(data, lon, lat, LMAX=LMAX)

INPUTS:
    data: data magnitude
    lon: longitude of points
    lat: latitude of points

OUTPUTS:
    clm: cosine spherical harmonic coefficients (geodesy normalization)
    slm: sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    UNITS: input data units
        1: grams of mass (default)
        2: gigatonnes of mass
        list: custom degree-dependent unit conversion factor
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

PROGRAM DEPENDENCIES:
    legendre.py: Computes associated Legendre polynomials for degree l
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCES:
    I. M. Longman, Journal of Geophysical Research, 67(2), 1962
        https://doi.org/10.1029/JZ067i002p00845
    W. E. Farrell, Reviews of Geophysics and Space Physics, 10(3), 1972
        https://doi.org/10.1029/RG010i003p00761
    H. N. Pollack, Journal of Geophysical Research, 78(11), 1973
        https://doi.org/10.1029/JB078i011p01760

UPDATE HISTORY:
    Updated 04/2023: allow love numbers to be None for custom units case
    Updated 03/2023: improve typing for variables in docstrings
    Updated 02/2023: set custom units as top option in if/else statements
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: added UNITS list option for converting from custom units
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Written 05/2020
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.legendre import legendre

def gen_point_load(data, lon, lat, LMAX=60, MMAX=None, UNITS=1, LOVE=None):
    """
    Calculates spherical harmonic coefficients for point masses
    :cite:p:`Longman:1962ev,Farrell:1972cm,Pollack:1973gi`

    Parameters
    ----------
    data: np.ndarray
        data magnitude
    lon: np.ndarray
        longitude of points
    lat: np.ndarray
        latitude of points
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    UNITS: int, default 1
        Input data units

            - ``1``: grams of mass (g)
            - ``2``: gigatonnes of mass (Gt)
            - list: custom degree-dependent unit conversion factor
    LOVE: tuple or NoneType, default None
        Input load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

    Returns
    -------
    clm: np.ndarray
        cosine spherical harmonic coefficients
    slm: np.ndarray
        sine spherical harmonic coefficients
    l: np.ndarray
        spherical harmonic degree to LMAX
    m: np.ndarray
        spherical harmonic order to MMAX
    """

    # upper bound of spherical harmonic orders (default == LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # number of input data points
    npts = len(data.flatten())
    # convert output longitude and latitude into radians
    phi = np.pi*lon.flatten()/180.0
    theta = np.pi*(90.0 - lat.flatten())/180.0

    # extract degree dependent factor for specific units
    factors = gravity_toolkit.units(lmax=LMAX)
    int_fact = np.zeros((npts))
    if isinstance(UNITS, (list, np.ndarray)):
        # custom units
        dfactor = np.copy(UNITS)
        int_fact[:] = 1.0
    elif (UNITS == 1):
        # Default Parameter: Input in grams (g)
        dfactor = factors.spatial(*LOVE).cmwe/(factors.rad_e**2)
        int_fact[:] = 1.0
    elif (UNITS == 2):
        # Input in gigatonnes (Gt)
        dfactor = factors.spatial(*LOVE).cmwe/(factors.rad_e**2)
        int_fact[:] = 1e15
    else:
        raise ValueError(f'Unknown units {UNITS}')
    # flattened form of data converted to units
    D = int_fact*data.flatten()

    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1, MMAX+1))
    Ylms.slm = np.zeros((LMAX+1, MMAX+1))
    # for each degree l
    for l in range(LMAX+1):
        m1 = np.min([l,MMAX]) + 1
        SPH = spherical_harmonic_matrix(l, D, phi, theta, dfactor[l])
        # truncate to spherical harmonic order and save to output
        Ylms.clm[l,:m1] = SPH.real[:m1]
        Ylms.slm[l,:m1] = SPH.imag[:m1]
    # return the output spherical harmonics object
    return Ylms

# calculate spherical harmonics of degree l evaluated at (theta,phi)
def spherical_harmonic_matrix(l, data, phi, theta, coeff):
    """
    Calculates the spherical harmonics for a particular degree
    evaluated from data at coordinates

    Parameters
    ----------
    l: int
        spherical harmonic degree
    data: np.ndarray
        data magnitude in grams
    phi: np.ndarray
        longitude of points in radians
    theta: np.ndarray
        colatitude of points in radians
    coeff: np.ndarray
        degree-dependent factor for converting units

    Returns
    -------
    Ylms: np.ndarray
        spherical harmonic coefficients in Eulerian form
    """
    # calculate normalized legendre polynomials (points, order)
    Pl = legendre(l, np.cos(theta), NORMALIZE=True).T
    # spherical harmonic orders up to degree l
    m = np.arange(0, l+1)
    # calculate Euler's of spherical harmonic order multiplied by azimuth phi
    mphi = np.exp(1j*np.dot(np.squeeze(phi)[:,np.newaxis], m[np.newaxis,:]))
    # reshape data to order
    D = np.kron(np.ones((1, l+1)), data[:,np.newaxis])
    # calculate spherical harmonics and multiply by coefficients and data
    Ylms = coeff*D*Pl*mphi
    # calculate the sum over all points and return harmonics for degree l
    return np.sum(Ylms, axis=0)
