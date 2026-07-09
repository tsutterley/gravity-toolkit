#!/usr/bin/env python
u"""
gen_stokes.py
Written by Tyler Sutterley (07/2026)

Converts data from the spatial domain to spherical harmonic coefficients

CALLING SEQUENCE:
    Ylms = gen_stokes(data, lon, lat, UNITS=1, LMIN=0, LMAX=60, LOVE=(hl,kl,ll))

INPUTS:
    data: data matrix
    lon: longitude array
    lat: latitude array

OUTPUTS:
    Ylms: harmonics object
        clm: fully-normalized cosine spherical harmonic coefficients
        slm: fully-normalied sine spherical harmonic coefficients
        l: spherical harmonic degree to LMAX
        m: spherical harmonic order to MMAX

OPTIONS:
    LMIN: Lower bound of Spherical Harmonic Degrees (default = 0)
    LMAX: Upper bound of Spherical Harmonic Degrees (default = 60)
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    UNITS: input data units
        1: cm of water thickness (default)
        2: Gigatonnes of mass
        3: kg/m^2
        list: custom degree-dependent unit conversion factor
    PLM: input Legendre polynomials
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: computes fully-normalized associated Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

UPDATE HISTORY:
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
    Updated 06/2025: copy latitude and longitude as float64 for numpy 2.0 stability
    Updated 04/2023: allow love numbers to be None for custom units case
    Updated 03/2023: improve typing for variables in docstrings
    Updated 02/2023: set custom units as top option in if/else statements
    Updated 01/2023: refactored associated legendre polynomials
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: added UNITS list option for converting from custom units
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Updated 04/2020: reading load love numbers outside of this function
        using the units class for converting to normalized spherical harmonics
        include degrees and orders in output dictionary for harmonics class
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2018: use copies of longitude and latitude to not modify inputs
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 05/2015: added parameter MMAX for LMAX != MMAX
    Updated 06/2014: changed message to sys.exit
    Updated 02/2014: minor update to if statements
    Updated 05/2013: added option to precompute plms
    Updated 05/2013: added linear interpolation of love numbers for LMAX > 696
    Updated 05/2013: transpose data to (LON,LAT) if originally (LAT,LON)
    Updated 04/2012: added lmin/lmax options
    Updated 02/2012: added DLON and DLAT options for different degree spacing
        revised structure of mathematics to improve computational efficiency
    Written 09/2011
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes

def gen_stokes(data, lon, lat, LMIN=0, LMAX=60, MMAX=None, UNITS=1,
    PLM=None, LOVE=None):
    r"""
    Converts data from the spatial domain to spherical harmonic
    coefficients :cite:p:`Wahr:1998hy`

    Parameters
    ----------
    data: np.ndarray
        data matrix
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
    UNITS: int, default 1
        Input data units

            - ``1``: cm water equivalent thickness (cm w.e., g/cm\ :sup:`2`)
            - ``2``: gigatonnes of mass (Gt)
            - ``3``:  mm water equivalent thickness (mm w.e., kg/m\ :sup:`2`)
            - list: custom degree-dependent unit conversion factor
    PLM: np.ndarray or NoneType, default None
        Input Legendre polynomials
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

    # converting LMIN and LMAX to integer
    LMIN = np.int64(LMIN)
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX

    # grid dimensions
    nlat = np.int64(len(lat))
    # Longitude in radians
    phi = np.radians(np.squeeze(lon.copy()))
    # reformatting longitudes to range 0:360 (if previously -180:180)
    phi = np.where(phi < 0, phi + 2.0*np.pi, phi)
    # colatitude in radians
    th = np.radians(90.0 - np.squeeze(lat.copy()))
    # grid step in radians
    dphi = np.abs(phi[1] - phi[0])
    dth = np.abs(th[1] - th[0])

    # reforming data to lonXlat if input latXlon
    sz = np.shape(data)
    data = data.T if (sz[0] == nlat) else np.copy(data)

    # extract degree dependent factor for specific units
    # calculate integration factors for theta and phi
    # Multiplying sin(th) with differentials of theta and phi
    # to calculate the integration factor at each latitude
    factors = gravity_toolkit.units(lmax=LMAX)
    int_fact = np.zeros((nlat))
    if isinstance(UNITS, (list, np.ndarray)):
        # custom units
        dfactor = np.copy(UNITS)
        int_fact[:] = np.sin(th)*dphi*dth
    elif (UNITS == 1):
        # Default Parameter: Input in cm w.e. (g/cm^2)
        dfactor = factors.spatial(*LOVE).cmwe
        int_fact[:] = np.sin(th)*dphi*dth
    elif (UNITS == 2):
        # Input in gigatonnes (Gt)
        dfactor = factors.spatial(*LOVE).cmwe
        # rad_e: Average Radius of the Earth [cm]
        int_fact[:] = 1e15/(factors.rad_e**2)
    elif (UNITS == 3):
        # Input in kg/m^2 (mm w.e.)
        dfactor = factors.spatial(*LOVE).mmwe
        int_fact[:] = np.sin(th)*dphi*dth
    else:
        raise ValueError(f'Unknown units {UNITS}')

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    mm = np.arange(MMAX+1)
    m_phi = np.exp(1j * np.einsum("m...,p...->mp...", mm, phi))

    # Calculating fully-normalized Legendre Polynomials
    # Output is plm[l,m,th]
    plm = np.zeros((LMAX+1, MMAX+1, nlat))
    # added option to precompute plms to improve computational speed
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # truncate legendre polynomials to degree and order
    plm = np.einsum("lmh...,h...->lmh...", PLM[:LMAX+1,:MMAX+1,:], int_fact)

    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    # Multiplying gridded data with sin/cos of m#phis
    # This will sum through all phis in the dot product
    # output [m,theta]
    d = np.einsum("mp...,ph...->mh...", m_phi, data)
    # Summing product of plms and data over all latitudes
    ylm = np.einsum("lmh...,mh...->lm...", plm, d)
    # Multiplying by factors to convert to fully normalized coefficients
    Ylms.clm = np.einsum("l...,lm...->lm...", dfactor, ylm.real)
    Ylms.slm = np.einsum("l...,lm...->lm...", dfactor, ylm.imag)

    # return the output spherical harmonics object
    return Ylms