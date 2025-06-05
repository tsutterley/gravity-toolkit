#!/usr/bin/env python
u"""
gen_disc_load.py
Written by Tyler Sutterley (06/2023)
Calculates gravitational spherical harmonic coefficients for a uniform disc load

CALLING SEQUENCE:
    Ylms = gen_disc_load(data, lon, lat, area, LMAX=60, MMAX=None)

INPUTS:
    data: data magnitude (Gt)
    lon: longitude of disc center
    lat: latitude of disc center
    area: area of disc (km^2)

OUTPUTS:
    clm: cosine spherical harmonic coefficients (geodesy normalization)
    slm: sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    UNITS: input data units
        1: cm of water thickness
        2: Gigatonnes of mass
        3: kg/m^2
        list: custom unit conversion factor
    PLM: input Legendre polynomials
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    legendre_polynomials.py: Computes fully normalized Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCES:
    Holmes and Featherstone, Journal of Geodesy, 76, 279-299, 2002
        https://doi.org/10.1007/s00190-002-0216-2
    I. M. Longman, Journal of Geophysical Research, 67(2), 1962
        https://doi.org/10.1029/JZ067i002p00845
    W. E. Farrell, Reviews of Geophysics and Space Physics, 10(3), 1972
        https://doi.org/10.1029/RG010i003p00761
    H. N. Pollack, Journal of Geophysical Research, 78(11), 1973
        https://doi.org/10.1029/JB078i011p01760
    T. Jacob et al., Journal of Geodesy, 86, 337-358, 2012
        https://doi.org/10.1007/s00190-011-0522-7

UPDATE HISTORY:
    Updated 06/2023: modified custom units case to not convert to cmwe
    Updated 03/2023: simplified unit degree factors using units class
        improve typing for variables in docstrings
    Updated 02/2023: set custom units as top option in if/else statements
    Updated 01/2023: refactored associated legendre polynomials
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: added UNITS option for converting from different inputs
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Updated 05/2020: vectorize calculation over degrees to improve compute time
        added option to precompute plms for disc load centers
    Updated 04/2020: reading load love numbers outside of this function
        include degrees and orders in output dictionary for harmonics class
        use units class for Earth parameters
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Updated 08/2017: Using Holmes and Featherstone relation for Plms
    Written 09/2016
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.legendre_polynomials import legendre_polynomials

def gen_disc_load(data, lon, lat, area, LMAX=60, MMAX=None, UNITS=2,
    PLM=None, LOVE=None):
    r"""
    Calculates spherical harmonic coefficients for a uniform disc load
    :cite:p:`Holmes:2002ff,Longman:1962ev,Farrell:1972cm,Pollack:1973gi,Jacob:2012eo`

    Parameters
    ----------
    data: np.ndarray
        data magnitude (Gt)
    lon: np.ndarray
        longitude of disc center
    lat: np.ndarray
        latitude of disc center
    area: float
        area of disc (km\ :sup:`2`)
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    UNITS: int, default 2
        Input data units

            - ``1``: cm water equivalent thickness (cm w.e., g/cm\ :sup:`2`)
            - ``2``: gigatonnes of mass (Gt)
            - ``3``:  mm water equivalent thickness (mm w.e., kg/m\ :sup:`2`)
            - list: custom unit conversion factor
    PLM: np.ndarray or NoneType, default None
        Legendre polynomials for ``cos(theta)`` (disc center)
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

    Returns
    -------
    clm: np.ndarray
        cosine spherical harmonic coefficients (geodesy normalization)
    slm: np.ndarray
        sine spherical harmonic coefficients (geodesy normalization)
    l: np.ndarray
        spherical harmonic degree to LMAX
    m: np.ndarray
        spherical harmonic order to MMAX
    """

    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # convert lon and lat to radians
    phi = lon*np.pi/180.0# Longitude in radians
    th = (90.0 - lat)*np.pi/180.0# Colatitude in radians

    # Earth Parameters
    factors = gravity_toolkit.units(lmax=LMAX)

    # convert input area into cm^2 and then divide by area of a half sphere
    # alpha will be 1 - the ratio of the input area with the half sphere
    alpha = (1.0 - 1e10*area/(2.0*np.pi*factors.rad_e**2))

    # Calculate factor to convert from input units into g/cm^2
    if isinstance(UNITS, (list, np.ndarray)):
        # custom units
        unit_conv = 1.0
        dfactor = np.copy(UNITS)
    elif (UNITS == 1):
        # Input data is in cm water equivalent (cmwe)
        unit_conv = 1.0
        # degree dependent factors to convert from coefficients
        # of mass into normalized geoid coefficients
        dfactor = 4.0*np.pi*factors.spatial(*LOVE).cmwe/(1.0 + 2.0*factors.l)
    elif (UNITS == 2):
        # Input data is in gigatonnes (Gt)
        # 1e15 converts from Gt to grams, 1e10 converts from km^2 to cm^2
        unit_conv = 1e15/(1e10*area)
        # degree dependent factors to convert from coefficients
        # of mass into normalized geoid coefficients
        dfactor = 4.0*np.pi*factors.spatial(*LOVE).cmwe/(1.0 + 2.0*factors.l)
    elif (UNITS == 3):
        # Input data is in kg/m^2
        # 1 kg = 1000 g
        # 1 m^2 = 100*100 cm^2 = 1e4 cm^2
        unit_conv = 0.1
        # degree dependent factors to convert from coefficients
        # of mass into normalized geoid coefficients
        dfactor = 4.0*np.pi*factors.spatial(*LOVE).cmwe/(1.0 + 2.0*factors.l)
    else:
        raise ValueError(f'Unknown units {UNITS}')

    # Calculating plms of the disc
    # allocating for constructed array
    pl_alpha = np.zeros((LMAX+1))
    # l=0 is a special case (P(-1) = 1, P(1) = cos(alpha))
    pl_alpha[0] = (1.0 - alpha)/2.0
    # for all other degrees: calculate the legendre polynomials up to LMAX+1
    pl_matrix,_ = legendre_polynomials(LMAX+1,alpha)
    for l in range(1, LMAX+1):# LMAX+1 to include LMAX
        # from Longman (1962) and Jacob et al (2012)
        # unnormalizing Legendre polynomials
        # sqrt(2*l - 1) == sqrt(2*(l-1) + 1)
        # sqrt(2*l + 3) == sqrt(2*(l+1) + 1)
        pl_lower = pl_matrix[l-1]/np.sqrt(2.0*l-1.0)
        pl_upper = pl_matrix[l+1]/np.sqrt(2.0*l+3.0)
        pl_alpha[l] = (pl_lower - pl_upper)/2.0

    # Calculate Legendre Polynomials using Holmes and Featherstone relation
    # this would be the plm for the center of the disc load
    # used to rotate the disc load to point lat/lon
    if PLM is None:
        plmout,_ = plm_holmes(LMAX, np.cos(th))
        # truncate precomputed plms to order
        plmout = np.squeeze(plmout[:,:MMAX+1,:])
    else:
        # truncate precomputed plms to degree and order
        plmout = PLM[:LMAX+1,:MMAX+1]

    # calculate array of m values ranging from 0 to MMAX (harmonic orders)
    # MMAX+1 as there are MMAX+1 elements between 0 and MMAX
    m = np.arange(MMAX+1)
    # Multiplying by the units conversion factor (unit_conv) to
    # convert from the input units into cmwe
    # Multiplying point mass data (converted to cmwe) with sin/cos of m*phis
    # data normally is 1 for a uniform 1cm water equivalent layer
    # but can be a mass point if reconstructing a spherical harmonic field
    # NOTE: NOT a matrix multiplication as data (and phi) is a single point
    dcos = unit_conv*data*np.cos(m*phi)
    dsin = unit_conv*data*np.sin(m*phi)

    # Multiplying by plm_alpha (F_l from Jacob 2012)
    plm = np.zeros((LMAX+1, MMAX+1))
    # Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1, MMAX+1))
    yslm = np.zeros((LMAX+1, MMAX+1))
    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1, MMAX+1))
    Ylms.slm = np.zeros((LMAX+1, MMAX+1))
    for m in range(0,MMAX+1):# MMAX+1 to include MMAX
        l = np.arange(m,LMAX+1)# LMAX+1 to include LMAX
        # rotate disc load to be centered at lat/lon
        plm[l,m] = plmout[l,m]*pl_alpha[l]
        # multiplying clm by cos(m*phi) and slm by sin(m*phi)
        # to get a field of spherical harmonics
        yclm[l,m] = plm[l,m]*dcos[m]
        yslm[l,m] = plm[l,m]*dsin[m]
        # multiplying by factors to convert to geoid coefficients
        Ylms.clm[l,m] = dfactor[l]*yclm[l,m]
        Ylms.slm[l,m] = dfactor[l]*yslm[l,m]

    # return the output spherical harmonics object
    return Ylms
