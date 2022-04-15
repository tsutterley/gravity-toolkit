#!/usr/bin/env python
u"""
clenshaw_summation.py
Written by Tyler Sutterley (04/2022)
Calculates the spatial field for a series of spherical harmonics for a
    sequence of ungridded points

CALLING SEQUENCE:
    spatial = clenshaw_summation(clm, slm, lon, lat, UNITS=1,
        LMAX=60, LOVE=(hl,kl,ll))

INPUTS:
    clm: cosine spherical harmonic coefficients
    slm: sine spherical harmonic coefficients
    lon: longitude of points
    lat: latitude of points

OPTIONS:
    RAD: Gaussian smoothing radius (km)
    UNITS: output data units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: mbar equivalent surface pressure
        6: cm of viscoelastic crustal uplift (GIA) [See Wahr 1995 or Wahr 2000]
        list: custom degree-dependent unit conversion factor
    LMAX: Upper bound of Spherical Harmonic Degrees
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)
    ASTYPE: floating point precision for calculating Clenshaw summation
    SCALE: scaling factor to prevent underflow in Clenshaw summation

OUTPUTS:
    spatial: spatial field for lon/lat

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    units.py: class for converting spherical harmonic data to specific units

REFERENCES:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        https://doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: added UNITS list option for converting to custom units
    Updated 09/2021: fix passing SCALE keyword argument to clenshaw_s_m
    Updated 06/2021: output equivalent pressure in pascals
    Updated 08/2020: parameterize float precision to improve computational time
    Updated 07/2020: added function docstrings
    Updated 04/2020: reading load love numbers outside of this function
        using the units class for converting normalized spherical harmonics
    Updated 03/2018: added option for output in equivalent pressure (UNITS=5)
        simplified love number extrapolation if LMAX is greater than 696
    Written 08/2017
"""
import numpy as np
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.units import units

def clenshaw_summation(clm, slm, lon, lat, RAD=0, UNITS=0, LMAX=0, LOVE=None,
    ASTYPE=np.longdouble, SCALE=1e-280):
    """
    Calculates the spatial field for a series of spherical harmonics for a
    sequence of ungridded points

    Parameters
    ----------
    clm: float
        cosine spherical harmonic coefficients
    slm: float
        sine spherical harmonic coefficients
    lon: float
        longitude of points
    lat: float
        latitude of points
    RAD: float, default 0.0
        Gaussian smoothing radius (km)
    UNITS: int, default 0
        Output data units

            - ``1``: cm water equivalent thickness (cm w.e., g/cm\ :sup:`2`)
            - ``2``: mm geoid height
            - ``3``: mm elastic crustal deformation
            - ``4``: microGal gravitational perturbation
            - ``5``: mbar equivalent surface pressure
            - ``6``: cm viscoelastic crustal uplift (GIA)
            - list: custom degree-dependent unit conversion factor
    LMAX: int, default 0
        Upper bound of Spherical Harmonic Degrees
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    ASTYPE: obj, default np.longdouble
        floating point precision for calculating Clenshaw summation
    SCALE: float, default 1e-280
        scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    spatial: float
        calculated spatial field for latitude and longitude

    References
    ----------
    .. [Davis2004] J. L. Davis et al.,
        "Climate‐driven deformation of the solid Earth from GRACE and GPS",
        *Geophysical Research Letters*, 31(L24605), (2004).
        `doi: 10.1029/2004GL021435 <https://doi.org/10.1029/2004GL021435>`_

    .. [Holmes2002] S. A. Holmes and W. E. Featherstone,
        "A unified approach to the Clenshaw summation and the recursive
        computation of very high degree and order normalised associated
        Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002).
        `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_

    .. [Tscherning1982] C. C. Tscherning and K. Poder,
        "Some Geodetic Applications of Clenshaw Summation",
        *Bollettino di Geodesia e Scienze*, 4, 349--375, (1982).

    .. [Wahr2000] J. Wahr, D. Wingham, and C. Bentley,
        "A method of combining ICESat and GRACE satellite data to constrain
        Antarctic mass balance", *Journal of Geophysical Research: Solid Earth*,
        105(B7), 16279--16294, (2000).
        `doi: 10.1029/2000JB900113 <https://doi.org/10.1029/2000JB900113>`_
    """

    #-- check if lat and lon are the same size
    if (len(lat) != len(lon)):
        raise ValueError('Incompatable vector dimensions (lon, lat)')

    #-- calculate colatitude and longitude in radians
    th = (90.0 - lat)*np.pi/180.0
    phi = np.squeeze(lon*np.pi/180.0)
    #-- calculate cos and sin of colatitudes
    t = np.cos(th)
    u = np.sin(th)

    #-- dimensions of theta and phi
    npts = len(th)

    #-- Gaussian Smoothing
    if (RAD != 0):
        wl = 2.0*np.pi*gauss_weights(RAD,LMAX)
    else:
        #-- else = 1
        wl = np.ones((LMAX+1))

    #-- Setting units factor for output
    #-- extract arrays of kl, hl, and ll Love Numbers
    factors = units(lmax=LMAX).harmonic(*LOVE)
    #-- dfactor computes the degree dependent coefficients
    if (UNITS == 0):
        #-- 0: keep original scale
        dfactor = factors.norm
    elif (UNITS == 1):
        #-- 1: cmH2O, centimeters water equivalent
        dfactor = factors.cmwe
    elif (UNITS == 2):
        #-- 2: mmGH, mm geoid height
        dfactor = factors.mmGH
    elif (UNITS == 3):
        #-- 3: mmCU, mm elastic crustal deformation
        dfactor = factors.mmCU
    elif (UNITS == 4):
        #-- 4: micGal, microGal gravity perturbations
        dfactor = factors.microGal
    elif (UNITS == 5):
        #-- 5: mbar, equivalent surface pressure
        dfactor = factors.mbar
    elif (UNITS == 6):
        #-- 6: cmVCU, cm viscoelastic  crustal uplift (GIA ONLY)
        dfactor = factors.cmVCU
    elif isinstance(UNITS,(list,np.ndarray)):
        #-- custom units
        dfactor = np.copy(UNITS)
    else:
        raise ValueError('Unknown units {0}'.format(UNITS))

    #-- calculate arrays for clenshaw summations over colatitudes
    s_m_c = np.zeros((npts,LMAX*2+2))
    for m in range(LMAX, -1, -1):
        #-- convolve harmonics with unit factors and smoothing
        s_m_c[:,2*m:2*m+2] = clenshaw_s_m(t, dfactor*wl, m, clm, slm,
            LMAX, SCALE=SCALE)

    #-- calculate cos(phi)
    cos_phi_2 = 2.0*np.cos(phi)
    #-- matrix of cos/sin m*phi summation
    cos_m_phi = np.zeros((npts,LMAX+2),dtype=ASTYPE)
    sin_m_phi = np.zeros((npts,LMAX+2),dtype=ASTYPE)
    #-- initialize matrix with values at lmax+1 and lmax
    cos_m_phi[:,LMAX+1] = np.cos(ASTYPE(LMAX + 1)*phi)
    sin_m_phi[:,LMAX+1] = np.sin(ASTYPE(LMAX + 1)*phi)
    cos_m_phi[:,LMAX] = np.cos(ASTYPE(LMAX)*phi)
    sin_m_phi[:,LMAX] = np.sin(ASTYPE(LMAX)*phi)
    #-- calculate summation for order LMAX
    s_m = s_m_c[:,2*LMAX]*cos_m_phi[:,LMAX] + s_m_c[:,2*LMAX+1]*sin_m_phi[:,LMAX]
    #-- iterate to calculate complete summation
    for m in range(LMAX-1, 0, -1):
        cos_m_phi[:,m] = cos_phi_2*cos_m_phi[:,m+1] - cos_m_phi[:,m+2]
        sin_m_phi[:,m] = cos_phi_2*sin_m_phi[:,m+1] - sin_m_phi[:,m+2]
        #-- calculate summation for order m
        a_m = np.sqrt((2.0*m+3.0)/(2.0*m+2.0))
        s_m = a_m*u*s_m + s_m_c[:,2*m]*cos_m_phi[:,m] + s_m_c[:,2*m+1]*sin_m_phi[:,m]
    #-- calculate spatial field
    spatial = np.sqrt(3.0)*u*s_m + s_m_c[:,0]
    #-- return the calculated spatial field
    return spatial

#-- PURPOSE: compute conditioned arrays for Clenshaw summation from the
#-- fully-normalized associated Legendre's function for an order m
def clenshaw_s_m(t, f, m, clm1, slm1, lmax, ASTYPE=np.longdouble, SCALE=1e-280):
    """
    Compute conditioned arrays for Clenshaw summation from the fully-normalized
    associated Legendre's function for an order m

    Parameters
    ----------
    t: elements ranging from -1 to 1, typically cos(th)
    f: degree dependent factors
    m: spherical harmonic order
    clm1: cosine spherical harmonics
    slm1: sine spherical harmonics
    lmax: maximum spherical harmonic degree

    Returns
    -------
    s_m_c: conditioned array for clenshaw summation
    """
    #-- allocate for output matrix
    N = len(t)
    s_m = np.zeros((N,2),dtype=ASTYPE)
    #-- scaling to prevent overflow
    clm = SCALE*clm1.astype(ASTYPE)
    slm = SCALE*slm1.astype(ASTYPE)
    #-- convert lmax and m to float
    lm = ASTYPE(lmax)
    mm = ASTYPE(m)
    if (m == lmax):
        s_m[:,0] = f[lmax]*clm[lmax,lmax]
        s_m[:,1] = f[lmax]*slm[lmax,lmax]
    elif (m == (lmax-1)):
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t
        s_m[:,0] = a_lm*f[lmax]*clm[lmax,lmax-1] + f[lmax-1]*clm[lmax-1,lmax-1]
        s_m[:,1] = a_lm*f[lmax]*slm[lmax,lmax-1] + f[lmax-1]*slm[lmax-1,lmax-1]
    elif ((m <= (lmax-2)) and (m >= 1)):
        s_mm_c_pre_2 = f[lmax]*clm[lmax,m]
        s_mm_s_pre_2 = f[lmax]*slm[lmax,m]
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t
        s_mm_c_pre_1 = a_lm*s_mm_c_pre_2 + f[lmax-1]*clm[lmax-1,m]
        s_mm_s_pre_1 = a_lm*s_mm_s_pre_2 + f[lmax-1]*slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = ASTYPE(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*t
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + f[l]*clm[l,m]
            s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + f[l]*slm[l,m]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_s_pre_2 = np.copy(s_mm_s_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
            s_mm_s_pre_1 = np.copy(s_mm_s)
        s_m[:,0] = np.copy(s_mm_c)
        s_m[:,1] = np.copy(s_mm_s)
    elif (m == 0):
        s_mm_c_pre_2 = f[lmax]*clm[lmax,0]
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))*t
        s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + f[lmax-1]*clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = ASTYPE(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0)*(ll+1.0)))*t
            b_lm=np.sqrt(((2.0*ll+5.0)*(ll+1.0)*(ll+1.0))/((ll+2.0)*(ll+2.0)*(2.0*ll+1.0)))
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + f[l]*clm[l,0]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
        s_m[:,0] = np.copy(s_mm_c)
    #-- return s_m rescaled with scalef
    return s_m/SCALE
