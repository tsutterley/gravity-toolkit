#!/usr/bin/env python
u"""
gen_spherical_cap.py
Written by Tyler Sutterley (06/2023)
Calculates gravitational spherical harmonic coefficients for a spherical cap

Creating a spherical cap with generating angle alpha is a 2 step process:
    1) obtain harmonics when cap is located at the north pole
    2) rotate the cap to an arbitrary latitude and longitude

CALLING SEQUENCE:
    Ylms = gen_spherical_cap(data, lon, lat, LMAX=LMAX, RAD_CAP=RAD_CAP)

INPUTS:
    data: data magnitude
    lon: longitude of spherical cap center
    lat: latitude of spherical cap center

OUTPUTS:
    clm: cosine spherical harmonic coefficients (geodesy normalization)
    slm: sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    RAD_CAP: spherical cap radius in degrees
    RAD_KM: spherical cap radius in kilometers
    AREA: spherical cap area in cm^2
    UNITS: input data units
        1: cm of water thickness (default)
        2: gigatonnes of mass
        3: kg/m^2
        list: custom unit conversion factor
    PLM: input Legendre polynomials
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully-normalized associated
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
    Updated 11/2021: added UNITS list option for converting from custom units
    Updated 07/2020: added function docstrings
    Updated 05/2020: vectorize calculation over degrees to improve compute time
    Updated 04/2020: reading load love numbers outside of this function
        include degrees and orders in output dictionary for harmonics class
        use units class for Earth parameters
    Updated 04/2018: Using Holmes and Featherstone relation for Plms
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Updated 07/2017: outputs of legendre_polynomials.py include derivatives now
    Updated 08/2015: added error handling for common exceptions
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
        added parameter RAD_KM for the cap radius in kilometers (alternative)
        added parameter UNITS to be similar to gen_stokes program and
        for reconstructing the spatial fields of recovered mascon time series
    Updated 03/2014: using Legendre polynomials code to generate the zonal harmonics
        rather than using the first array of plm_mohlenkamp.
        added some updates to comments
    Updated 10/2013: major reorganization for computational efficiency
    Updated 10/2013: added option for entering cap radius instead of AREA
    Updated 05/2013: added option to precompute plms
    Updated 05/2013: adapted for python (changed UNGRID and PTMASS options to Y/N)
    Updated 06/2012: major revision to code organizzation
    Written 04/2012
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.legendre_polynomials import legendre_polynomials

def gen_spherical_cap(data, lon, lat, LMAX=60, MMAX=None,
    AREA=0, RAD_CAP=0, RAD_KM=0, UNITS=1, PLM=None, LOVE=None):
    r"""
    Calculates spherical harmonic coefficients for a spherical cap
    :cite:p:`Holmes:2002ff,Longman:1962ev,Farrell:1972cm,Pollack:1973gi,Jacob:2012eo`

    Parameters
    ----------
    data: np.ndarray
        data magnitude
    lon: np.ndarray
        longitude of spherical cap center
    lat: np.ndarray
        latitude of spherical cap center
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    AREA: Float, default 0
        Area of spherical cap (cm\ :sup:`2`)
    UNITS: int, default 1
        Input data units

            - ``1``: cm water equivalent thickness (cm w.e., g/cm\ :sup:`2`)
            - ``2``: gigatonnes of mass (Gt)
            - ``3``: mm water equivalent thickness (mm w.e., kg/m\ :sup:`2`)
            - list: custom unit conversion factor
    PLM: np.ndarray, default 0
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

    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # convert lon and lat to radians
    phi = lon*np.pi/180.0# Longitude in radians
    th = (90.0 - lat)*np.pi/180.0# Colatitude in radians

    # Earth Parameters
    factors = gravity_toolkit.units(lmax=LMAX)

    # Converting input area into an equivalent spherical cap radius
    # Following Jacob et al. (2012) Equation 4 and 5
    # alpha is the vertical semi-angle subtending a cone at the
    # center of the earth
    if (RAD_CAP != 0):
        # if given spherical cap radius in degrees
        # converting to radians
        alpha = RAD_CAP*np.pi/180.0
    elif (AREA != 0):
        # if given spherical cap area in cm^2
        # radius in centimeters
        radius_cm = np.sqrt(AREA/np.pi)
        # Calculating angular radius of spherical cap
        alpha = (radius_cm/factors.rad_e)
    elif (RAD_KM != 0):
        # if given spherical cap radius in kilometers
        # Calculating angular radius of spherical cap
        alpha = (1e5*RAD_KM)/factors.rad_e
    else:
        raise ValueError('Input RAD_CAP, AREA or RAD_KM of spherical cap')

    # Calculate factor to convert from input units
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
        # calculate spherical cap area from angular radius
        area = np.pi*(alpha*factors.rad_e)**2
        # the 1.e15 converts from gigatons/cm^2 to cm of water
        # 1 g/cm^3 = 1000 kg/m^3 = density water
        # 1 Gt = 1 Pg = 1.e15 g
        unit_conv = 1.e15/area
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

    # Calculating plms of the spherical caps
    # From Longman et al. (1962)
    # pl_alpha = F(alpha) from Jacob 2011
    # pl_alpha is purely zonal and depends only on the size of the cap
    # allocating for constructed array
    pl_alpha = np.zeros((LMAX+1))
    # l=0 is a special case (P(-1) = 1, P(1) = cos(alpha))
    pl_alpha[0] = (1.0 - np.cos(alpha))/2.0
    # for all other degrees: calculate the legendre polynomials up to LMAX+1
    pl_matrix,_ = legendre_polynomials(LMAX+1,np.cos(alpha))
    for l in range(1, LMAX+1):# LMAX+1 to include LMAX
        # from Longman (1962) and Jacob et al (2012)
        # unnormalizing Legendre polynomials
        # sqrt(2*l - 1) == sqrt(2*(l-1) + 1)
        # sqrt(2*l + 3) == sqrt(2*(l+1) + 1)
        pl_lower = pl_matrix[l-1]/np.sqrt(2.0*l-1.0)
        pl_upper = pl_matrix[l+1]/np.sqrt(2.0*l+3.0)
        pl_alpha[l] = (pl_lower - pl_upper)/2.0

    # Calculating Legendre Polynomials
    # added option to precompute plms to improve computational speed
    # this would be the plm for the center of the spherical cap
    # used to rotate the spherical cap to point lat/lon
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
        # rotate spherical cap to be centered at lat/lon
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
