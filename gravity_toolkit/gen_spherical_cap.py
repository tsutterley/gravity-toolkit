#!/usr/bin/env python
u"""
gen_spherical_cap.py
Written by Tyler Sutterley (04/2022)
Calculates gravitational spherical harmonic coefficients for a spherical cap

Spherical cap derivation from Longman (1962), Farrell (1972), Pollack (1973)
    and Jacob (2012)

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
    plm_holmes.py: Computes fully-normalized associated Legendre polynomials
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
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.legendre_polynomials import legendre_polynomials

def gen_spherical_cap(data, lon, lat, LMAX=60, MMAX=None,
    AREA=0, RAD_CAP=0, RAD_KM=0, UNITS=1, PLM=None, LOVE=None):
    """
    Calculates spherical harmonic coefficients for a spherical cap

    Parameters
    ----------
    data: float
        data magnitude
    lon: float
        longitude of spherical cap center
    lat: float
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
    PLM: float, default 0
        Input Legendre polynomials
    LOVE: tuple or NoneType, default None
        Input load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

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

    References
    ----------
    .. [Holmes2002] S. A. Holmes and W. E. Featherstone,
        "A unified approach to the Clenshaw summation and the recursive
        computation of very high degree and order normalised associated
        Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002).
        `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_
    .. [Longman1962] I. M. Longman, "A Green's function for determining
        the deformation of the Earth under surface mass loads: 1. Theory",
        *Journal of Geophysical Research*, 67(2), (1962).
        `doi: 10.1029/JZ067i002p00845 <https://doi.org/10.1029/JZ067i002p00845>`_
    .. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads",
        *Reviews of Geophysics and Space Physics*, 10(3), (1972).
        `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_
    .. [Pollack1973] H. N. Pollack, "Spherical harmonic representation of the
        gravitational potential of a point mass, a spherical cap, and a
        spherical rectangle", *Journal of Geophysical Research*, 78(11), (1973).
        `doi: 10.1029/JB078i011p01760 <https://doi.org/10.1029/JB078i011p01760>`_
    .. [Jacob2012] T. Jacob et al., "Estimating geoid height change in North America:
        past, present and future", *Journal of Geodesy*, 86, 337-358, (2012).
        `doi: 10.1007/s00190-011-0522-7 <https://doi.org/10.1007/s00190-011-0522-7>`_
    """

    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- Earth Parameters
    factors = gravity_toolkit.units(lmax=LMAX)
    rho_e = factors.rho_e#-- Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e#-- Average Radius of the Earth [cm]

    #-- convert lon and lat to radians
    phi = lon*np.pi/180.0#-- Longitude in radians
    th = (90.0 - lat)*np.pi/180.0#-- Colatitude in radians

    #-- Converting input area into an equivalent spherical cap radius
    #-- Following Jacob et al. (2012) Equation 4 and 5
    #-- alpha is the vertical semi-angle subtending a cone at the
    #-- center of the earth
    if (RAD_CAP != 0):
        #-- if given spherical cap radius in degrees
        #-- converting to radians
        alpha = RAD_CAP*np.pi/180.0
    elif (AREA != 0):
        #-- if given spherical cap area in cm^2
        #-- radius in centimeters
        radius_cm = np.sqrt(AREA/np.pi)
        #-- Calculating angular radius of spherical cap
        alpha = (radius_cm/rad_e)
    elif (RAD_KM != 0):
        #-- if given spherical cap radius in kilometers
        #-- Calculating angular radius of spherical cap
        alpha = (1e5*RAD_KM)/rad_e
    else:
        raise ValueError('Input RAD_CAP, AREA or RAD_KM of spherical cap')

    #-- Calculate factor to convert from input units into cmH2O equivalent
    #-- Default input is for inputs already in cmH2O (unit_conv = 1)
    if (UNITS == 1):
        #-- Input data is in cm water equivalent (cmH2O)
        unit_conv = 1.0
    elif (UNITS == 2):
        #-- Input data is in gigatonnes (Gt)
        #-- calculate spherical cap area from angular radius
        area = np.pi*(alpha*rad_e)**2
        #-- the 1.e15 converts from gigatons/cm^2 to cm of water
        #-- 1 g/cm^3 = 1000 kg/m^3 = density water
        #-- 1 Gt = 1 Pg = 1.e15 g
        unit_conv = 1.e15/area
    elif (UNITS == 3):
        #-- Input data is in kg/m^2
        #-- 1 kg = 1000 g
        #-- 1 m^2 = 100*100 cm^2 = 1e4 cm^2
        unit_conv = 0.1
    elif isinstance(UNITS,(list,np.ndarray)):
        #-- custom units
        unit_conv = np.copy(UNITS)
    else:
        raise ValueError('Unknown units {0}'.format(UNITS))

    #-- Coefficient for calculating Stokes coefficients for a spherical cap
    #-- From Jacob et al (2012), Farrell (1972) and Longman (1962)
    coeff = 3.0/(rad_e*rho_e)

    #-- extract arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = LOVE

    #-- calculate array of l values ranging from 0 to LMAX (harmonic degrees)
    #-- LMAX+1 as there are LMAX+1 elements between 0 and LMAX
    l = np.arange(LMAX+1)
    #-- calculate SH degree dependent factors to convert from coefficients
    #-- of mass into normalized geoid coefficients
    #-- NOTE: these are not the normal factors for converting to geoid due
    #-- to the square of the denominator
    #-- kl[l] is the Load Love Number of degree l
    dfactor = (1.0 + kl[l])/((1.0 + 2.0*l)**2)

    #-- Calculating plms of the spherical caps
    #-- From Longman et al. (1962)
    #-- pl_alpha = F(alpha) from Jacob 2011
    #-- pl_alpha is purely zonal and depends only on the size of the cap
    #-- allocating for constructed array
    pl_alpha = np.zeros((LMAX+1))
    #-- l=0 is a special case (P(-1) = 1, P(1) = cos(alpha))
    pl_alpha[0] = (1.0 - np.cos(alpha))/2.0
    #-- for all other degrees: calculate the legendre polynomials up to LMAX+1
    pl_matrix,_ = legendre_polynomials(LMAX+1,np.cos(alpha))
    for l in range(1, LMAX+1):#-- LMAX+1 to include LMAX
        #-- from Longman (1962) and Jacob et al (2012)
        #-- unnormalizing Legendre polynomials
        #-- sqrt(2*l - 1) == sqrt(2*(l-1) + 1)
        #-- sqrt(2*l + 3) == sqrt(2*(l+1) + 1)
        pl_lower = pl_matrix[l-1]/np.sqrt(2.0*l-1.0)
        pl_upper = pl_matrix[l+1]/np.sqrt(2.0*l+3.0)
        pl_alpha[l] = (pl_lower - pl_upper)/2.0

    #-- Calculating Legendre Polynomials
    #-- added option to precompute plms to improve computational speed
    #-- this would be the plm for the center of the spherical cap
    #-- used to rotate the spherical cap to point lat/lon
    if PLM is None:
        plmout,dplm = plm_holmes(LMAX,np.cos(th))
        #-- truncate precomputed plms to order
        plmout = np.squeeze(plmout[:,:MMAX+1,:])
    else:
        #-- truncate precomputed plms to degree and order
        plmout = PLM[:LMAX+1,:MMAX+1]

    #-- calculate array of m values ranging from 0 to MMAX (harmonic orders)
    #-- MMAX+1 as there are MMAX+1 elements between 0 and MMAX
    m = np.arange(MMAX+1)
    #-- Multiplying by the units conversion factor (unit_conv) to
    #-- convert from the input units into cmH2O equivalent
    #-- Multiplying point mass data (converted to cmH2O) with sin/cos of m*phis
    #-- data normally is 1 for a uniform 1cm water equivalent layer
    #-- but can be a mass point if reconstructing a spherical harmonic field
    #-- NOTE: NOT a matrix multiplication as data (and phi) is a single point
    dcos = unit_conv*data*np.cos(m*phi)
    dsin = unit_conv*data*np.sin(m*phi)

    #-- Multiplying by plm_alpha (F_l from Jacob 2012)
    plm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
        l = np.arange(m,LMAX+1)#-- LMAX+1 to include LMAX
        #-- rotate spherical cap to be centered at lat/lon
        plm[l,m] = plmout[l,m]*pl_alpha[l]
        #-- multiplying clm by cos(m*phi) and slm by sin(m*phi)
        #-- to get a field of spherical harmonics
        yclm[l,m] = plm[l,m]*dcos[m]
        yslm[l,m] = plm[l,m]*dsin[m]
        #-- multiplying by coefficients to convert to geoid coefficients
        Ylms.clm[l,m] = coeff*dfactor[l]*yclm[l,m]
        Ylms.slm[l,m] = coeff*dfactor[l]*yslm[l,m]

    #-- return the output spherical harmonics object
    return Ylms
