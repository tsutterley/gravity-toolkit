#!/usr/bin/env python
u"""
gen_disc_load.py (04/2020)
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
    PLM: input Legendre polynomials (for improving computational time)
        input PLMs are the polynomials for cos(theta) (disc center)
        not for the spherical cap alpha plms
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    legendre_polynomials.py: Computes fully normalized Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units

REFERENCE:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)

UPDATE HISTORY:
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
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.legendre_polynomials import legendre_polynomials
from gravity_toolkit.units import units

def gen_disc_load(data,lon,lat,area,LMAX=60,MMAX=None,PLM=None,LOVE=None):

    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- Earth Parameters
    factors = units(lmax=LMAX)
    rho_e = factors.rho_e#-- Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e#-- Average Radius of the Earth [cm]

    #-- convert lon and lat to radians
    phi = lon*np.pi/180.0#-- Longitude in radians
    th = (90.0 - lat)*np.pi/180.0#-- Colatitude in radians

    #-- convert input area into cm^2 and then divide by area of a half sphere
    #-- alpha will be 1 - the ratio of the input area with the half sphere
    alpha = (1.0 - 1e10*area/(2.0*np.pi*rad_e**2))

    #-- Input data is in gigatonnes (Gt)
    #-- 1e15 converts from Gt to grams, 1e10 converts from km^2 to cm^2
    unit_conv = 1e15/(1e10*area)

    #-- Coefficient for calculating Stokes coefficients for a disc load
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

    #-- Calculating plms of the disc
    #-- allocating for constructed array
    pl_alpha = np.zeros((LMAX+1))
    #-- l=0 is a special case (P(-1) = 1, P(1) = cos(alpha))
    pl_alpha[0] = (1.0 - alpha)/2.0
    #-- for all other degrees: calculate the legendre polynomials up to LMAX+1
    pl_matrix, dpl_matrix = legendre_polynomials(LMAX+1,alpha)
    for l in range(1, LMAX+1):#-- LMAX+1 to include LMAX
        #-- from Longman (1962) and Jacob et al (2012)
        #-- unnormalizing Legendre polynomials
        #-- sqrt(2*l - 1) == sqrt(2*(l-1) + 1)
        #-- sqrt(2*l + 3) == sqrt(2*(l+1) + 1)
        pl_lower = pl_matrix[l-1]/np.sqrt(2.0*l-1.0)
        pl_upper = pl_matrix[l+1]/np.sqrt(2.0*l+3.0)
        pl_alpha[l] = (pl_lower - pl_upper)/2.0

    #-- Calculate Legendre Polynomials using Holmes and Featherstone relation
    #-- this would be the plm for the center of the disc load
    #-- used to rotate the disc load to point lat/lon
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
    Ylms = {}
    Ylms['l'] = np.arange(LMAX+1)
    Ylms['m'] = np.arange(MMAX+1)
    Ylms['clm'] = np.zeros((LMAX+1,MMAX+1))
    Ylms['slm'] = np.zeros((LMAX+1,MMAX+1))
    for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
        l = np.arange(m,LMAX+1)#-- LMAX+1 to include LMAX
        #-- rotate disc load to be centered at lat/lon
        plm[l,m] = plmout[l,m]*pl_alpha[l]
        #-- multiplying clm by cos(m*phi) and slm by sin(m*phi)
        #-- to get a field of spherical harmonics
        yclm[l,m] = plm[l,m]*dcos[m]
        yslm[l,m] = plm[l,m]*dsin[m]
        #-- multiplying by coefficients to convert to geoid coefficients
        Ylms['clm'][l,m] = coeff*dfactor[l]*yclm[l,m]
        Ylms['slm'][l,m] = coeff*dfactor[l]*yslm[l,m]

    #-- return the output spherical harmonics
    return Ylms
