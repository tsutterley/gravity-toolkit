#!/usr/bin/env python
u"""
gen_stokes.py
Written by Tyler Sutterley (04/2020)

Returns a series of spherical harmonics for an input spatial field

CALLING SEQUENCE:
    Ylms = gen_stokes(data, lon, lat, UNITS=1, LMIN=0, LMAX=60, LOVE=(hl,kl,ll))

INPUTS:
    data: data matrix
    lon: longitude array
    lat: latitude array

OUTPUTS:
    clm: Cosine spherical harmonic coefficients (geodesy normalization)
    slm: Sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

OPTIONS:
    UNITS: input data units
        1: cm of water thickness (default)
        2: Gtons of mass
        3: kg/m^2
    LMIN: Lower bound of Spherical Harmonic Degrees (default = 0)
    LMAX: Upper bound of Spherical Harmonic Degrees (default = 60)
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    PLM: input Legendre polynomials (for improving computational time)
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: computes fully-normalized associated Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units

UPDATE HISTORY:
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
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.units import units

def gen_stokes(data, lon, lat, UNITS=1, LMIN=0, LMAX=60, MMAX=None,
    PLM=None, LOVE=None):

    #-- converting LMIN and LMAX to integer
    LMIN = np.int(LMIN)
    LMAX = np.int(LMAX)
    #-- upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if (MMAX is None) else MMAX

    #-- grid dimensions
    nlon = np.int(len(lon))
    nlat = np.int(len(lat))

    #-- grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    #-- longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    #-- colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    #-- reformatting longitudes to range 0:360 (if previously -180:180)
    lon = np.squeeze(lon.copy())
    if np.any(lon < 0):
        lon_ind, = np.nonzero(lon < 0)
        lon[lon_ind] += 360.0
    #-- Longitude in radians
    phi = lon[np.newaxis,:]*np.pi/180.0
    #-- Colatitude in radians
    th = (90.0 - np.squeeze(lat.copy()))*np.pi/180.0

    #-- reforming data to lonXlat if input latXlon
    sz = np.shape(data)
    data = data.T if (sz[0] == nlat) else np.copy(data)

    #-- SH Degree dependent factors to convert into geodesy normalized SH's
    #-- use splat operator to extract arrays of kl, hl, and ll Love Numbers
    factors = units(lmax=LMAX).spatial(*LOVE)

    #-- extract degree dependent factor for specific units
    #-- calculate integration factors for theta and phi
    #-- Multiplying sin(th) with differentials of theta and phi
    #-- to calculate the integration factor at each latitude
    int_fact = np.zeros((nlat))
    if (UNITS == 1):
        #-- Default Parameter: Input in cm w.e. (g/cm^2)
        dfactor = factors.cmwe
        int_fact[:] = np.sin(th)*dphi*dth
    elif (UNITS == 2):
        #-- Input in gigatonnes (Gt)
        dfactor = factors.cmwe
        #-- rad_e: Average Radius of the Earth [cm]
        int_fact[:] = 1.e15/(factors.rad_e**2)
    elif (UNITS == 3):
        #-- Input in kg/m^2 (mm w.e.)
        dfactor = factors.mmwe
        int_fact[:] = np.sin(th)*dphi*dth
    else:
        #-- default is cm w.e. (g/cm^2)
        dfactor = factors.cmwe
        int_fact[:] = np.sin(th)*dphi*dth

    #-- Calculating cos/sin of phi arrays
    #-- output [m,phi]
    #-- LMAX+1 as there are LMAX+1 elements between 0 and LMAX
    m = np.arange(MMAX+1)
    ccos = np.cos(np.dot(m[:,np.newaxis],phi))
    ssin = np.sin(np.dot(m[:,np.newaxis],phi))

    #-- Calculates fully-normalized Legendre Polynomials with plm_holmes.py
    #-- Output is plm[l,m,th]
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    #-- added option to precompute plms to improve computational speed
    if PLM is None:
        #-- if plms are not pre-computed: calculate Legendre polynomials
        PLM,dPLM = plm_holmes(LMAX,np.cos(th))

    #-- Multiplying by integration factors [sin(theta)*dtheta*dphi]
    #-- truncate legendre polynomials to spherical harmonic order MMAX
    for j in range(0,nlat):
        plm[:,m,j] = PLM[:,m,j]*int_fact[j]

    #-- Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing output spherical harmonic matrices
    clm = np.zeros((LMAX+1,MMAX+1))
    slm = np.zeros((LMAX+1,MMAX+1))
    #-- Multiplying gridded data with sin/cos of m#phis
    #-- This will sum through all phis in the dot product
    #-- output [m,theta]
    dcos = np.dot(ccos,data)
    dsin = np.dot(ssin,data)
    for l in range(LMIN,LMAX+1):#-- equivalent to LMIN:LMAX
        mm = np.min([MMAX,l])#-- truncate to MMAX if specified (if l > MMAX)
        m = np.arange(0,mm+1)#-- mm+1 elements between 0 and mm
        #-- Summing product of plms and data over all latitudes
        #-- axis=1 signifies the direction of the summation
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        #-- Multiplying by factors to convert to geodesy normalized coefficients
        clm[l,m] = dfactor[l]*yclm[l,m]
        slm[l,m] = dfactor[l]*yslm[l,m]

    return {'clm':clm, 'slm':slm, 'l':np.arange(LMAX+1), 'm':np.arange(MMAX+1)}
