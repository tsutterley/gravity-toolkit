#!/usr/bin/env python
u"""
gen_harmonics.py
Written by Tyler Sutterley (04/2022)
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
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
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
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.fourier_legendre import fourier_legendre

def gen_harmonics(data, lon, lat, **kwargs):
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
    PLM: float, default 0
        Fully normalized associated Legendre polynomials or
        Fourier coefficients of Legendre polynomials
    METHOD: str, default 'integration'
        Conversion method for calculating harmonics

            - ``'integration'``: for global grids
            - ``'fourier'``: for regional or global grids

    Returns
    -------
    clm: float
        cosine spherical harmonic coefficients (4-pi normalized)
    slm: float
        sine spherical harmonic coefficients (4-pi normalized)
    l: int
        spherical harmonic degree to LMAX
    m: int
        spherical harmonic order to MMAX
    """
    #-- set default keyword arguments
    kwargs.setdefault('LMAX',60)
    kwargs.setdefault('MMAX',None)
    kwargs.setdefault('PLM',0)
    kwargs.setdefault('METHOD','integration')
    #-- upper bound of spherical harmonic orders (default = LMAX)
    if kwargs['MMAX'] is None:
        kwargs['MMAX'] = np.copy(kwargs['LMAX'])
    #-- convert latitude and longitude to float if integers
    lon = lon.astype(np.float64)
    lat = lat.astype(np.float64)
    #-- reforming data to lonXlat if input latXlon
    sz = np.shape(data)
    dinput = np.transpose(data) if (sz[0] == len(lat)) else np.copy(data)
    #-- convert spatial field into spherical harmonics
    if (kwargs['METHOD'].lower() == 'integration'):
        Ylms = integration(dinput, lon, lat, **kwargs)
    elif (kwargs['METHOD'].lower() == 'fourier'):
        Ylms = fourier(dinput, lon, lat, **kwargs)
    #-- return the output spherical harmonics object
    return Ylms

def integration(data, lon, lat, LMAX=60, MMAX=None, PLM=0, **kwargs):
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

    #-- dimensions of the longitude and latitude arrays
    nlon = np.int64(len(lon))
    nlat = np.int64(len(lat))
    #-- grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    #-- longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    #-- colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    #-- reformatting longitudes to range 0:360 (if previously -180:180)
    if np.count_nonzero(lon < 0):
        lon[lon < 0] += 360.0
    #-- calculate longitude and colatitude arrays in radians
    phi = np.reshape(lon,(1,nlon))*np.pi/180.0#-- reshape to 1xnlon
    th = (90.0 - np.squeeze(lat))*np.pi/180.0#-- remove singleton dimensions

    #-- Calculating cos/sin of phi arrays (output [m,phi])
    #-- LMAX+1 as there are LMAX+1 elements between 0 and LMAX
    m = np.arange(MMAX+1)[:, np.newaxis]
    ccos = np.cos(np.dot(m,phi))
    ssin = np.sin(np.dot(m,phi))

    #-- Multiplying sin(th) with differentials of theta and phi
    #-- to calculate the integration factor at each latitude
    int_fact = np.sin(th)*dphi*dth
    coeff = 1.0/(4.0*np.pi)

    #-- Calculate polynomials using Holmes and Featherstone (2002) relation
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    if (np.ndim(PLM) == 0):
        plmout,dplm = plm_holmes(LMAX,np.cos(th))
    else:
        #-- use precomputed plms to improve computational speed
        #-- or to use a different recursion relation for polynomials
        plmout = PLM

    #-- Multiply plms by integration factors [sin(theta)*dtheta*dphi]
    #-- truncate plms to maximum spherical harmonic order if MMAX < LMAX
    m = np.arange(MMAX+1)
    for j in range(0,nlat):
        plm[:,m,j] = plmout[:,m,j]*int_fact[j]

    #-- Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    #-- Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    #-- Multiplying gridded data with sin/cos of m#phis (output [m,theta])
    #-- This will sum through all phis in the dot product
    dcos = np.dot(ccos,data)
    dsin = np.dot(ssin,data)
    for l in range(0,LMAX+1):
        mm = np.min([MMAX,l])#-- truncate to MMAX if specified (if l > MMAX)
        m = np.arange(0,mm+1)#-- mm+1 elements between 0 and mm
        #-- Summing product of plms and data over all latitudes
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        #-- convert to output normalization (4-pi normalized harmonics)
        Ylms.clm[l,m] = coeff*yclm[l,m]
        Ylms.slm[l,m] = coeff*yslm[l,m]

    #-- return the output spherical harmonics object
    return Ylms

def fourier(data, lon, lat, LMAX=60, MMAX=None, PLM=0, **kwargs):
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

    #-- dimensions of the longitude and latitude arrays
    nlon = np.int64(len(lon))
    nlat = np.int64(len(lat))
    #-- remove singleton dimensions and convert to radians
    phi = (np.squeeze(lon)*np.pi/180.0)
    #-- Colatitude in radians
    theta = ((90.0 - np.squeeze(lat))*np.pi/180.0)

    #-- MMAX+1 to include MMAX
    mm = np.arange(MMAX+1)[:, np.newaxis]
    #-- Calculate cos and sin coefficients of signal
    ccos = np.cos(np.dot(mm,phi[np.newaxis,:]))
    ssin = np.sin(np.dot(mm,phi[np.newaxis,:]))
    dcos = np.dot(ccos,data)
    dsin = np.dot(ssin,data)

    #-- Normalize fourier coefficients
    dcos[0,:] = dcos[0,:]/nlon
    dcos[1:MMAX+1,:] = 2.0*dcos[1:MMAX+1,:]/nlon
    dsin[0,:] = dsin[0,:]/nlon
    dsin[1:MMAX+1,:] = 2.0*dsin[1:MMAX+1,:]/nlon

    #-- Calculate cos and sin coefficients of theta component
    #-- Because the function is defined on (0,pi)
    #-- it can be expanded in just cosine terms.
    #-- this routine assumes that 0 and pi are not included
    theta_cc = np.zeros((MMAX+1,MMAX+1))
    theta_sc = np.zeros((MMAX+1,MMAX+1))
    m_even = np.arange(0,MMAX+1,2)
    m_odd = np.arange(1,MMAX,2)
    n_even = len(m_even)
    n_odd = len(m_odd)

    if np.isclose([theta[0],theta[nlat-1]],[0.0,np.pi]).all():
        #-- non-endpoints
        nt = np.dot(mm,theta[1:nlat-1][np.newaxis,:])
        theta_cc[m_even,:] = 2.0*np.dot(dcos[m_even,1:nlat-1],np.cos(nt).T)
        theta_sc[m_even,:] = 2.0*np.dot(dsin[m_even,1:nlat-1],np.cos(nt).T)
        theta_cc[m_odd,:] = 2.0*np.dot(dcos[m_odd,1:nlat-1],np.sin(nt).T)
        theta_sc[m_odd,:] = 2.0*np.dot(dsin[m_odd,1:nlat-1],np.sin(nt).T)

        #-- endpoints
        theta_cc[m_even,:] += np.dot((dcos[m_even,0]*np.cos(theta[0]) +
            dcos[m_even,nlat-1]*np.cos(theta[nlat-1]))[:,np.newaxis], mm.T)
        theta_sc[m_even,:] += np.dot((dsin[m_even,0]*np.cos(theta[0]) +
            dsin[m_even,nlat-1]*np.cos(theta[nlat-1]))[:,np.newaxis], mm.T)
        theta_cc[m_odd,:] += np.dot((dcos[m_odd,0]*np.sin(theta[0]) +
            dcos[m_odd,nlat-1]*np.sin(theta[nlat-1]))[:,np.newaxis], mm.T)
        theta_sc[m_odd,:] += np.dot((dsin[m_odd,0]*np.sin(theta[0]) +
            dsin[m_odd,nlat-1]*np.sin(theta[nlat-1]))[:,np.newaxis], mm.T)

    elif not np.isclose([theta[0],theta[nlat-1]],[0.0,np.pi]).any():
        nt = np.dot(mm,theta[np.newaxis,:])
        theta_cc[m_even,:] = 2.0*np.dot(dcos[m_even,:],np.cos(nt).T)
        theta_sc[m_even,:] = 2.0*np.dot(dsin[m_even,:],np.cos(nt).T)
        theta_cc[m_odd,:] = 2.0*np.dot(dcos[m_odd,:],np.sin(nt).T)
        theta_sc[m_odd,:] = 2.0*np.dot(dsin[m_odd,:],np.sin(nt).T)
    else:
        raise ValueError('Latitude coordinates incompatible')

    #-- Normalize theta fourier coefficients
    theta_cc[:,0] = theta_cc[:,0]/(2.0*nlat)
    theta_cc[:,1:MMAX+1] = theta_cc[:,1:MMAX+1]/nlat
    theta_sc[:,0] = theta_sc[:,0]/(2.0*nlat)
    theta_sc[:,1:MMAX+1] = theta_sc[:,1:MMAX+1]/nlat

    #-- Correct normalization for the incomplete coverage of the sphere
    delphi = np.abs(phi[1]-phi[0])
    deltheta = np.abs(theta[1]-theta[0])
    norm = nlon*delphi/(2.0*np.pi)*nlat*deltheta/np.pi
    theta_cc = theta_cc*norm
    theta_sc = theta_sc*norm

    #-- Calculate cos and sin coefficients of Legendre functions
    #-- Expand m = even terms in a cosine series
    #-- Expand m = odd terms in a sine series
    #-- Both are stride 2
    if (np.ndim(PLM) == 0):
        plm = fourier_legendre(LMAX,MMAX)
    else:
        #-- use precomputed plms to improve computational speed
        plm = PLM

    #-- Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))

    #-- Sum theta fourier coefficients
    #-- temp is the integral of cos(n theta) cos(k theta) dcos(theta)
    #-- over the interval 0 to pi
    #-- n and k must have like parities

    #-- m = even terms
    k_even = np.zeros((n_even,n_even))
    for n in range(0,MMAX+2,2):
        k_even[:,n//2] = 0.5*(1.0/(1.0-m_even-n) + 1.0/(1.0+m_even-n) +
            1.0/(1.0-m_even+n) + 1.0/(1.0+m_even+n))

    k_odd = np.zeros((n_odd,n_odd))
    for n in range(1,MMAX+1,2):
        k_odd[:,(n-1)//2] = 0.5*(1.0/(1-m_odd-n) + 1.0/(1+m_odd-n) +
            1.0/(1-m_odd+n) + 1.0/(1+m_odd+n))

    #-- calculate spherical harmonics for m == even terms
    l_even = np.arange(0,LMAX+1,2)
    l_odd = np.arange(1,LMAX,2)
    for m in range(0,MMAX+2,2):
        temp = np.dot(plm[l_even,m,m_even[:,np.newaxis]].T,k_even)
        Ylms.clm[l_even,m] = np.dot(theta_cc[m,m_even[:,np.newaxis]].T,temp.T)
        Ylms.slm[l_even,m] = np.dot(theta_sc[m,m_even[:,np.newaxis]].T,temp.T)
        temp = np.dot(plm[l_odd,m,m_odd[:,np.newaxis]].T,k_odd)
        Ylms.clm[l_odd,m] = np.dot(theta_cc[m,m_odd[:,np.newaxis]].T,temp.T)
        Ylms.slm[l_odd,m] = np.dot(theta_sc[m,m_odd[:,np.newaxis]].T,temp.T)

    #-- m = odd terms
    k_even = np.zeros((n_even,n_even))
    for n in range(0,MMAX+2,2):
        k_even[:,n//2] = 0.5*(-1.0/(1-m_even-n) + 1.0/(1.0+m_even-n) +
            1.0/(1.0-m_even+n) - 1.0/(1.0+m_even+n))

    k_odd = np.zeros((n_odd,n_odd))
    for n in range(1,MMAX+1,2):
        k_odd[:,(n-1)//2] = 0.5*(-1.0/(1-m_odd-n) + 1.0/(1.0+m_odd-n) +
            1.0/(1.0-m_odd+n) - 1.0/(1.0+m_odd+n))

    #-- calculate spherical harmonics for m == odd terms
    l_even = np.arange(2,LMAX+1,2)#-- do not in include l=0
    l_odd = np.arange(1,LMAX,2)
    for m in range(1,MMAX+1,2):
        temp = np.dot(plm[l_even,m,m_even[:,np.newaxis]].T,k_even)
        Ylms.clm[l_even,m] = np.dot(theta_cc[m,m_even[:,np.newaxis]].T,temp.T)
        Ylms.slm[l_even,m] = np.dot(theta_sc[m,m_even[:,np.newaxis]].T,temp.T)
        temp = np.dot(plm[l_odd,m,m_odd[:,np.newaxis]].T,k_odd)
        Ylms.clm[l_odd,m] = np.dot(theta_cc[m,m_odd[:,np.newaxis]].T,temp.T)
        Ylms.slm[l_odd,m] = np.dot(theta_sc[m,m_odd[:,np.newaxis]].T,temp.T)

    #-- Divide by Plm normalization
    Ylms.clm[:,0] /= 2.0
    Ylms.slm[:,0] /= 2.0
    Ylms.clm[:,1:MMAX+1] /= 4.0
    Ylms.slm[:,1:MMAX+1] /= 4.0

    #-- return the output spherical harmonics object
    return Ylms
