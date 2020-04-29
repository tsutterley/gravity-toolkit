#!/usr/bin/env python
u"""
plm_mohlenkamp.py
Written by Tyler Sutterley (05/2015)

Computes fully-normalized associated Legendre Polynomials
    for an array of x values
Uses Martin Mohlenkamp's recursion relation derived from the
    Szego (1939) Recurrence formula for Jacobi Polynomials (Pg 71)

    A User's Guide to Spherical Harmonics
    http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf

With this algorithm, the associated Legendre Functions are
    constructed as an amplitude times a Jacobi Polynomial
    P[l,m](cos(theta)) = (sin(theta)^2)*J[l-m,m,m](cos(theta))

CALLING SEQUENCE:
    plm = plm_mohlenkamp(LMAX, np.cos(theta))

INPUTS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    x: typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    plm: Legendre polynomials (geodesy normalization)

OPTIONS:
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

NOTES:
    Modified and updated from IDL plm_x.pro coded by Sean Swenson
    Difference from martin function in geoid_mk.mac.f:
        plms from plm_mohlenkamp are normalized inside the function
        plms from martin are normalized outside the function
    For large spherical harmonic degrees this recurrence relation
        is poorly conditioned
    For spherical harmonic orders above ~1000 can cause overflows

UPDATE HISTORY:
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 09/2013
"""
import numpy as np

def plm_mohlenkamp(LMAX, x, MMAX=None):

    #-- Verify LMAX as integer
    LMAX = np.int(LMAX)
    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- size of the x array
    #-- x is an array
    if (np.ndim(x) > 0):
        sx = np.shape(x)[0]
    else:#-- x is a single value
        sx =1

    #-- Initialize the output Legendre polynomials
    plm=np.zeros((LMAX+1,MMAX+1,sx))
    #-- Jacobi polynomial for the recurrence relation
    jlmm=np.zeros((LMAX+1,MMAX+1,sx))
    #-- for x=cos(th): rsin= sin(th)
    rsin=np.sqrt(1.0 - x**2)

    #-- for all spherical harmonic orders of interest
    for mm in range(0,MMAX+1):#-- equivalent to 0:MMAX
        #-- Initialize the recurrence relation
        #-- J-1,m,m Term == 0
        #-- J0,m,m Term
        if (mm > 0):
            #-- j ranges from 1 to mm for the product
            j = np.arange(0,mm)+1.0
            jlmm[0,mm,:] = np.prod(np.sqrt(1.0 + 1.0/(2.0*j)))/np.sqrt(2.0)
        else: #-- if mm == 0: jlmm = 1/sqrt(2)
            jlmm[0,mm,:] = 1.0/np.sqrt(2.0)
        #-- Jk,m,m Terms
        for k in range (1,(LMAX+1)):#-- computation for SH degrees
            #-- Initialization begins at -1
            #-- this is to make the formula parallel the function written in
            #-- Martin Mohlenkamp's Guide to Spherical Harmonics
            #-- Jacobi General Terms
            if (k == 1):#-- for degree 1 terms
                jlmm[k,mm,:] = 2.0*x * jlmm[k-1,mm,:] * \
                    np.sqrt(1.0 + (mm - 0.5)/k) * \
                    np.sqrt(1.0 - (mm - 0.5)/(k + 2.0*mm))
            else:#-- for all other spherical harmonic degrees
                jlmm[k,mm,:] = 2.0*x * jlmm[k-1,mm,:] * \
                    np.sqrt(1.0 + (mm - 0.5)/k) * \
                    np.sqrt(1.0 - (mm - 0.5)/(k + 2.0*mm)) - \
                    jlmm[k-2,mm,:] * np.sqrt(1.0 + 4.0/(2.0*k + 2.0*mm - 3.0)) * \
                    np.sqrt(1.0 - (1.0/k)) * np.sqrt(1.0 - 1.0/(k + 2.0*mm))
        #-- Normalization is geodesy convention
        for l in range(mm,LMAX+1): #-- equivalent to mm:LMAX
            if (mm == 0):#-- Geodesy normalization (m=0) == sqrt(2)*sin(th)^0
                #-- rsin^mm term is dropped as rsin^0 = 1
                plm[l,mm,:] = np.sqrt(2.0)*jlmm[l-mm,mm,:]
            else:#-- Geodesy normalization all others == 2*sin(th)^mm
                plm[l,mm,:] = 2.0*(rsin**mm)*jlmm[l-mm,mm,:]
    return plm
