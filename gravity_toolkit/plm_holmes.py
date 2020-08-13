#!/usr/bin/env python
u"""
plm_holmes.py
Written by Tyler Sutterley (08/2020)

Computes fully-normalized associated Legendre Polynomials
    for a vector of x values (can also be singular)

Uses Holmes and Featherstone (2002) recursion relation

This recursion relation is better conditioned for high
    degree and order than the Martin Mohlenkamp relation
It is stable up to very high degree and order (at least 3000).

This is achieved by individually computing the sectorials to P_m,m
and then iterating up to P_l,m divided by P_m,m and the scale factor 1e280.
Eventually, the result is multiplied again with these to terms.

CALLING SEQUENCE:
    plm,dplm = plm_holmes(LMAX, np.cos(theta))

INPUTS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    x: elements ranging from -1 to 1
        typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    plms: Legendre polynomials of x (geodesy normalization)
    dplms: first differentials of Legendre polynomials of x

OPTIONS:
    ASTYPE: output variable type (e.g. np.float128).  Default is np.float64

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order
    normalised associated Legendre functions" Journal of Geodesy,
    76: 279-299, 2002. https://doi.org/10.1007/s00190-002-0216-2

    Geoid Cookbook: http://mitgcm.org/~mlosch/geoidcookbook.pdf

UPDATE HISTORY:
    Updated 08/2020: prevent zero divisions by changing u==0 to eps of data type
    Updated 07/2020: added function docstrings
    Updated 10/2018: using future division for python3 Compatibility
    Updated 07/2017: output first differential of legendre polynomials
    Written 05/2015
"""
from __future__ import division
import numpy as np

def plm_holmes(LMAX, x, ASTYPE=np.float):
    """
    Computes fully-normalized associated Legendre Polynomials and
    their first derivative using Holmes and Featherstone relation

    Arguments
    ---------
    LMAX: Upper bound of Spherical Harmonic Degrees
    x: elements ranging from -1 to 1

    Keyword arguments
    -----------------
    ASTYPE: output variable data type

    Returns
    -------
    plms: fully-normalized Legendre polynomials
    dplms: first differentials of Legendre polynomials
    """

    if (np.ndim(x) > 0):
        #-- length of the x array
        jm = np.shape(x)[0]
    else:
        jm = 1

    LMAX = np.int(LMAX)
    #-- removing singleton dimensions of x
    x = np.squeeze(x).astype(ASTYPE)
    #-- scaling factor
    scalef = 1.0e-280

    #-- allocate for multiplicative factors, and plms
    f1 = np.zeros(((LMAX+1)*(LMAX+2)//2),dtype=ASTYPE)
    f2 = np.zeros(((LMAX+1)*(LMAX+2)//2),dtype=ASTYPE)
    p = np.zeros(((LMAX+1)*(LMAX+2)//2,jm),dtype=ASTYPE)
    plm = np.zeros((LMAX+1,LMAX+1,jm),dtype=ASTYPE)
    dplm = np.zeros((LMAX+1,LMAX+1,jm),dtype=ASTYPE)

    #-- Precompute multiplicative factors used in recursion relationships
    #-- Note that prefactors are not used for the case when m=l and m=l-1,
    #-- as a different recursion is used for these two values.
    k = 2#-- k = l*(l+1)/2 + m
    for l in range(2, LMAX+1):
        k += 1
        f1[k] = np.sqrt(2.0*l-1.0)*np.sqrt(2.0*l+1.0)/np.float128(l)
        f2[k] = np.float128(l-1.0)*np.sqrt(2.0*l+1.0)/(np.sqrt(2.0*l-3.0)*np.float128(l))
        for m in range(1, l-1):
            k += 1
            f1[k] = np.sqrt(2.0*l+1.0)*np.sqrt(2.0*l-1.0)/(np.sqrt(l+m)*np.sqrt(l-m))
            f2[k] = np.sqrt(2.0*l+1.0)*np.sqrt(l-m-1.0)*np.sqrt(l+m-1.0)/ \
                (np.sqrt(2.0*l-3.0)*np.sqrt(l+m)*np.sqrt(l-m))
        k += 2

    #-- u is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    #-- for x=cos(th): u=sin(th)
    u = np.sqrt(1.0 - x**2)
    #-- update where u==0 to eps of data type to prevent invalid divisions
    u[u == 0] = np.finfo(u.dtype).eps

    #-- Calculate P(l,0). These are not scaled.
    p[0,:] = 1.0
    p[1,:]  = np.sqrt(3.0)*x
    k = 1
    for l in range(2, LMAX+1):
        k += l
        p[k,:] = f1[k]*x*p[k-l,:] - f2[k]*p[k-2*l+1,:]

    #-- Calculate P(m,m), P(m+1,m), and P(l,m)
    pmm = np.sqrt(2.0)*scalef
    rescalem = 1.0/scalef
    kstart = 0

    for m in range(1, LMAX):
        rescalem = rescalem * u
        #-- Calculate P(m,m)
        kstart += m+1
        pmm = pmm * np.sqrt(2*m+1)/np.sqrt(2*m)
        p[kstart,:] = pmm
        #-- Calculate P(m+1,m)
        k = kstart+m+1
        p[k,:] = x*np.sqrt(2*m+3)*pmm
        #-- Calculate P(l,m)
        for l in range(m+2, LMAX+1):
            k += l
            p[k,:] = x*f1[k]*p[k-l,:] - f2[k]*p[k-2*l+1,:]
            p[k-2*l+1,:] = p[k-2*l+1,:] * rescalem
        #-- rescale
        p[k,:] = p[k,:] * rescalem
        p[k-LMAX,:] = p[k-LMAX,:] * rescalem

    #-- Calculate P(LMAX,LMAX)
    rescalem = rescalem * u
    kstart += m+2
    p[kstart,:] = pmm * np.sqrt(2*LMAX+1) / np.sqrt(2*LMAX) * rescalem
    #-- reshape Legendre polynomials to output dimensions
    for m in range(LMAX+1):
        for l in range(m,LMAX+1):
            lm = (l*(l+1))//2 + m
            plm[l,m,:] = p[lm,:]
            #-- calculate first derivatives
            if (l == m):
                dplm[l,m,:] = np.float128(m)*(x/u)*plm[l,m,:]
            else:
                flm = np.sqrt(((l**2.0 - m**2.0)*(2.0*l + 1.0))/(2.0*l - 1.0))
                dplm[l,m,:]= (1.0/u)*(l*x*plm[l,m,:] - flm*plm[l-1,m,:])

    #-- return the legendre polynomials and their first derivative
    return plm,dplm
