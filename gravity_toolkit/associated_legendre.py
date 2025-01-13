#!/usr/bin/env python
u"""
associated_legendre.py
Written by Tyler Sutterley (03/2023)

Computes fully-normalized associated Legendre Polynomials

UPDATE HISTORY:
    Updated 03/2023: improve typing for variables in docstrings
    Updated 01/2023: refactored associated legendre polynomials
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 09/2020: verify dimensions of input x variable
    Updated 08/2020: prevent zero divisions by changing u==0 to eps of data type
    Updated 07/2020: added function docstrings
    Updated 10/2018: using future division for python3 Compatibility
    Updated 07/2017: output first differential of legendre polynomials
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Updated 09/2013: new format for file headers
    Written 03/2013
"""
from __future__ import division
import numpy as np

def associated_legendre(LMAX, x,
        method='holmes',
        MMAX=None,
        astype=np.float64
    ):
    """
    Computes fully-normalized associated Legendre Polynomials and their
    first derivative

    Parameters
    ----------
    LMAX: int
        maximum degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    method: str, default 'holmes'
        Method for computing the associated Legendre polynomials

            - ``'columbo'``
            - ``'holmes'``
            - ``'mohlenkamp'``
    MMAX: int or NoneType, default None
        maximum order of Associated Legendre polynomials
    astype: np.dtype, default np.float64
        output variable data type

    Returns
    -------
    plms: np.ndarray
        fully-normalized Legendre polynomials
    dplms: np.ndarray
        first derivative of Legendre polynomials
    """
    if (method.lower() == 'colombo'):
        return plm_colombo(LMAX, x, MMAX=MMAX, astype=astype)
    elif (method.lower() == 'holmes'):
        return plm_holmes(LMAX, x, MMAX=MMAX, astype=astype)
    elif (method.lower() == 'mohlenkamp'):
        return plm_mohlenkamp(LMAX, x, MMAX=MMAX, astype=astype)
    raise ValueError(f'Unknown method {method}')

def plm_colombo(LMAX, x,
        MMAX=None,
        astype=np.float64
    ):
    """
    Computes fully-normalized associated Legendre Polynomials and their
    first derivative using a Standard forward column method :cite:p:`Colombo:1981vh`

    Parameters
    ----------
    LMAX: int
        maximum degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    MMAX: int or NoneType, default None
        maximum order of Associated Legendre polynomials
    astype: np.dtype, default np.float64
        output variable data type

    Returns
    -------
    plm: np.ndarray
        fully-normalized Legendre polynomials
    dplm: np.ndarray
        first derivative of Legendre polynomials
    """

    # removing singleton dimensions of x
    x = np.atleast_1d(x).flatten().astype(astype)
    # length of the x array
    jm = len(x)
    # verify data type of spherical harmonic truncation
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # allocating for the plm matrix and differentials
    plm = np.zeros((LMAX+1,LMAX+1,jm))
    dplm = np.zeros((LMAX+1,LMAX+1,jm))

    # u is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    # for x=cos(th): u=sin(th)
    u = np.sqrt(1.0 - x**2)
    # update where u==0 to eps of data type to prevent invalid divisions
    u[u == 0] = np.finfo(u.dtype).eps

    # Calculating the initial polynomials for the recursion
    plm[0,0,:] = 1.0
    plm[1,0,:] = np.sqrt(3.0)*x
    plm[1,1,:] = np.sqrt(3.0)*u
    # calculating first derivatives for harmonics of degree 1
    dplm[1,0,:] = (1.0/u)*(x*plm[1,0,:] - np.sqrt(3)*plm[0,0,:])
    dplm[1,1,:] = (x/u)*plm[1,1,:]
    for l in range(2, LMAX+1):
        for m in range(0, l):# Zonal and Tesseral harmonics (non-sectorial)
            # Computes the non-sectorial terms from previously computed
            # sectorial terms.
            alm = np.sqrt(((2.0*l-1.0)*(2.0*l+1.0))/((l-m)*(l+m)))
            blm = np.sqrt(((2.0*l+1.0)*(l+m-1.0)*(l-m-1.0))/((l-m)*(l+m)*(2.0*l-3.0)))
            # if (m == l-1): plm[l-2,m,:] will be 0
            plm[l,m,:] = alm*x*plm[l-1,m,:] - blm*plm[l-2,m,:]
            # calculate first derivatives
            flm = np.sqrt(((l**2.0 - m**2.0)*(2.0*l + 1.0))/(2.0*l - 1.0))
            dplm[l,m,:] = (1.0/u)*(l*x*plm[l,m,:] - flm*plm[l-1,m,:])

        # Sectorial harmonics
        # The sectorial harmonics serve as seed values for the recursion
        # starting with P00 and P11 (outside the loop)
        plm[l,l,:] = u*np.sqrt((2.0*l+1.0)/(2.0*l))*np.squeeze(plm[l-1,l-1,:])
        # calculate first derivatives for sectorial harmonics
        dplm[l,l,:] = np.longdouble(l)*(x/u)*plm[l,l,:]

    # return the legendre polynomials and their first derivative
    # truncating orders to MMAX
    return plm[:,:MMAX+1,:], dplm[:,:MMAX+1,:]

def plm_holmes(LMAX, x,
        MMAX=None,
        astype=np.float64
    ):
    """
    Computes fully-normalized associated Legendre Polynomials and their
    first derivative using the recursion relation from :cite:p:`Holmes:2002ff`

    Parameters
    ----------
    LMAX: int
        maximum degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    MMAX: int or NoneType, default None
        maximum order of Associated Legendre polynomials
    astype: np.dtype, default np.float64
        output variable data type

    Returns
    -------
    plm: np.ndarray
        fully-normalized Legendre polynomials
    dplm: np.ndarray
        first derivative of Legendre polynomials
    """

    # removing singleton dimensions of x
    x = np.atleast_1d(x).flatten().astype(astype)
    # length of the x array
    jm = len(x)
    # verify data type of spherical harmonic truncation
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)
    # scaling factor
    scalef = 1.0e-280

    # allocate for multiplicative factors, and plms
    f1 = np.zeros(((LMAX+1)*(LMAX+2)//2), dtype=astype)
    f2 = np.zeros(((LMAX+1)*(LMAX+2)//2), dtype=astype)
    p = np.zeros(((LMAX+1)*(LMAX+2)//2,jm), dtype=astype)
    plm = np.zeros((LMAX+1,LMAX+1,jm), dtype=astype)
    dplm = np.zeros((LMAX+1,LMAX+1,jm), dtype=astype)

    # Precompute multiplicative factors used in recursion relationships
    # Note that prefactors are not used for the case when m=l and m=l-1,
    # as a different recursion is used for these two values.
    k = 2# k = l*(l+1)/2 + m
    for l in range(2, LMAX+1):
        k += 1
        f1[k] = np.sqrt(2.0*l-1.0)*np.sqrt(2.0*l+1.0)/np.longdouble(l)
        f2[k] = np.longdouble(l-1.0)*np.sqrt(2.0*l+1.0)/(np.sqrt(2.0*l-3.0)*np.longdouble(l))
        for m in range(1, l-1):
            k += 1
            f1[k] = np.sqrt(2.0*l+1.0)*np.sqrt(2.0*l-1.0)/(np.sqrt(l+m)*np.sqrt(l-m))
            f2[k] = np.sqrt(2.0*l+1.0)*np.sqrt(l-m-1.0)*np.sqrt(l+m-1.0)/ \
                (np.sqrt(2.0*l-3.0)*np.sqrt(l+m)*np.sqrt(l-m))
        k += 2

    # u is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    # for x=cos(th): u=sin(th)
    u = np.sqrt(1.0 - x**2)
    # update where u==0 to eps of data type to prevent invalid divisions
    u[u == 0] = np.finfo(u.dtype).eps

    # Calculate P(l,0). These are not scaled.
    p[0,:] = 1.0
    p[1,:]  = np.sqrt(3.0)*x
    k = 1
    for l in range(2, LMAX+1):
        k += l
        p[k,:] = f1[k]*x*p[k-l,:] - f2[k]*p[k-2*l+1,:]

    # Calculate P(m,m), P(m+1,m), and P(l,m)
    pmm = np.sqrt(2.0)*scalef
    rescalem = 1.0/scalef
    kstart = 0

    for m in range(1, LMAX):
        rescalem = rescalem * u
        # Calculate P(m,m)
        kstart += m+1
        pmm = pmm * np.sqrt(2*m+1)/np.sqrt(2*m)
        p[kstart,:] = pmm
        # Calculate P(m+1,m)
        k = kstart+m+1
        p[k,:] = x*np.sqrt(2*m+3)*pmm
        # Calculate P(l,m)
        for l in range(m+2, LMAX+1):
            k += l
            p[k,:] = x*f1[k]*p[k-l,:] - f2[k]*p[k-2*l+1,:]
            p[k-2*l+1,:] = p[k-2*l+1,:] * rescalem
        # rescale
        p[k,:] = p[k,:] * rescalem
        p[k-LMAX,:] = p[k-LMAX,:] * rescalem

    # Calculate P(LMAX,LMAX)
    rescalem = rescalem * u
    kstart += m+2
    p[kstart,:] = pmm * np.sqrt(2*LMAX+1) / np.sqrt(2*LMAX) * rescalem
    # reshape Legendre polynomials to output dimensions
    for m in range(LMAX+1):
        for l in range(m,LMAX+1):
            lm = (l*(l+1))//2 + m
            plm[l,m,:] = p[lm,:]
            # calculate first derivatives
            if (l == m):
                dplm[l,m,:] = np.longdouble(m)*(x/u)*plm[l,m,:]
            else:
                flm = np.sqrt(((l**2.0 - m**2.0)*(2.0*l + 1.0))/(2.0*l - 1.0))
                dplm[l,m,:]= (1.0/u)*(l*x*plm[l,m,:] - flm*plm[l-1,m,:])

    # return the legendre polynomials and their first derivative
    # truncating orders to MMAX
    return plm[:,:MMAX+1,:], dplm[:,:MMAX+1,:]

def plm_mohlenkamp(LMAX, x,
        MMAX=None,
        astype=np.float64
    ):
    """
    Computes fully-normalized associated Legendre Polynomials and their
    first derivative using the recursion relation from :cite:p:`Mohlenkamp:2016vv`

    Derived from :cite:p:`Szego:1939tn` recurrence formula for Jacobi Polynomials

    Parameters
    ----------
    LMAX: int
        maximum degree of Legendre polynomials
    x: np.ndarray
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    MMAX: int or NoneType, default None
        maximum order of Associated Legendre polynomials
    astype: np.dtype, default np.float64
        output variable data type

    Returns
    -------
    plm: np.ndarray
        fully-normalized Legendre polynomials
    dplm: np.ndarray
        first derivative of Legendre polynomials
    """

    # Verify LMAX as integer
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # removing singleton dimensions of x
    x = np.atleast_1d(x).flatten()
    # length of the x array
    sx = len(x)

    # Initialize the output Legendre polynomials
    plm = np.zeros((LMAX+1, MMAX+1, sx), dtype=astype)
    dplm = np.zeros((LMAX+1, LMAX+1, sx), dtype=astype)
    # Jacobi polynomial for the recurrence relation
    jlmm = np.zeros((LMAX+1, MMAX+1, sx))
    # for x=cos(th): u= sin(th)
    u = np.sqrt(1.0 - x**2)
    # update where u==0 to eps of data type to prevent invalid divisions
    u[u == 0] = np.finfo(u.dtype).eps

    # for all spherical harmonic orders of interest
    for mm in range(0,MMAX+1):# equivalent to 0:MMAX
        # Initialize the recurrence relation
        # J-1,m,m Term == 0
        # J0,m,m Term
        if (mm > 0):
            # j ranges from 1 to mm for the product
            j = np.arange(0,mm)+1.0
            jlmm[0,mm,:] = np.prod(np.sqrt(1.0 + 1.0/(2.0*j)))/np.sqrt(2.0)
        else: # if mm == 0: jlmm = 1/sqrt(2)
            jlmm[0,mm,:] = 1.0/np.sqrt(2.0)
        # Jk,m,m Terms
        for k in range(1, LMAX+1):# computation for SH degrees
            # Initialization begins at -1
            # this is to make the formula parallel the function written in
            # Martin Mohlenkamp's Guide to Spherical Harmonics
            # Jacobi General Terms
            if (k == 1):# for degree 1 terms
                jlmm[k,mm,:] = 2.0*x * jlmm[k-1,mm,:] * \
                    np.sqrt(1.0 + (mm - 0.5)/k) * \
                    np.sqrt(1.0 - (mm - 0.5)/(k + 2.0*mm))
            else:# for all other spherical harmonic degrees
                jlmm[k,mm,:] = 2.0*x * jlmm[k-1,mm,:] * \
                    np.sqrt(1.0 + (mm - 0.5)/k) * \
                    np.sqrt(1.0 - (mm - 0.5)/(k + 2.0*mm)) - \
                    jlmm[k-2,mm,:] * np.sqrt(1.0 + 4.0/(2.0*k + 2.0*mm - 3.0)) * \
                    np.sqrt(1.0 - (1.0/k)) * np.sqrt(1.0 - 1.0/(k + 2.0*mm))
        # Normalization is geodesy convention
        for l in range(mm,LMAX+1): # equivalent to mm:LMAX
            if (mm == 0):# Geodesy normalization (m=0) == sqrt(2)*sin(th)^0
                # u^mm term is dropped as u^0 = 1
                plm[l,mm,:] = np.sqrt(2.0)*jlmm[l-mm,mm,:]
            else:# Geodesy normalization all others == 2*sin(th)^mm
                plm[l,mm,:] = 2.0*(u**mm)*jlmm[l-mm,mm,:]
            # calculate first derivatives
            if (l == mm):
                dplm[l,mm,:] = np.longdouble(mm)*(x/u)*plm[l,mm,:]
            else:
                flm = np.sqrt(((l**2.0 - mm**2.0)*(2.0*l + 1.0))/(2.0*l - 1.0))
                dplm[l,mm,:]= (1.0/u)*(l*x*plm[l,mm,:] - flm*plm[l-1,mm,:])
    # return the legendre polynomials and their first derivative
    return plm, dplm
