#!/usr/bin/env python
u"""
legendre_polynomials.py
Written by Tyler Sutterley (04/2022)

Computes fully normalized Legendre polynomials for an array of x values
    and their first derivative
Calculates Legendre polynomials for zonal harmonics (order 0)

CALLING SEQUENCE:
    Pl,dPl = legendre_polynomials(lmax, np.cos(theta))

INPUTS:
    lmax: maximum degree of Legrendre polynomials
    x: elements ranging from -1 to 1
        typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    pl: Legendre polynomials (geodesy normalization)
    dpl: first derivative of Legendre polynomials

OPTIONS:
    ASTYPE: output variable type.  Default is np.float64

NOTES:
    ptemp is a dummy array of length (0:lmax) storing unnormalized pl values

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 09/2020: verify dimensions of x variable
    Updated 08/2020: prevent zero divisions by changing u==0 to eps of data type
    Updated 07/2020: added function docstrings
    Updated 07/2017: added first derivative of Legendre polynomials (dpl)
        added option ASTYPE to output as different variable types e.g. np.float64
    Written 03/2013
"""
import numpy as np

def legendre_polynomials(lmax,x,ASTYPE=np.float64):
    """
    Computes fully-normalized Legendre polynomials and their first derivative
    following [HofmannWellenhof2006]_

    Parameters
    ----------
    lmax: int
        maximum degree of Legrendre polynomials
    x: float
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    ASTYPE: obj, default np.float64
        output variable data type

    Returns
    -------
    pl: float
        fully-normalized Legendre polynomials
    dpl: float
        first derivative of Legendre polynomials

    References
    ----------
    .. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz,
        *Physical Geodesy*, 2nd Edition, 403 pp., (2006).
        `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
    """
    #-- verify dimensions
    x = np.atleast_1d(x).flatten().astype(ASTYPE)
    #-- size of the x array
    nx = len(x)
    #-- verify data type of spherical harmonic truncation
    lmax = np.int64(lmax)
    #-- output matrix of normalized legendre polynomials
    pl = np.zeros((lmax+1,nx),dtype=ASTYPE)
    #-- output matrix of First derivative of Legendre polynomials
    dpl = np.zeros((lmax+1,nx),dtype=ASTYPE)
    #-- dummy matrix for the recurrence relation
    ptemp = np.zeros((lmax+1,nx),dtype=ASTYPE)

    #-- u is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    #-- for x=cos(th): u=sin(th)
    u = np.sqrt(1.0 - x**2)
    #-- update where u==0 to eps of data type to prevent invalid divisions
    u[u == 0] = np.finfo(u.dtype).eps

    #-- Initialize the recurrence relation
    ptemp[0,:] = 1.0
    ptemp[1,:] = x
    #-- Normalization is geodesy convention
    pl[0,:] = ptemp[0,:]
    pl[1,:] = np.sqrt(3.0)*ptemp[1,:]
    for l in range(2,lmax+1):
        ptemp[l,:] = (((2.0*l)-1.0)/l)*x*ptemp[l-1,:] - ((l-1.0)/l)*ptemp[l-2,:]
        #-- Normalization is geodesy convention
        pl[l,:] = np.sqrt((2.0*l)+1.0)*ptemp[l,:]

    #-- First derivative of Legendre polynomials
    for l in range(1,lmax+1):
        fl = np.sqrt(((l**2.0) * (2.0*l + 1.0)) / (2.0*l - 1.0))
        dpl[l,:] = (1.0/u)*(l*x*pl[l,:] - fl*pl[l-1,:])

    #-- return the legendre polynomials and their first derivative
    return (pl, dpl)
