#!/usr/bin/env python
u"""
legendre.py
Written by Tyler Sutterley (04/2022)
Computes associated Legendre functions of degree l evaluated for elements x
l must be a scalar integer and x must contain real values ranging -1 <= x <= 1
Parallels the MATLAB legendre function

Based on Fortran program by Robert L. Parker, Scripps Institution of
Oceanography, Institute for Geophysics and Planetary Physics, UCSD. 1993

INPUTS:
    l: degree of Legrendre polynomials
    x: elements ranging from -1 to 1
        typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    Pl: legendre polynomials of degree l for orders 0 to l

OPTIONS:
    NORMALIZE: output Fully Normalized Associated Legendre Functions

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

REFERENCES:
    M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
        Dover Publications, 1965, Ch. 8.
    J. A. Jacobs, "Geomagnetism", Academic Press, 1987, Ch.4.

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: modify normalization to prevent high degree overflows
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: modify case with underflow
    Updated 09/2020: verify dimensions of x variable
    Updated 07/2020: added function docstrings
    Updated 05/2020: added normalization option for output polynomials
    Updated 03/2019: calculate twocot separately to avoid divide warning
    Written 08/2016
"""
import numpy as np

def legendre(l, x, NORMALIZE=False):
    """
    Computes associated Legendre functions of degree ``l``
    following [Abramowitz1965]_ and [Jacobs1987]_

    Parameters
    ----------
    l: int
        degree of Legrendre polynomials
    x: float
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    NORMALIZE: bool, default False
        Fully-normalize the Legendre Functions

    Returns
    -------
    Pl: legendre polynomials of degree ``l``

    References
    ----------
    .. [Abramowitz1965] M. Abramowitz and I. A. Stegun,
        *Handbook of Mathematical Functions*, 1046 pp., (1965).

    .. [Jacobs1987] J. A. Jacobs, *Geomagnetism*,
        Volume 1, 1st Edition, 832 pp., (1987).
    """
    #-- verify integer
    l = np.int64(l)
    #-- verify dimensions
    x = np.atleast_1d(x).flatten()
    #-- size of the x array
    nx = len(x)

    #-- for the l = 0 case
    if (l == 0):
        Pl = np.ones((1,nx), dtype=np.float64)
        return Pl

    #-- for all other degrees greater than 0
    rootl = np.sqrt(np.arange(0,2*l+1))#-- +1 to include 2*l
    #-- s is sine of colatitude (cosine of latitude) so that 0 <= s <= 1
    s = np.sqrt(1.0 - x**2)#-- for x=cos(th): s=sin(th)
    P = np.zeros((l+3,nx), dtype=np.float64)

    #-- Find values of x,s for which there will be underflow
    sn = (-s)**l
    tol = np.sqrt(np.finfo(np.float64).tiny)
    count = np.count_nonzero((s > 0) & (np.abs(sn) <= tol))
    if (count > 0):
        ind, = np.nonzero((s > 0) & (np.abs(sn) <= tol))
        #-- Approximate solution of x*ln(x) = Pl
        v = 9.2 - np.log(tol)/(l*s[ind])
        w = 1.0/np.log(v)
        m1 = 1+l*s[ind]*v*w*(1.0058+ w*(3.819 - w*12.173))
        m1 = np.where(l < np.floor(m1), l, np.floor(m1)).astype(np.int64)
        #-- Column-by-column recursion
        for k,mm1 in enumerate(m1):
            col = ind[k]
            #-- Calculate twocot for underflow case
            twocot = -2.0*x[col]/s[col]
            P[mm1-1:l+1,col] = 0.0
            #-- Start recursion with proper sign
            tstart = np.finfo(np.float64).eps
            P[mm1-1,col] = np.sign(np.fmod(mm1,2)-0.5)*tstart
            if (x[col] < 0):
                P[mm1-1,col] = np.sign(np.fmod(l+1,2)-0.5)*tstart
            #-- Recur from m1 to m = 0, accumulating normalizing factor.
            sumsq = tol.copy()
            for m in range(mm1-2,-1,-1):
                P[m,col] = ((m+1)*twocot*P[m+1,col] - \
                    rootl[l+m+2]*rootl[l-m-1]*P[m+2,col]) / \
                    (rootl[l+m+1]*rootl[l-m])
                sumsq += P[m,col]**2
            #-- calculate scale
            scale = 1.0/np.sqrt(2.0*sumsq - P[0,col]**2)
            P[0:mm1+1,col] = scale*P[0:mm1+1,col]

    #-- Find the values of x,s for which there is no underflow, and (x != +/-1)
    count = np.count_nonzero((x != 1) & (np.abs(sn) >= tol))
    if (count > 0):
        nind, = np.nonzero((x != 1) & (np.abs(sn) >= tol))
        #-- Calculate twocot for normal case
        twocot = -2.0*x[nind]/s[nind]
        #-- Produce normalization constant for the m = l function
        d = np.arange(2,2*l+2,2)
        c = np.prod(1.0 - 1.0/d)
        #-- Use sn = (-s)**l (written above) to write the m = l function
        P[l,nind] = np.sqrt(c)*sn[nind]
        P[l-1,nind] = P[l,nind]*twocot*l/rootl[-1]

        #-- Recur downwards to m = 0
        for m in range(l-2,-1,-1):
            P[m,nind] = (P[m+1,nind]*twocot*(m+1) - \
                P[m+2,nind]*rootl[l+m+2]*rootl[l-m-1]) / \
                (rootl[l+m+1]*rootl[l-m])

    #-- calculate Pl from P
    Pl = np.copy(P[0:l+1,:])

    #-- Polar argument (x == +/-1)
    count = np.count_nonzero(s == 0)
    if (count > 0):
        s0, = np.nonzero(s == 0)
        Pl[0,s0] = x[s0]**l

    #-- calculate Fully Normalized Associated Legendre functions
    if NORMALIZE:
        norm = np.zeros((l+1))
        norm[0] = np.sqrt(2.0*l+1)
        m = np.arange(1,l+1)
        norm[1:] = (-1)**m*np.sqrt(2.0*(2.0*l+1.0))
        Pl *= np.kron(np.ones((1,nx)), norm[:,np.newaxis])
    else:
        #-- Calculate the unnormalized Legendre functions by multiplying each row
        #-- by: sqrt((l+m)!/(l-m)!) == sqrt(prod(n-m+1:n+m))
        #-- following Abramowitz and Stegun
        for m in range(1,l):
            Pl[m,:] *= np.prod(rootl[l-m+1:l+m+1])
        #-- sectoral case (l = m) should be done separately to handle 0!
        Pl[l,:] *= np.prod(rootl[1:])

    return Pl
