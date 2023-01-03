#!/usr/bin/env python
u"""
plm_mohlenkamp.py
Written by Tyler Sutterley (04/2022)

Computes fully-normalized associated Legendre Polynomials
    for an array of x values
Uses Martin Mohlenkamp's recursion relation derived from the
    Szego (1939) Recurrence formula for Jacobi Polynomials (Pg 71)

With this algorithm, the associated Legendre Functions are
    constructed as an amplitude times a Jacobi Polynomial
    P[l,m](cos(theta)) = (sin(theta)^2)*J[l-m,m,m](cos(theta))

CALLING SEQUENCE:
    plm = plm_mohlenkamp(LMAX, np.cos(theta))

INPUTS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    x: elements ranging from -1 to 1
        typically cos(theta), where theta is the colatitude in radians

OUTPUT:
    plm: Legendre polynomials (geodesy normalization)

OPTIONS:
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

NOTES:
    Modified and updated from IDL plm_x.pro coded by Sean Swenson
    Difference from martin function in geoid_mk.mac.f:
        plms from plm_mohlenkamp are normalized inside the function
        plms from martin are normalized outside the function
    For large spherical harmonic degrees this recurrence relation
        is poorly conditioned
    For spherical harmonic orders above ~1000 can cause overflows

REFERENCES:
    Martin Mohlenkamp, "A User's Guide to Spherical Harmonics"
    http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 09/2020: verify dimensions of input x variable
    Updated 07/2020: added function docstrings
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 09/2013
"""
import warnings
import numpy as np

def plm_mohlenkamp(LMAX, x, MMAX=None):
    """
    Computes fully-normalized associated Legendre Polynomials
    using Martin Mohlenkamp's recursion relation [Mohlenkamp2016]_

    Derived from [Szego1939]_ recurrence formula for Jacobi Polynomials

    Parameters
    ----------
    LMAX: int
        maximum degree of Legrendre polynomials
    x: float
        elements ranging from -1 to 1

        Typically ``cos(theta)``, where ``theta`` is the colatitude in radians
    MMAX: int or NoneType, default None
        maximum order of Associated Legrendre polynomials

    Returns
    -------
    plms: float
        fully-normalized Legendre polynomials
    dplms: float
        first derivative of Legendre polynomials

    References
    ----------
    .. [Mohlenkamp2016] M. J. Mohlenkamp,
        "A User's Guide to Spherical Harmonics", (2016).
        `[pdf] <http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf>`_
    .. [Szego1939] Gabor Szeg\ |ouml|\ , "Orthogonal Polynomials", 440 pp., (1939).
        `[pdf] <https://people.math.osu.edu/nevai.1/AT/SZEGO/szego=szego1975=ops=OCR.pdf>`_

    .. |ouml|    unicode:: U+00F6 .. LATIN SMALL LETTER O WITH DIAERESIS
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.associated_legendre instead",
        DeprecationWarning)

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
    plm=np.zeros((LMAX+1,MMAX+1,sx))
    # Jacobi polynomial for the recurrence relation
    jlmm=np.zeros((LMAX+1,MMAX+1,sx))
    # for x=cos(th): rsin= sin(th)
    rsin=np.sqrt(1.0 - x**2)

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
                # rsin^mm term is dropped as rsin^0 = 1
                plm[l,mm,:] = np.sqrt(2.0)*jlmm[l-mm,mm,:]
            else:# Geodesy normalization all others == 2*sin(th)^mm
                plm[l,mm,:] = 2.0*(rsin**mm)*jlmm[l-mm,mm,:]
    return plm
