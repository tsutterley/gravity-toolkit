#!/usr/bin/env python
u"""
destripe_harmonics.py
Original Fortran program remove_errors.f written by Isabella Velicogna
Adapted by Chia-Wei Hsu (05/2018)
Updated by Tyler Sutterley (04/2022)

Filters spherical harmonic coefficients for correlated "striping" errors

ALGORITHM:
    clm1, slm1 are the spherical harmonic coefficients
        after the mean field has been removed
    Smooth values over l, for l=even and l=odd separately
        by fitting a quadratic function to every 7 points
    Remove those smoothed values

CALLING SEQUENCE:
    Ylms = destripe_harmonics(clm,slm,LMAX=60)
    Wclm = WYlms['clm']
    Wslm = WYlms['slm']

INPUTS:
    clm1: cosine spherical harmonic coefficients (matrix 2 dims)
    slm1: sine spherical harmonic coefficients (matrix 2 dims)
        clm1 and slm1 are matrix with 2 dimensions
        the dimensions are in the following order [l,m]

OUTPUTS:
    Wclm: filtered cosine spherical harmonic coefficients
    Wslm: filtered sine spherical harmonic coefficients

OPTIONS:
    LMIN: Lower bound of Spherical Harmonic Degrees (default = 2)
    LMAX: Upper bound of Spherical Harmonic Degrees (default = 60)
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    ROUND: use round to find nearest even (True) or use floor (False)
    NARROW: Clm=Slm=0 if number of points is less than window size (False)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    S C Swenson and J Wahr, "Post-processing removal of correlated errors in
        GRACE data", Geophysical Research Letters, 33(L08402), 2006
        https://doi.org/10.1029/2005GL025285

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 07/2020: added function docstrings
    Updated 03/2020: Updated for public release
    Updated 05/2018: using __future__ print and updated flags comments
    Updated 08/2015: changed from sys.exit to raise ValueError
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Updated 02/2014: generalization for GRACE GUI and other routines
"""
from __future__ import print_function
import numpy as np

def destripe_harmonics(clm1, slm1, LMIN=2, LMAX=60, MMAX=None,
    ROUND=True, NARROW=False):
    """
    Filters spherical harmonic coefficients for correlated striping errors

    Parameters
    ----------
    clm1: float
        cosine spherical harmonic coefficients
    slm1: float
        sine spherical harmonic coefficients
    LMIN: int, default 2
        Lower bound of Spherical Harmonic Degrees
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    ROUND: bool, default True
        use round to find nearest even
    NARROW: bool, default True
        set harmonics to 0 if less than window size

    Returns
    -------
    clm: float
        filtered cosine spherical harmonic coefficients
    slm: float
        filtered sine spherical harmonic coefficients

    References
    ----------
    .. [Swenson2006] S. Swenson and J. Wahr,
        "Post-processing removal of correlated errors in GRACE data",
        *Geophysical Research Letters*, 33(L08402), (2006).
        `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_
    """

    # tests if spherical harmonics have been imported
    if (clm1.shape[0] == 1) or (slm1.shape[0] == 1):
        raise ValueError('Input harmonics need to be matrices')

    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # output filtered coefficients (copy to not modify input)
    Wclm = clm1.copy()
    Wslm = slm1.copy()
    # matrix size declarations
    clmeven = np.zeros((LMAX), dtype=np.float64)
    slmeven = np.zeros((LMAX), dtype=np.float64)
    clmodd = np.zeros((LMAX+1), dtype=np.float64)
    slmodd = np.zeros((LMAX+1), dtype=np.float64)
    clmsm = np.zeros((LMAX+1,MMAX+1), dtype=np.float64)
    slmsm = np.zeros((LMAX+1,MMAX+1), dtype=np.float64)

    # start of the smoothing over orders (m)
    for m in range(int(MMAX+1)):
        smooth = np.exp(-np.float64(m)/10.0)*15.0
        if ROUND:
            # round(smooth) to nearest even instead of int(smooth)
            nsmooth = np.around(smooth)
        else:
            # Sean's method for finding nsmooth (use floor of smooth)
            nsmooth = np.int64(smooth)

        if (nsmooth < 2):
            # Isabella's method of picking nsmooth sets minimum to 2
            nsmooth = np.int64(2)

        rmat = np.zeros((3,3), dtype=np.float64)
        lll = np.arange(np.float64(nsmooth)*2.+1.)-np.float64(nsmooth)
        # create design matrix to have the following form:
        #    [    1     ll     ll^2   ]
        #    [    ll    ll^2   ll^3   ]
        #    [    ll^2  ll^3   ll^4   ]
        for i,ill in enumerate(lll):
            rmat[0,0] += 1.0
            rmat[0,1] += ill
            rmat[0,2] += ill**2

            rmat[1,0] += ill
            rmat[1,1] += ill**2
            rmat[1,2] += ill**3

            rmat[2,0] += ill**2
            rmat[2,1] += ill**3
            rmat[2,2] += ill**4

        # put the even and odd l's into their own arrays
        ieven = -1
        iodd = -1
        leven = np.zeros((LMAX), dtype=np.int64)
        lodd = np.zeros((LMAX), dtype=np.int64)

        for l in range(int(m),int(LMAX+1)):
            # check if degree is odd or even
            if np.remainder(l,2).astype(bool):
                iodd += 1
                lodd[iodd] = l
                clmodd[iodd] = clm1[l,m].copy()
                slmodd[iodd] = slm1[l,m].copy()
            else:
                ieven += 1
                leven[ieven] = l
                clmeven[ieven] = clm1[l,m].copy()
                slmeven[ieven] = slm1[l,m].copy()

        # smooth, by fitting a quadratic polynomial to 7 points at a time
        # deal with even stokes coefficients
        l1 = 0
        l2 = ieven

        if (l1 > (l2-2*nsmooth)):
            for l in range(l1,l2+1):
                if NARROW:
                    # Sean's method
                    # Clm=Slm=0 if number of points is less than window size
                    clmsm[leven[l],m] = 0.0
                    slmsm[leven[l],m] = 0.0
                else:
                    # Isabella's method
                    # Clm and Slm passed through unaltered
                    clmsm[leven[l],m] = clm1[leven[l],m].copy()
                    slmsm[leven[l],m] = slm1[leven[l],m].copy()
        else:
            for l in range(int(l1+nsmooth),int(l2-nsmooth+1)):
                rhsc = np.zeros((3), dtype=np.float64)
                rhss = np.zeros((3), dtype=np.float64)
                for ll in range(int(-nsmooth),int(nsmooth+1)):
                    rhsc[0] += clmeven[l+ll]
                    rhsc[1] += clmeven[l+ll]*np.float64(ll)
                    rhsc[2] += clmeven[l+ll]*np.float64(ll**2)
                    rhss[0] += slmeven[l+ll]
                    rhss[1] += slmeven[l+ll]*np.float64(ll)
                    rhss[2] += slmeven[l+ll]*np.float64(ll**2)

                # fit design matrix to coefficients
                # to get beta parameters
                bhsc = np.linalg.lstsq(rmat,rhsc.T,rcond=-1)[0]
                bhss = np.linalg.lstsq(rmat,rhss.T,rcond=-1)[0]

                # all other l is assigned as bhsc
                clmsm[leven[l],m] = bhsc[0].copy()
                # all other l is assigned as bhss
                slmsm[leven[l],m] = bhss[0].copy()

                if (l == (l1+nsmooth)):
                    # deal with l=l1+nsmooth
                    for ll in range(int(-nsmooth),0):
                        clmsm[leven[l+ll],m] = bhsc[0]+bhsc[1]*np.float64(ll) + \
                            bhsc[2]*np.float64(ll**2)
                        slmsm[leven[l+ll],m] = bhss[0]+bhss[1]*np.float64(ll) + \
                            bhss[2]*np.float64(ll**2)

                if (l == (l2-nsmooth)):
                    # deal with l=l2-nsmnooth
                    for ll in range(1,int(nsmooth+1)):
                        clmsm[leven[l+ll],m] = bhsc[0]+bhsc[1]*np.float64(ll) + \
                            bhsc[2]*np.float64(ll**2)
                        slmsm[leven[l+ll],m] = bhss[0]+bhss[1]*np.float64(ll) + \
                            bhss[2]*np.float64(ll**2)

        # deal with odd stokes coefficients
        l1 = 0
        l2 = iodd

        if (l1 > (l2-2*nsmooth)):
            for l in range(l1,l2+1):
                if NARROW:
                    # Sean's method
                    # Clm=Slm=0 if number of points is less than window size
                    clmsm[lodd[l],m] = 0.0
                    slmsm[lodd[l],m] = 0.0
                else:
                    # Isabella's method
                    # Clm and Slm passed through unaltered
                    clmsm[lodd[l],m] = clm1[lodd[l],m].copy()
                    slmsm[lodd[l],m] = slm1[lodd[l],m].copy()
        else:
            for l in range(int(l1+nsmooth),int(l2-nsmooth+1)):
                rhsc = np.zeros((3), dtype=np.float64)
                rhss = np.zeros((3), dtype=np.float64)
                for ll in range(int(-nsmooth),int(nsmooth+1)):
                    rhsc[0] += clmodd[l+ll]
                    rhsc[1] += clmodd[l+ll]*np.float64(ll)
                    rhsc[2] += clmodd[l+ll]*np.float64(ll**2)
                    rhss[0] += slmodd[l+ll]
                    rhss[1] += slmodd[l+ll]*np.float64(ll)
                    rhss[2] += slmodd[l+ll]*np.float64(ll**2)

                # fit design matrix to coefficients
                # to get beta parameters
                bhsc = np.linalg.lstsq(rmat,rhsc.T,rcond=-1)[0]
                bhss = np.linalg.lstsq(rmat,rhss.T,rcond=-1)[0]

                # all other l is assigned as bhsc
                clmsm[lodd[l],m] = bhsc[0].copy()
                # all other l is assigned as bhss
                slmsm[lodd[l],m] = bhss[0].copy()

                if (l == (l1+nsmooth)):
                    # deal with l=l1+nsmooth
                    for ll in range(int(-nsmooth),0):
                        clmsm[lodd[l+ll],m] = bhsc[0]+bhsc[1]*np.float64(ll) + \
                            bhsc[2]*np.float64(ll**2)
                        slmsm[lodd[l+ll],m] = bhss[0]+bhss[1]*np.float64(ll) + \
                            bhss[2]*np.float64(ll**2)

                if (l == (l2-nsmooth)):
                    # deal with l=l2-nsmnooth
                    for ll in range(1,int(nsmooth+1)):
                        clmsm[lodd[l+ll],m] = bhsc[0]+bhsc[1]*np.float64(ll) + \
                            bhsc[2]*np.float64(ll**2)
                        slmsm[lodd[l+ll],m] = bhss[0]+bhss[1]*np.float64(ll) + \
                            bhss[2]*np.float64(ll**2)

        # deal with m greater than or equal to 5
        for l in range(int(m),int(LMAX+1)):
            if (m >= 5):
                # remove smoothed clm/slm from original spherical harmonics
                Wclm[l,m] -= clmsm[l,m]
                Wslm[l,m] -= slmsm[l,m]

    return {'clm':Wclm,'slm':Wslm}
