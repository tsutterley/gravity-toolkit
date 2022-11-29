#!/usr/bin/env python
u"""
gauss_weights.py
Original IDL code gauss_weights.pro written by Sean Swenson
Adapted by Tyler Sutterley (04/2022)

Computes the Gaussian weights as a function of degree
A normalized version of Jekeli's Gaussian averaging function

Christopher Jekeli (1981)
Alternative Methods to Smooth the Earth's Gravity Field
http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf

CALLING SEQUENCE:
    wl = gauss_weights(hw, LMAX)

INPUTS:
    hw: Gaussian smoothing radius in kilometers
        Radius r corresponds to the distance at which the weight
        drops to half its peak value at the shortest wavelength
    LMAX: Maximum degree of spherical harmonic coefficients

OPTIONS:
    CUTOFF: minimum value for tail of Gaussian averaging function

OUTPUTS:
    wl: degree dependent weighting function

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

NOTES:
    Differences from recurs function in combine.mac.f:
        weighting from gauss_weights is normalized outside of the function
            wt = 2.0*pi*gauss_weights(rad,LMAX)
        weighting from recurs is normalized inside of the function
        call recurs(alpha,bcoef) calculates bcoef up to LMAX 150 (=wt[0:150])
            alpha = alog(2.)/(1.-cos(rad/6371.))

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: added option for setting minimum value threshold
    Updated 07/2020: added function docstrings
    Updated 06/2015: adjusted threshold from 1e-9 to 1e-10
    Updated 12/2014: updated comments and header text updating full reference
    Updated 02/2014: changed variables from ints to floats to prevent truncation
    Written 05/2013
"""
import numpy as np

def gauss_weights(hw, LMAX, CUTOFF=1e-10):
    """
    Computes the Gaussian weights as a function of degree using
    a normalized form from [Jekeli1981]_

    Parameters
    ----------
    hw: float
        Gaussian smoothing radius in kilometers
    LMAX: int
        Maximum degree of spherical harmonic coefficients
    CUTOFF: float, default 1e-10
        minimum value for tail of Gaussian averaging function

    Returns
    -------
    wl: float
        degree dependent weighting function

    References
    ----------
    .. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth
        the Earth's Gravity Field", NASA Grant No. NGR 36-008-161,
        OSURF Proj. No. 783210, 48 pp., (1981).
    """
    # allocate for output weights
    wl = np.zeros((LMAX+1))
    # radius of the Earth in km
    rad_e = 6371.0
    if (hw < CUTOFF):
        # distance is smaller than cutoff
        wl[:]=1.0/(2.0*np.pi)
    else:
        # calculate gaussian weights using recursion
        b = np.log(2.0)/(1.0 - np.cos(hw/rad_e))
        # weight for degree 0
        wl[0] = 1.0/(2.0*np.pi)
        # weight for degree 1
        wl[1] = wl[0]*((1.0+np.exp(-2.0*b))/(1.0-np.exp(-2.0*b))-1.0/b)
        # valid flag
        valid = True
        # spherical harmonic degree
        l = 2
        # while valid (within cutoff)
        # and spherical harmonic degree is less than LMAX
        while (valid and (l <= LMAX)):
            # calculate weight with recursion
            wl[l] = (1.0-2.0*l)/b*wl[l-1]+wl[l-2]
            # weight is less than cutoff
            if (wl[l] < CUTOFF):
                # set all weights to cutoff
                wl[l:LMAX+1] = CUTOFF
                # set valid flag
                valid = False
            # add 1 to l
            l += 1
    # return the gaussian weights
    return wl
