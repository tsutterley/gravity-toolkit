#!/usr/bin/env python
u"""
gen_averaging_kernel.py
Original IDL code gen_wclms_me.pro written by Sean Swenson
Adapted by Tyler Sutterley (06/2023)

Generates averaging kernel coefficients which minimize the total error

CALLING SEQUENCE:
    Wlms = gen_averaging_kernel(gclm,gslm,eclm,eslm,sigma,hw,
        LMIN=0, LMAX=60, UNITS=0, LOVE=(hl,kl,ll))

INPUTS:
    gclm: cosine spherical harmonics of exact averaging kernel
    gslm: sine spherical harmonics of exact averaging kernel
    eclm: measurement error in the cosine harmonics
    eslm: measurement error in the sine harmonics
    sigma: variance of the surface mass signal
    hw: Gaussian radius of the kernel in kilometers

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    UNITS: units of input spherical harmonics
        0: fully-normalized
        1: mass coefficients (cmwe)
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

OUTPUTS:
    clm: cosine coefficients of the averaging kernel
    slm: sine coefficients of the averaging kernel

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    units.py: class for converting spherical harmonic data to specific units

REFERENCES:
    Swenson and Wahr, "Methods for inferring regional surface-mass anomalies
        from Gravity Recovery and Climate Experiment (GRACE) measurements of
        time-variable gravity," Journal of Geophysical Research: Solid Earth,
        107(B9), (2002). https://doi.org/10.1029/2001JB000576

UPDATE HISTORY:
    Updated 06/2023: added option for setting minimum value threshold
        use harmonics class for spherical harmonic operations
    Updated 04/2023: allow love numbers to be None for mass units case
    Updated 03/2023: improve typing for variables in docstrings
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2021: using units module for Earth parameters
    Updated 04/2020: reading load love numbers outside of this function
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 05/2013
"""
import numpy as np
import gravity_toolkit.units

def gen_averaging_kernel(gclm, gslm, eclm, eslm, sigma, hw,
    LMAX=60, MMAX=None, CUTOFF=1e-15, UNITS=0, LOVE=None):
    r"""
    Generates averaging kernel coefficients which minimize the
    total error following :cite:t:`Swenson:2002hs`

    Uses a normalized form of the Gaussian averaging function
    from :cite:p:`Jekeli:1981vj`

    Parameters
    ----------
    gclm: np.ndarray
        cosine spherical harmonics of exact averaging kernel
    gslm: np.ndarray
        sine spherical harmonics of exact averaging kernel
    eclm: np.ndarray
        measurement error in the cosine harmonics
    eslm: np.ndarray
        measurement error in the sine harmonics
    sigma: float
        variance of the surface mass signal
    hw: float
        Gaussian radius of the kernel in kilometers
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    CUTOFF: float, default 1e-15
        minimum value for tail of Gaussian averaging function
    UNITS: int, default 0
        Input data units

            - ``0``: fully-normalized
            - ``1``: mass coefficients (cm w.e., g/cm\ :sup:`2`)
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

    Returns
    -------
    clm: np.ndarray
        cosine coefficients of the averaging kernel
    slm: np.ndarray
        sine coefficients of the averaging kernel
    """
    # upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # Earth Parameters
    factors = gravity_toolkit.units(lmax=LMAX)
    # extract arrays of kl, hl, and ll Love Numbers
    if (UNITS == 0):
        # Input coefficients are fully-normalized
        dfactor = factors.harmonic(*LOVE).cmwe
    elif (UNITS == 1):
        # Inputs coefficients are mass (cmwe)
        dfactor = np.ones((LMAX+1))
    # average radius of the earth (km)
    rad_e = factors.rad_e/1e5

    # allocate for gaussian function
    gl = np.zeros((LMAX+1))
    # calculate gaussian weights using recursion
    b = np.log(2.0)/(1.0-np.cos(hw/rad_e))
    # weight for degree 0
    gl[0] = (1.0-np.exp(-2.0*b))/b
    # weight for degree 1
    gl[1] = (1.0+np.exp(-2.0*b))/b - (1.0-np.exp(-2.0*b))/b**2
    # valid flag
    valid = True
    # spherical harmonic degree
    l = 2
    # generate Legendre coefficients of Gaussian correlation function
    while (valid and (l <= LMAX)):
        gl[l] = (1.0 - 2.0*l)/b*gl[l-1] + gl[l-2]
        # check validity
        if (gl[l] < CUTOFF):
            gl[l:LMAX+1] = CUTOFF
            valid = False
        # add to counter for spherical harmonic degree
        l += 1

    # Convert sigma to correlation function amplitude
    area = np.copy(gclm[0,0])
    temp_0 = np.zeros((LMAX+1))
    for l in range(0,LMAX+1):# equivalent to 0:LMAX
        mm = np.min([MMAX,l])# find min of MMAX and l
        m = np.arange(0,mm+1)# create m array 0:l or 0:MMAX
        temp_0[l] = (gl[l]/2.0)*np.sum(gclm[l,m]**2 + gslm[l,m]**2)

    # divide by the square of the area under the kernel
    temp = np.sum(temp_0)/area**2
    # signal variance
    sigma_0 = sigma/np.sqrt(temp)

    # Compute averaging kernel coefficients
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1, MMAX+1))
    Ylms.slm = np.zeros((LMAX+1, MMAX+1))
    # for each spherical harmonic degree
    for l in range(0,LMAX+1):# equivalent to 0:lmax
        # inverse of smoothed signal variance in output units
        ldivg = (dfactor[l]**2)/(gl[l]*sigma_0**2)
        # for each valid spherical harmonic order
        mm = np.min([MMAX,l])
        for m in range(0,mm+1):
            temp = 1.0 + 2.0*ldivg*eclm[l,m]**2
            Ylms.clm[l,m] = gclm[l,m]/temp
            temp = 1.0 + 2.0*ldivg*eslm[l,m]**2
            Ylms.slm[l,m] = gslm[l,m]/temp

    # return kernels divided by the area under the kernel
    return Ylms.scale(1.0/area)
