#!/usr/bin/env python
u"""
gen_averaging_kernel.py
Original IDL code gen_wclms_me.pro written by Sean Swenson
Adapted by Tyler Sutterley (04/2022)

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
    units.py: class for converting spherical harmonic data to specific units

REFERENCES:
    Swenson and Wahr, "Methods for inferring regional surface-mass anomalies
        from Gravity Recovery and Climate Experiment (GRACE) measurements of
        time-variable gravity," Journal of Geophysical Research: Solid Earth,
        107(B9), (2002). https://doi.org/10.1029/2001JB000576

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2021: using units module for Earth parameters
    Updated 04/2020: reading load love numbers outside of this function
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 05/2013
"""
import numpy as np
import gravity_toolkit.units

def gen_averaging_kernel(gclm, gslm, eclm, eslm, sigma, hw,
    LMAX=60, MMAX=None, UNITS=0, LOVE=None):
    """
    Generates averaging kernel coefficients which minimize the
    total error following [Swenson2002]_

    Uses a normalized form of the Gaussian averaging function
    from [Jekeli1981]_

    Parameters
    ----------
    gclm: float
        cosine spherical harmonics of exact averaging kernel
    gslm: float
        sine spherical harmonics of exact averaging kernel
    eclm: float
        measurement error in the cosine harmonics
    eslm: float
        measurement error in the sine harmonics
    sigma: float
        variance of the surface mass signal
    hw: float
        Gaussian radius of the kernel in kilometers
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    UNITS: int, default 0
        Input data units

            - ``0``: fully-normalized
            - ``1``: mass coefficients (cm w.e., g/cm\ :sup:`2`)
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

    Returns
    -------
    clm: float
        cosine coefficients of the averaging kernel
    slm: float
        sine coefficients of the averaging kernel

    References
    ----------
    .. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth
        the Earth's Gravity Field", NASA Grant No. NGR 36-008-161,
        OSURF Proj. No. 783210, 48 pp., (1981).

    .. [Swenson2002] S. Swenson and J. Wahr, "Methods for inferring regional
        surface‐mass anomalies from Gravity Recovery and Climate Experiment
        (GRACE) measurements of time‐variable gravity", *Journal of
        Geophysical Research: Solid Earth*, 107(B9), 2193, (2002).
        `doi: 10.1029/2001JB000576 <https://doi.org/10.1029/2001JB000576>`_
    """
    #-- upper bound of spherical harmonic orders (default = LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    #-- Earth Parameters
    #-- extract arrays of kl, hl, and ll Love Numbers
    dfactor = gravity_toolkit.units(lmax=LMAX).harmonic(*LOVE)
    #-- average radius of the earth (km)
    rad_e = dfactor.rad_e/1e5

    #-- allocate for gaussian function
    gl = np.zeros((LMAX+1))
    #-- calculate gaussian weights using recursion
    b = np.log(2.0)/(1.0-np.cos(hw/rad_e))
    #-- weight for degree 0
    gl[0] = (1.0-np.exp(-2.0*b))/b
    #-- weight for degree 1
    gl[1] = (1.0+np.exp(-2.0*b))/b - (1.0-np.exp(-2.0*b))/b**2
    #-- valid flag
    valid = True
    #-- spherical harmonic degree
    l = 2
    #-- generate Legendre coefficients of Gaussian correlation function
    while (valid and (l <= LMAX)):
        gl[l] = (1.0 - 2.0*l)/b*gl[l-1] + gl[l-2]
        #-- check validity
        if (np.abs(gl[l]) < 1.0e-15):
            gl[l:LMAX+1] = 1.0e-15
            valid = False
        #-- add to counter for spherical harmonic degree
        l += 1

    #-- Convert sigma to correlation function amplitude
    area = np.copy(gclm[0,0])
    temp_0 = np.zeros((LMAX+1))
    for l in range(0,LMAX+1):#-- equivalent to 0:LMAX
        mm = np.min([MMAX,l])#-- find min of MMAX and l
        m = np.arange(0,mm+1)#-- create m array 0:l or 0:MMAX
        temp_0[l] = (gl[l]/2.0)*np.sum(gclm[l,m]**2 + gslm[l,m]**2)

    #-- divide by the square of the area under the kernel
    temp = np.sum(temp_0)/area**2
    #-- signal variance
    sigma_0 = sigma/np.sqrt(temp)

    #-- Compute averaging kernel coefficients
    wclm = np.zeros((LMAX+1,MMAX+1))
    wslm = np.zeros((LMAX+1,MMAX+1))
    #-- for each spherical harmonic degree
    for l in range(0,LMAX+1):#-- equivalent to 0:lmax
        if (UNITS == 0):
            #-- Input coefficients are fully-normalized
            cmwe = dfactor.cmwe[l]
            ldivg = (cmwe**2)/(gl[l]*sigma_0**2)
        elif (UNITS == 1):
            #-- Inputs coefficients are mass (cmwe)
            ldivg = 1.0/(gl[l]*sigma_0**2)
        #-- for each valid spherical harmonic order
        mm = np.min([MMAX,l])
        for m in range(0,mm+1):
            temp = 1.0 + 2.0*ldivg*eclm[l,m]**2
            wclm[l,m] = gclm[l,m]/temp
            temp = 1.0 + 2.0*ldivg*eslm[l,m]**2
            wslm[l,m] = gslm[l,m]/temp

    #-- return kernels divided by the area under the kernel
    return {'clm':wclm/area, 'slm':wslm/area}
