#!/usr/bin/env python
"""
sea_level_equation.py (07/2026)
Solves the sea level equation with the option of including polar motion feedback
Uses a Clenshaw summation to calculate the spherical harmonic summation

CALLING SEQUENCE:
    sea_level = sea_level_equation(loadClm, loadSlm, glon, glat, landmask,
        LMAX=LMAX, LOVE=(hl,kl,ll), BODY_TIDE_LOVE=0, FLUID_LOVE=0, POLAR=True,
        ITERATIONS=6, FILL_VALUE=0)

INPUTS:
    loadClm: input set of cosine spherical harmonics
    loadSlm: input set of sine spherical harmonics
    glon: longitude of the land-sea mask used in the sea level solver
    glat: latitude of the land-sea mask used in the sea level solver
    land_function: land-sea mask for use with the sea level solver (land=1)

OUTPUTS:
    sea_level: spatial field calculated using sea level solver

OPTIONS:
    LMAX: Maximum spherical harmonic degree
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)
    BODY_TIDE_LOVE: Treatment of the body tide Love number
        0: Wahr (1981) and Wahr (1985) values from PREM
        1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
        list or tuple: custom values (k2b,h2b)
    FLUID_LOVE: Treatment of the fluid Love number
        0: Han and Wahr (1989) fluid love number
        1: Munk and MacDonald (1960) secular love number
        2: Munk and MacDonald (1960) fluid love number
        3: Lambeck (1980) fluid love number
        list or tuple: custom value (klf)
    DENSITY: Density of water [g/cm^3]
    POLAR: Include polar feedback
    ITERATIONS: maximum number of iterations for the solver
    PLM: input Legendre polynomials
    FILL_VALUE: value used over land points
    SCALE: scaling factor to prevent underflow in Clenshaw summation

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gen_harmonics.py: Computes spherical harmonic coefficients from a grid
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCES:
    W. E. Farrell, "Deformation of the Earth by surface loads",
        Reviews of Geophysics, 10(3), 761--797, (1972).
        https://doi.org/10.1029/RG010i003p00761
    W. E. Farrell and J. A. Clark, "On Postglacial Sea Level", Geophysical
        Journal of the Royal Astronomical Society, 46(3), 647--667, (1976).
        https://doi.org/10.1111/j.1365-246X.1976.tb01252.x
    D. Han and J. Wahr, "Post-Glacial Rebound Analysis for a Rotating Earth",
        Slow Deformation and Transmission of Stress in the Earth, 49, (1989).
        https://doi.org/10.1029/GM049p0001
    S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw
        summation and the recursive computation of very high degree and
        order normalised associated Legendre functions",
        Journal of Geodesy, 76, 279--299, (2002).
        https://doi.org/10.1007/s00190-002-0216-2
    R. A. Kendall, J. X. Mitrovica, and G. A. Milne,
        "On post-glacial sea level -- II. Numerical formulation and
        comparative results on spherically symmetric models",
        Geophysical Journal International, 161(3), 679--706, (2005).
        https://doi.org/10.1111/j.1365-246X.2005.02553.x
    K. Lambeck, The Earth's Variable Rotation:
        Geophysical Causes and Consequences, First Edition, (1980).
    J. X. Mitrovica and G. A. Milne,
        "On post-glacial sea level: I. General theory",
        Geophysical Journal International, 154(2), 253--267, (2003).
        https://doi.org/10.1046/j.1365-246X.2003.01942.x
    W. H. Munk and G. J. F. MacDonald, The Rotation of the Earth:
        A Geophysical Discussion, First Edition, (1960).
    C. C. Tscherning and K. Poder, "Some Geodetic Applications of Clenshaw
        Summation", Bollettino di Geodesia e Scienze, 4, 349--375, (1982).
    J. M. Wahr, "Body tides on an elliptical, rotating, elastic and
        oceanless Earth", Geophysical Journal of the Royal Astronomical
        Society, 64(3), 677--703, (1981).
        https://doi.org/10.1111/j.1365-246X.1981.tb02690.x
    J. M. Wahr, "Deformation induced by polar motion", Journal of
        Geophysical Research: Solid Earth, 90(B11), 9363--9368, (1985).
        https://doi.org/10.1029/JB090iB11p09363

UPDATE HISTORY:
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
    Updated 06/2025: added option to set the density of sea water (g/cm^3)
    Updated 03/2023: improve typing for variables in docstrings
    Updated 01/2023: refactored associated legendre polynomials
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: using python logging for handling verbose output
        can set custom values for BODY_TIDE_LOVE and FLUID_LOVE
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 08/2020: parameterize float precision to improve computational time
    Updated 07/2020: added function docstrings
    Updated 06/2020: added verbose output of processing run
    Updated 04/2020: reading load love numbers outside of this function
        added option to precompute Legendre Polynomials
    Updated 09/2019: parameters for Munk and MacDonald (1960) fluid love number
    Forked 10/2018 from run_sea_level_equation.py
    Updated 06/2018: using python3 compatible octal and input
    Updated 04/2018: added option to set the gravitational load Love number
        added extrapolation of load love numbers for LMAX > 696
    Updated 03/2018: updated header text.  --format option sets data format
    Updated 02/2018: spherical harmonic order same as spherical harmonic degree
    Updated 09/2017: ocean mask with GEUS Greenland coastlines.  added comments
        more calculations shifted from spatial domain to spherical harmonics
    Updated 08/2017: different landsea mask (small fixes for West Antarctica)
        Clenshaw summation of spherical harmonics to run at 0.5x0.5 degrees
        iterate until convergence of spherical harmonic modulus residuals
    Updated 07/2017: outputs of legendre_polynomials.py include derivatives now
    Updated 04/2017: finished writing algorithms. should be in working condition
        set the permissions mode of the output files with --mode
    Written 09/2016
"""

import logging
import numpy as np
from gravity_toolkit.gen_harmonics import gen_harmonics
from gravity_toolkit.associated_legendre import plm_holmes
from gravity_toolkit.units import units


# PURPOSE: Computes Sea Level Fingerprints including polar motion feedback
def sea_level_equation(
    loadClm,
    loadSlm,
    glon,
    glat,
    land_function,
    LMAX=0,
    LOVE=None,
    BODY_TIDE_LOVE=0,
    FLUID_LOVE=0,
    DENSITY=1.0,
    POLAR=True,
    ITERATIONS=6,
    PLM=None,
    FILL_VALUE=0,
    SCALE=1e-280,
    **kwargs,
):
    r"""
    Solves the sea level equation with the option of including
    polar motion feedback :cite:p:`Farrell:1976hm,Kendall:2005ds,Mitrovica:2003cq`

    Uses a Clenshaw summation to calculate the spherical harmonic
    summation :cite:p:`Holmes:2002ff,Tscherning:1982tu`

    Parameters
    ----------
    loadClm: np.ndarray
        Cosine spherical harmonic coefficients
    loadSlm: np.ndarray
        Sine spherical harmonic coefficients
    glon: np.ndarray
        longitude of the land-sea mask
    glat: np.ndarray
        latitude of the land-sea mask
    land_function: np.ndarray
        land-sea mask with land=1
    LMAX: int, default 0
        Maximum spherical harmonic degree
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    BODY_TIDE_LOVE: int, default 0
        Treatment of the body tide Love number

            - ``0``: :cite:p:`Wahr:1981ea` and :cite:p:`Wahr:1985gr` values from PREM
            - ``1``: :cite:p:`Farrell:1972cm` values from Gutenberg-Bullen oceanic mantle model
            - list or tuple: custom values ``(k2b,h2b)``
    FLUID_LOVE: int, default 0
        Treatment of the fluid Love number

            - ``0``: :cite:p:`Han:1989kj` fluid love number
            - ``1``: :cite:p:`Munk:1960uk` secular love number
            - ``2``: :cite:p:`Munk:1960uk` fluid love number
            - ``3``: :cite:p:`Lambeck:1980um`  fluid love number
            - list or tuple: custom value ``(klf)``
    DENSITY: float, default 1.0
        Density of water [g/cm\ :sup:`3`]
    POLAR: bool, default True
        Include polar feedback
    ITERATIONS: int, default 6
        Maximum number of iterations for the solver
    PLM: np.ndarray or NoneType, default None
        Legendre polynomials
    FILL_VALUE: float, default 0
        Invalid value used over land points
    SCALE: float, default 1e-280
        Scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    sea_level: np.ndarray
        spatial field calculated using sea level solver
    """

    # dimensions of land function
    nphi, nth = np.shape(land_function)
    # calculate colatitude and longitude in radians
    th = np.radians(90.0 - glat)
    phi = np.radians(np.squeeze(glon))
    # calculate ocean function from land function
    ocean_function = 1.0 - land_function
    # indices of the ocean function
    ii, jj = np.nonzero(ocean_function)

    # extract arrays of kl, hl, and ll Love Numbers
    hl, kl, ll = LOVE
    # density of water [g/cm^3]
    rho_water = np.float64(DENSITY)
    # Earth Parameters
    factors = units(lmax=LMAX)
    # Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    # Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    # different treatments of the body tide Love numbers of degree 2
    if isinstance(BODY_TIDE_LOVE, (list, tuple)):
        # use custom defined values
        k2b, h2b = BODY_TIDE_LOVE
    elif BODY_TIDE_LOVE == 0:
        # Wahr (1981) and Wahr (1985) values from PREM
        k2b = 0.298
        h2b = 0.604
    elif BODY_TIDE_LOVE == 1:
        # Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
        k2b = 0.3055
        h2b = 0.6149

    # different treatments of the fluid Love number of gravitational potential
    if isinstance(FLUID_LOVE, (list, tuple)):
        # use custom defined value
        (klf,) = FLUID_LOVE
    elif FLUID_LOVE == 0:
        # Han and Wahr (1989) fluid love number
        # klf = 3.0*G*(C-A)/(rad_e**5*omega**2)
        # klf = 3.0*G*H0*A/(rad_e**5*omega**2)
        G = 6.6740e-11  # gravitational constant [m^3/(kg*s^2)]
        Re = 6.371e6  # mean radius of the Earth [m]
        A_moi = 8.0077e37  # mean equatorial moment of inertia [kg m^2]
        omega = 7.292115e-5  # mean rotation rate of the Earth [radians/s]
        H0 = 0.00328475  # dynamical ellipticity (C_moi-A_moi)/A_moi
        klf = 3.0 * G * H0 * A_moi * (Re**-5) * (omega**-2)
        klf = 0.00328475 / 0.00348118
    if FLUID_LOVE == 1:
        # Munk and MacDonald (1960) secular love number with IERS and PREM values
        GM = 3.98004418e14  # geocentric gravitational constant [m^3/s^2]
        Re = 6.371e6  # mean radius of the Earth [m]
        omega = 7.292115e-5  # mean rotation rate of the Earth [radians/s]
        C_moi = 0.33068  # reduced polar moment of inertia (C/Ma^2)
        H = 1.0 / 305.51  # precessional constant (C_moi-A_moi)/C_moi
        klf = 3.0 * GM * H * C_moi / (Re**3 * omega**2)
    elif FLUID_LOVE == 2:
        # Munk and MacDonald (1960) fluid love number with IERS and WGS84 values
        flat = 1.0 / 298.257223563  # flattening of the WGS84 ellipsoid
        Re = 6.371e6  # mean radius of the Earth [m]
        omega = 7.292115e-5  # mean rotation rate of the Earth [radians/s]
        ge = 9.80665  # standard gravity (mean gravitational acceleration) [m/s^2]
        klf = 2.0 * flat * ge / (omega**2 * Re) - 1.0
    elif FLUID_LOVE == 3:
        # Fluid love number from Lambeck (1980)
        # klf = 3.0*(C-A)*G/(omega**2*rad_e**5) = 3.0*GM*C20/(omega**2*rad_e**3)
        G = 6.672e-11  # gravitational constant [m^3/(kg*s^2)]
        M = 5.974e24  # mass of the Earth [kg]
        R = 6.378140e6  # equatorial radius of the Earth [m]
        Re = 6.3710121e6  # mean radius of the Earth [m]
        omega = 7.292115e-5  # mean rotation rate of the Earth [radians/s]
        A_moi = 0.3295 * M * R**2  # mean equatorial moment of inertia [kg m^2]
        H = 0.003275  # precessional constant (C_moi-A_moi)/C_moi
        C_moi = -A_moi / (H - 1.0)  # mean polar moment of inertia [kg m^2]
        klf = 3.0 * (C_moi - A_moi) * G * (omega**-2) * (Re**-5)
        klf = 0.942

    # calculate coefh and coefp for each degree and order
    # see equation 11 from Tamisiea et al (2010)
    coefh = np.zeros((LMAX + 1, LMAX + 1))
    coefp = np.zeros((LMAX + 1, LMAX + 1))
    for l in range(LMAX + 1):
        m = np.arange(0, l + 1)
        # tilt factor for degree l
        gamma_l = 1.0 + kl[l] - hl[l]
        # coefh and coefp will be the same for all orders except for degree 2
        # and order 1 (if POLAR motion feedback is included)
        coefh[l, m] = 3.0 * rho_water * gamma_l / rho_e / np.float64(2 * l + 1)
        coefp[l, m] = gamma_l / (kl[l] + 1.0)
        # if degree 2 and POLAR parameter is set
        if (l == 2) and POLAR:
            # tilt factor for body tides
            gamma_2b = 1.0 + k2b - h2b
            # calculate coefficient for polar motion feedback and add to coefs
            # For small perturbations in rotation vector: driving potential
            # will be dominated by degree two and order one polar wander
            # effects (quadrantal geometry effects) (Kendall et al., 2005)
            coefpmf = gamma_2b * (1.0 + kl[l]) / (klf - k2b)
            # add effects of polar motion feedback to order 1 coefficients
            coefh[l, 1] += (
                3.0 * rho_water * coefpmf / rho_e / np.float64(2 * l + 1)
            )
            coefp[l, 1] += coefpmf / (kl[l] + 1.0)

    # added option to precompute plms to improve computational speed
    if PLM is None:
        # calculate Legendre polynomials using Holmes and Featherstone relation
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))
    # calculate sin of colatitudes
    gth, gphi = np.meshgrid(th, phi)
    u = np.sin(gth[ii, jj])
    # indices of spherical harmonics for calculating eps
    l1, m1 = np.tril_indices(LMAX + 1)

    # total mass of the surface mass load [g] from harmonics
    tmass = 4.0 * np.pi * (rad_e**3.0) * rho_e * loadClm[0, 0] / 3.0
    # convert ocean function into a series of spherical harmonics
    ocean_Ylms = gen_harmonics(ocean_function, glon, glat, LMAX=LMAX, PLM=PLM)
    # total area of ocean calculated by integrating the ocean function
    ocean_area = 4.0 * np.pi * ocean_Ylms.clm[0, 0]

    # uniform distribution as initial guess of the ocean change following
    # Mitrovica and Peltier (1991) doi:10.1029/91JB01284
    # sea level height change
    sea_height = -tmass / rho_water / rad_e**2 / ocean_area

    # if verbose output: print ocean area and uniform sea level height
    logging.info(f'Total Ocean Area: {ocean_area:0.10g}')
    logging.info(f'Uniform Ocean Height: {sea_height:0.10g}')

    # allocate for output sea level field
    sea_level = np.empty((nphi, nth))
    # complex load spherical harmonics
    loadYlms = loadClm - 1j * loadSlm
    # distribute sea height over ocean harmonics
    height_Ylms = ocean_Ylms * sea_height
    # calculating cos(m*phi) and sin(m*phi) using Euler's formula
    mm = np.arange(0, LMAX + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->pm...', mm, phi))

    # iterate solutions until convergence or reaching total iterations
    n_iter = 1
    # use maximum eps values from Mitrovica and Peltier (1991)
    # Milne and Mitrovica (1998) doi:10.1046/j.1365-246X.1998.1331455.x
    eps = np.inf
    eps_max = 1e-4
    while (eps > eps_max) and (n_iter <= ITERATIONS):
        # zero out the sea level field for this iteration
        sea_level[:, :] = 0.0
        # calculate combined spherical harmonics for Clenshaw summation
        Ylm1 = coefh * height_Ylms.ilm + rad_e * coefp * loadYlms

        # initate summation
        s_m = 0.0
        # iterate to calculate complete summation
        for m in range(LMAX, 0, -1):
            # calculate summation for order m
            a_m = np.sqrt((2.0 * m + 3.0) / (2.0 * m + 2.0))
            cs_m = _clenshaw(np.cos(th), m, Ylm1, LMAX, SCALE=SCALE)
            g = np.einsum('h...,p...->ph...', cs_m, m_phi[:, m])
            # update summation and discard imaginary component
            s_m = a_m * u * s_m + g[ii, jj].real
        # add the l=0/m=0 term
        cs_m = _clenshaw(np.cos(th), 0, Ylm1, LMAX, SCALE=SCALE)
        gs_m = np.kron(np.ones((nphi, 1)), cs_m.real)
        # calculate new sea level for iteration
        sea_level[ii, jj] = np.sqrt(3.0) * u * s_m + gs_m[ii, jj]

        # calculate spherical harmonic field for iteration
        Ylms = gen_harmonics(sea_level, glon, glat, LMAX=LMAX, PLM=PLM)
        # total sea level height for iteration
        # integrated total rmass will differ as sea_level is only over ocean
        # whereas the crustal and gravitational effects are global
        rmass = 4.0 * np.pi * Ylms.clm[0, 0]
        # mass anomaly converted to ocean height to ensure mass conservation
        # (this is the gravitational perturbation (Delta Phi)/g)
        sea_height = (-tmass / rho_water / rad_e**2 - rmass) / ocean_area

        # if verbose output: print iteration, mass and anomaly for convergence
        logging.info(f'Iteration: {n_iter:d}')
        logging.info(f'Integrated Ocean Height: {rmass:0.10g}')
        logging.info(f'Difference from Initial Height: {sea_height:0.10g}')

        # geoid component is split into two parts (Kendall 2005)
        # this part is the spatially uniform shift in the geoid that is
        # constrained by invoking conservation of mass of the surface load
        # Equation 48 of Mitrovica and Peltier (1991)
        # add difference to total sea level field to force mass conservation
        sea_level += sea_height * ocean_function[:, :]
        Ylms += ocean_Ylms * sea_height
        # calculate eps to determine if solution is appropriately converged
        mod1 = np.hypot(height_Ylms.clm, height_Ylms.slm)
        mod2 = np.hypot(Ylms.clm, Ylms.slm)
        eps = np.abs(np.sum(mod2[l1, m1] - mod1[l1, m1]) / np.sum(mod1[l1, m1]))
        # save height harmonics for use in the next iteration
        height_Ylms = Ylms.copy()
        # add 1 to n_iter
        n_iter += 1

    # calculate final total mass for sanity check
    omass = 4.0 * np.pi * (rad_e**2.0) * rho_water * height_Ylms.clm[0, 0]
    # if verbose output: sanity check of masses
    logging.info(f'Original Total Ocean Mass: {-tmass / 1e15:0.10g}')
    logging.info(f'Final Iterated Ocean Mass: {omass / 1e15:0.10g}')

    # set final invalid points to fill value if applicable
    if FILL_VALUE != 0:
        ii, jj = np.nonzero(land_function)
        sea_level[ii, jj] = FILL_VALUE

    # return the sea level spatial field
    return sea_level


# PURPOSE: compute Clenshaw summation of the fully normalized associated
# Legendre's function for constant order m
def _clenshaw(t, m, Ylm1, lmax, SCALE=1e-280):
    """
    Compute conditioned arrays for Clenshaw summation from the fully-normalized
    associated Legendre's function for an order m

    Parameters
    ----------
    t: np.ndarray
        elements ranging from -1 to 1, typically cos(th)
    m: int
        spherical harmonic order
    Ylm1: np.ndarray
        complex form of spherical harmonics
    lmax: int
        maximum spherical harmonic degree
    SCALE: float, default 1e-280
        scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    s_m_c: np.ndarray
        conditioned array for clenshaw summation
    """
    # allocate for output matrix
    N = len(t)
    s_m = np.zeros((N), dtype=np.clongdouble)
    # scaling to prevent overflow
    ylm = SCALE * Ylm1.astype(np.clongdouble)
    # convert lmax and m to float
    lm = np.float64(lmax)
    mm = np.float64(m)
    if m == lmax:
        s_m[:] = np.copy(ylm[lmax, lmax])
    elif m == (lmax - 1):
        a_lm = (
            np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
            * t
        )
        s_m[:] = a_lm * ylm[lmax, lmax - 1] + ylm[lmax - 1, lmax - 1]
    elif (m <= (lmax - 2)) and (m >= 1):
        s_mm_minus_2 = np.copy(ylm[lmax, m])
        a_lm = (
            np.sqrt(
                ((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / ((lm - mm) * (lm + mm))
            )
            * t
        )
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, m]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.float64(l)
            a_lm = (
                np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0 - mm) * (ll + 1.0 + mm))
                )
                * t
            )
            b_lm = np.sqrt(
                ((2.0 * ll + 5.0) * (ll + mm + 1.0) * (ll - mm + 1.0))
                / ((ll + 2.0 - mm) * (ll + 2.0 + mm) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, m]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        s_m[:] = np.copy(s_mm_l)
    elif m == 0:
        s_mm_minus_2 = np.copy(ylm[lmax, 0])
        a_lm = np.sqrt(((2.0 * lm - 1.0) * (2.0 * lm + 1.0)) / (lm * lm)) * t
        s_mm_minus_1 = a_lm * s_mm_minus_2 + ylm[lmax - 1, 0]
        for l in range(lmax - 2, m - 1, -1):
            ll = np.float64(l)
            a_lm = (
                np.sqrt(
                    ((2.0 * ll + 1.0) * (2.0 * ll + 3.0))
                    / ((ll + 1.0) * (ll + 1.0))
                )
                * t
            )
            b_lm = np.sqrt(
                ((2.0 * ll + 5.0) * (ll + 1.0) * (ll + 1.0))
                / ((ll + 2.0) * (ll + 2.0) * (2.0 * ll + 1.0))
            )
            s_mm_l = a_lm * s_mm_minus_1 - b_lm * s_mm_minus_2 + ylm[l, 0]
            s_mm_minus_2 = np.copy(s_mm_minus_1)
            s_mm_minus_1 = np.copy(s_mm_l)
        s_m[:] = np.copy(s_mm_l)
    # return rescaled s_m
    return s_m / SCALE
