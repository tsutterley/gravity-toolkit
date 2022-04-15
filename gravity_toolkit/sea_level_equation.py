#!/usr/bin/env python
u"""
sea_level_equation.py (04/2022)
Solves the sea level equation with the option of including polar motion feedback
Uses a Clenshaw summation to calculate the spherical harmonic summation

CALLING SEQUENCE:
    sea_level = sea_level_equation(loadClm, loadSlm, glon, glat, landmask,
        LMAX=LMAX, LOVE=(hl,kl,ll), BODY_TIDE_LOVE=0, FLUID_LOVE=0, POLAR=True,
        INTERATIONS=6, FILL_VALUE=0)

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
    POLAR: Include polar feedback
    ITERATIONS: maximum number of iterations for the solver
    PLM: input Legendre polynomials
    FILL_VALUE: value used over land points
    ASTYPE: floating point precision for calculating Clenshaw summation
    SCALE: scaling factor to prevent underflow in Clenshaw summation

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
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
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.units import units

#-- PURPOSE: Computes Sea Level Fingerprints including polar motion feedback
def sea_level_equation(loadClm, loadSlm, glon, glat, land_function, LMAX=0,
    LOVE=None, BODY_TIDE_LOVE=0, FLUID_LOVE=0, POLAR=True, ITERATIONS=6,
    PLM=None, FILL_VALUE=0, ASTYPE=np.longdouble, SCALE=1e-280, **kwargs):
    """
    Solves the sea level equation with the option of including
    polar motion feedback [Farrell1976]_ [Kendall2005]_ [Mitrovica2003]_

    Uses a Clenshaw summation to calculate the spherical harmonic
    summation [Holmes2002]_ [Tscherning1982]_

    Parameters
    ----------
    loadClm: float
        Cosine spherical harmonic coefficients
    loadSlm: float
        Sine spherical harmonic coefficients
    glon: float
        longitude of the land-sea mask
    glat: float
        latitude of the land-sea mask
    land_function: int
        land-sea mask with land=1
    LMAX: int, default 0
        Maximum spherical harmonic degree
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    BODY_TIDE_LOVE: int, default 0
        Treatment of the body tide Love number

            - ``0``: [Wahr1981]_ and [Wahr1985]_ values from PREM
            - ``1``: [Farrell1972]_ values from Gutenberg-Bullen oceanic mantle model
            - list or tuple: custom values ``(k2b,h2b)``
    FLUID_LOVE: int, default 0
        Treatment of the fluid Love number

            - ``0``: [Han1989]  fluid love number
            - ``1``: [Munk1960]_ secular love number
            - ``2``: [Munk1960]_ fluid love number
            - ``3``: [Lambeck1980]  fluid love number
            - list or tuple: custom value ``(klf)``
    POLAR: bool, default True
        Include polar feedback
    ITERATIONS: int, default 6
        Maximum number of iterations for the solver
    PLM: float or NoneType, default None
        Legendre polynomials
    FILL_VALUE: float, default 0
        Invalid value used over land points
    ASTYPE: obj, default np.longdouble
        Floating point precision for calculating Clenshaw summation
    SCALE: float, default 1e-280
        Scaling factor to prevent underflow in Clenshaw summation

    Returns
    -------
    sea_level: float
        spatial field calculated using sea level solver

    References
    ----------
    .. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads",
        *Reviews of Geophysics*, 10(3), 761--797, (1972).
        `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_
    .. [Farrell1976] W. E. Farrell and J. A. Clark, "On Postglacial Sea Level",
        *Geophysical Journal of the Royal Astronomical Society*, 46(3), 647--667,
        (1976). `doi: 10.1111/j.1365-246X.1976.tb01252.x <https://doi.org/10.1111/j.1365-246X.1976.tb01252.x>`_
    .. [Han1989] D. Han and J. Wahr, "Post-Glacial Rebound Analysis for a
        Rotating Earth", *Slow Deformation and Transmission of Stress in the Earth*,
        49, (1989). `doi: 10.1029/GM049p0001 <https://doi.org/10.1029/GM049p0001>`_
    .. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach
        to the Clenshaw summation and the recursive computation of very high
        degree and order normalised associated Legendre functions",
        *Journal of Geodesy*, 76, 279--299, (2002).
        `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_
    .. [Kendall2005] R. A. Kendall, J. X. Mitrovica, and G. A. Milne,
        "On post-glacial sea level -- II. Numerical formulation and comparative
        results on spherically symmetric models", *Geophysical Journal International*,
        161(3), 679--706, (2005).
        `doi: 10.1111/j.1365-246X.2005.02553.x <https://doi.org/10.1111/j.1365-246X.2005.02553.x>`_
    .. [Lambeck1980] K. Lambeck, *The Earth's Variable Rotation:
        Geophysical Causes and Consequences*, First Edition, (1980).
    .. [Mitrovica2003] J. X. Mitrovica and G. A. Milne,
        "On post-glacial sea level: I. General theory",
        *Geophysical Journal International*, 154(2), 253--267, (2003).
        `doi: 10.1046/j.1365-246X.2003.01942.x <https://doi.org/10.1046/j.1365-246X.2003.01942.x>`_
    .. [Munk1960] W. H. Munk and G. J. F. MacDonald,
        *The Rotation of the Earth: A Geophysical Discussion*, First Edition, (1960).
    .. [Tscherning1982] C. C. Tscherning and K. Poder,
        "Some Geodetic Applications of Clenshaw Summation",
        *Bollettino di Geodesia e Scienze*, 4, 349--375, (1982).
    .. [Wahr1981] J. M. Wahr, "Body tides on an elliptical, rotating,
        elastic and oceanless Earth", *Geophysical Journal of the Royal
        Astronomical Society*, 64(3), 677--703, (1981).
        `doi: 10.1111/j.1365-246X.1981.tb02690.x <https://doi.org/10.1111/j.1365-246X.1981.tb02690.x>`_
    .. [Wahr1985] J. M. Wahr, "Deformation induced by polar motion",
        *Journal of Geophysical Research: Solid Earth*, 90(B11), 9363--9368,
        (1985). `doi: 10.1029/JB090iB11p09363 <https://doi.org/10.1029/JB090iB11p09363>`_
    """

    #-- dimensions of land function
    nphi,nth = np.shape(land_function)
    #-- calculate colatitude and longitude in radians
    th = (90.0 - glat)*np.pi/180.0
    phi = np.squeeze(glon*np.pi/180.0)
    #-- calculate ocean function from land function
    ocean_function = 1.0 - land_function
    #-- indices of the ocean function
    ii,jj = np.nonzero(ocean_function)

    #-- extract arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = LOVE
    #-- density of water [g/cm^3]
    rho_water = 1.0
    #-- Earth Parameters
    factors = units(lmax=LMAX)
    #-- Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    #-- Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    #-- different treatments of the body tide Love numbers of degree 2
    if isinstance(BODY_TIDE_LOVE,(list,tuple)):
        #-- use custom defined values
        k2b,h2b = BODY_TIDE_LOVE
    elif (BODY_TIDE_LOVE == 0):
        #-- Wahr (1981) and Wahr (1985) values from PREM
        k2b = 0.298
        h2b = 0.604
    elif (BODY_TIDE_LOVE == 1):
        #-- Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
        k2b = 0.3055
        h2b = 0.6149

    #-- different treatments of the fluid Love number of gravitational potential
    if isinstance(FLUID_LOVE,(list,tuple)):
        #-- use custom defined value
        klf, = FLUID_LOVE
    elif (FLUID_LOVE == 0):
        #-- Han and Wahr (1989) fluid love number
        #-- klf = 3.0*G*(C-A)/(rad_e**5*omega**2)
        G = 6.6740e-11#-- gravitational constant [m^3/(kg*s^2)]
        Re = 6.371e6#-- mean radius of the Earth [m]
        A_moi = 8.0077e+37#-- mean equatorial moment of inertia [kg m^2]
        omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
        ef = 0.00328475#-- dynamical ellipticity (C_moi-A_moi)/A_moi
        C_moi = A_moi*(1.0 + ef)#-- mean polar moment of inertia [kg m^2]
        klf = 3.0*G*(C_moi-A_moi)*(Re**-5)*(omega**-2)
        klf = 0.00328475/0.00348118
    if (FLUID_LOVE == 1):
        #-- Munk and MacDonald (1960) secular love number with IERS and PREM values
        GM = 3.98004418e14#-- geocentric gravitational constant [m^3/s^2]
        Re = 6.371e6#-- mean radius of the Earth [m]
        omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
        C_moi = 0.33068#-- reduced polar moment of inertia (C/Ma^2)
        H = 1.0/305.51#-- precessional constant (C_moi-A_moi)/C_moi
        klf = 3.0*GM*H*C_moi/(Re**3*omega**2)
    elif (FLUID_LOVE == 2):
        #-- Munk and MacDonald (1960) fluid love number with IERS and WGS84 values
        flat = 1.0/298.257223563#-- flattening of the WGS84 ellipsoid
        Re = 6.371e6#-- mean radius of the Earth [m]
        omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
        ge = 9.80665#-- standard gravity (mean gravitational acceleration) [m/s^2]
        klf = 2.0*flat*ge/(omega**2*Re) - 1.0
    elif (FLUID_LOVE == 3):
        #-- Fluid love number from Lambeck (1980)
        #-- klf = 3.0*(C-A)*G/(omega**2*rad_e**5) = 3.0*GM*C20/(omega**2*rad_e**3)
        G = 6.672e-11#-- gravitational constant [m^3/(kg*s^2)]
        M = 5.974e+24#-- mass of the Earth [kg]
        R = 6.378140e6#-- equatorial radius of the Earth [m]
        Re = 6.3710121e6#-- mean radius of the Earth [m]
        omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
        A_moi = 0.3295*M*R**2#-- mean equatorial moment of inertia [kg m^2]
        H = 0.003275#-- precessional constant (C_moi-A_moi)/C_moi
        C_moi = -A_moi/(H-1.0)#-- mean polar moment of inertia [kg m^2]
        klf = 3.0*(C_moi-A_moi)*G*(omega**-2)*(Re**-5)
        klf = 0.942

    #-- calculate coefh and coefp for each degree and order
    #-- see equation 11 from Tamisiea et al (2010)
    coefh = np.zeros((LMAX+1,LMAX+1))
    coefp = np.zeros((LMAX+1,LMAX+1))
    for l in range(LMAX+1):
        #-- coefh and coefp will be the same for all orders except for degree 2
        #-- and order 1 (if POLAR motion feedback is included)
        m = np.arange(0,l+1)
        coefh[l,m] = 3.0*rho_water*(1.0 + kl[l] - hl[l])/rho_e/np.float64(2*l+1)
        coefp[l,m] = (1.0 + kl[l] - hl[l])/(kl[l] + 1.0)
        #-- if degree 2 and POLAR parameter is set
        if (l == 2) and POLAR:
            #-- calculate coefficient for polar motion feedback and add to coefs
            #-- For small perturbations in rotation vector: driving potential
            #-- will be dominated by degree two and order one polar wander
            #-- effects (quadrantal geometry effects) (Kendall et al., 2005)
            coefpmf = (1.0 + k2b - h2b)*(1.0 + kl[l])/(klf - k2b)
            #-- add effects of polar motion feedback to order 1 coefficients
            coefh[l,1] += 3.0*rho_water*coefpmf/rho_e/np.float64(2*l+1)
            coefp[l,1] += coefpmf/(kl[l] + 1.0)

    #-- added option to precompute plms to improve computational speed
    if PLM is None:
        #-- calculate Legendre polynomials using Holmes and Featherstone relation
        PLM,dPLM = plm_holmes(LMAX,np.cos(th))
    #-- calculate sin of colatitudes
    gth,gphi = np.meshgrid(th, phi)
    u = np.sin(gth[ii,jj])
    #-- indices of spherical harmonics for calculating eps
    l1,m1 = np.tril_indices(LMAX+1)

    #-- total mass of the surface mass load [g] from harmonics
    tmass = 4.0*np.pi*(rad_e**3.0)*rho_e*loadClm[0,0]/3.0
    #-- convert ocean function into a series of spherical harmonics
    ocean_Ylms = gen_harmonics(ocean_function,glon,glat,LMAX=LMAX,PLM=PLM)
    #-- total area of ocean calculated by integrating the ocean function
    ocean_area = 4.0*np.pi*ocean_Ylms.clm[0,0]

    #-- uniform distribution as initial guess of the ocean change following
    #-- Mitrovica and Peltier (1991) doi:10.1029/91JB01284
    #-- sea level height change
    sea_height = -tmass/rho_water/rad_e**2/ocean_area

    #-- if verbose output: print ocean area and uniform sea level height
    logging.info('Total Ocean Area: {0:0.10g}'.format(ocean_area))
    logging.info('Uniform Ocean Height: {0:0.10g}'.format(sea_height))

    #-- distribute sea height over ocean harmonics
    height_Ylms = ocean_Ylms.scale(sea_height)
    #-- iterate solutions until convergence or reaching total iterations
    n_iter = 1
    #-- use maximum eps values from Mitrovica and Peltier (1991)
    #-- Milne and Mitrovica (1998) doi:10.1046/j.1365-246X.1998.1331455.x
    eps = np.inf
    eps_max = 1e-4
    while (eps > eps_max) and (n_iter <= ITERATIONS):
        #-- allocate for sea level field of iteration
        sea_level = np.zeros((nphi,nth))
        #-- calculate combined spherical harmonics for Clenshaw summation
        clm1 = coefh*height_Ylms.clm + rad_e*coefp*loadClm
        slm1 = coefh*height_Ylms.slm + rad_e*coefp*loadSlm
        #-- calculate clenshaw summations over colatitudes
        s_m_c = np.zeros((nth,LMAX*2+2))
        for m in range(LMAX, -1, -1):
            s_m_c[:,2*m:2*m+2] = clenshaw_s_m(np.cos(th), m, clm1, slm1, LMAX,
                ASTYPE=ASTYPE, SCALE=SCALE)

        #-- calculate cos(phi)
        cos_phi_2 = 2.0*np.cos(phi)
        #-- matrix of cos/sin m*phi summation
        cos_m_phi = np.zeros((nphi,LMAX+2),dtype=ASTYPE)
        sin_m_phi = np.zeros((nphi,LMAX+2),dtype=ASTYPE)
        #-- initialize matrix with values at lmax+1 and lmax
        cos_m_phi[:,LMAX+1] = np.cos(ASTYPE(LMAX + 1)*phi)
        sin_m_phi[:,LMAX+1] = np.sin(ASTYPE(LMAX + 1)*phi)
        cos_m_phi[:,LMAX] = np.cos(ASTYPE(LMAX)*phi)
        sin_m_phi[:,LMAX] = np.sin(ASTYPE(LMAX)*phi)
        #-- calculate summation
        gc=np.multiply(s_m_c[np.newaxis,:,2*LMAX],cos_m_phi[:,np.newaxis,LMAX])
        gs=np.multiply(s_m_c[np.newaxis,:,2*LMAX+1],sin_m_phi[:,np.newaxis,LMAX])
        s_m = gc[ii,jj] + gs[ii,jj]
        #-- iterate to calculate complete summation
        for m in range(LMAX-1, 0, -1):
            cos_m_phi[:,m] = cos_phi_2*cos_m_phi[:,m+1] - cos_m_phi[:,m+2]
            sin_m_phi[:,m] = cos_phi_2*sin_m_phi[:,m+1] - sin_m_phi[:,m+2]
            a_m = np.sqrt((2.0*m+3.0)/(2.0*m+2.0))
            gc=np.multiply(s_m_c[np.newaxis,:,2*m],cos_m_phi[:,np.newaxis,m])
            gs=np.multiply(s_m_c[np.newaxis,:,2*m+1],sin_m_phi[:,np.newaxis,m])
            s_m = a_m*u*s_m + gc[ii,jj] + gs[ii,jj]
        #-- calculate new sea level for iteration
        gsmc,gcmp = np.meshgrid(s_m_c[:,0],cos_m_phi[:,0])
        sea_level[ii,jj] = np.sqrt(3.0)*u*s_m + gsmc[ii,jj]

        #-- calculate spherical harmonic field for iteration
        Ylms = gen_harmonics(sea_level, glon, glat, LMAX=LMAX, PLM=PLM)
        #-- total sea level height for iteration
        #-- integrated total rmass will differ as sea_level is only over ocean
        #-- whereas the crustal and gravitational effects are global
        rmass = 4.0*np.pi*Ylms.clm[0,0]
        #-- mass anomaly converted to ocean height to ensure mass conservation
        #-- (this is the gravitational perturbation (Delta Phi)/g)
        sea_height = (-tmass/rho_water/rad_e**2 - rmass)/ocean_area

        #-- if verbose output: print iteration, mass and anomaly for convergence
        logging.info('Iteration: {0:d}'.format(n_iter))
        logging.info('Integrated Ocean Height: {0:0.10g}'.format(rmass))
        logging.info('Difference from Initial Height: {0:0.10g}'.format(sea_height))

        #-- geoid component is split into two parts (Kendall 2005)
        #-- this part is the spatially uniform shift in the geoid that is
        #-- constrained by invoking conservation of mass of the surface load
        #-- Equation 48 of Mitrovica and Peltier (1991)
        #-- add difference to total sea level field to force mass conservation
        sea_level += sea_height*ocean_function[:,:]
        uniform_Ylms = ocean_Ylms.scale(sea_height)
        Ylms.add(uniform_Ylms)
        #-- calculate eps to determine if solution is appropriately converged
        mod1 = np.sqrt(height_Ylms.clm**2 + height_Ylms.slm**2)
        mod2 = np.sqrt(Ylms.clm**2 + Ylms.slm**2)
        eps = np.abs(np.sum(mod2[l1,m1] - mod1[l1,m1])/np.sum(mod1[l1,m1]))
        #-- save height harmonics for use in the next iteration
        height_Ylms = Ylms.copy()
        #-- add 1 to n_iter
        n_iter += 1

    #-- calculate final total mass for sanity check
    omass = 4.0*np.pi*(rad_e**2.0)*rho_water*height_Ylms.clm[0,0]
    #-- if verbose output: sanity check of masses
    logging.info('Original Total Ocean Mass: {0:0.10g}'.format(-tmass/1e15))
    logging.info('Final Iterated Ocean Mass: {0:0.10g}'.format(omass/1e15))

    #-- set final invalid points to fill value if applicable
    if (FILL_VALUE != 0):
        ii,jj = np.nonzero(land_function)
        sea_level[ii,jj] = FILL_VALUE

    #-- return the sea level spatial field
    return sea_level

#-- PURPOSE: compute Clenshaw summation of the fully normalized associated
#-- Legendre's function for constant order m
def clenshaw_s_m(t, m, clm1, slm1, lmax, ASTYPE=np.longdouble, SCALE=1e-280):
    #-- allocate for output matrix
    N = len(t)
    s_m = np.zeros((N,2),dtype=ASTYPE)
    #-- scaling to prevent overflow
    clm = SCALE*clm1.astype(ASTYPE)
    slm = SCALE*slm1.astype(ASTYPE)
    #-- convert lmax and m to float
    lm = ASTYPE(lmax)
    mm = ASTYPE(m)
    if (m == lmax):
        s_m[:,0] = np.copy(clm[lmax,lmax])
        s_m[:,1] = np.copy(slm[lmax,lmax])
    elif (m == (lmax-1)):
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t
        s_m[:,0] = a_lm*clm[lmax,lmax-1] + clm[lmax-1,lmax-1]
        s_m[:,1] = a_lm*slm[lmax,lmax-1] + slm[lmax-1,lmax-1]
    elif ((m <= (lmax-2)) and (m >= 1)):
        s_mm_c_pre_2 = np.copy(clm[lmax,m])
        s_mm_s_pre_2 = np.copy(slm[lmax,m])
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/((lm-mm)*(lm+mm)))*t
        s_mm_c_pre_1 = a_lm*s_mm_c_pre_2 + clm[lmax-1,m]
        s_mm_s_pre_1 = a_lm*s_mm_s_pre_2 + slm[lmax-1,m]
        for l in range(lmax-2, m-1, -1):
            ll = ASTYPE(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0-mm)*(ll+1.0+mm)))*t
            b_lm=np.sqrt(((2.*ll+5.)*(ll+mm+1.)*(ll-mm+1.))/((ll+2.-mm)*(ll+2.+mm)*(2.*ll+1.)))
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,m]
            s_mm_s = a_lm * s_mm_s_pre_1 - b_lm * s_mm_s_pre_2 + slm[l,m]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_s_pre_2 = np.copy(s_mm_s_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
            s_mm_s_pre_1 = np.copy(s_mm_s)
        s_m[:,0] = np.copy(s_mm_c)
        s_m[:,1] = np.copy(s_mm_s)
    elif (m == 0):
        s_mm_c_pre_2 = np.copy(clm[lmax,0])
        a_lm = np.sqrt(((2.0*lm-1.0)*(2.0*lm+1.0))/(lm*lm))*t
        s_mm_c_pre_1 = a_lm * s_mm_c_pre_2 + clm[lmax-1,0]
        for l in range(lmax-2, m-1, -1):
            ll = ASTYPE(l)
            a_lm=np.sqrt(((2.0*ll+1.0)*(2.0*ll+3.0))/((ll+1.0)*(ll+1.0)))*t
            b_lm=np.sqrt(((2.0*ll+5.0)*(ll+1.0)*(ll+1.0))/((ll+2.0)*(ll+2.0)*(2.0*ll+1.0)))
            s_mm_c = a_lm * s_mm_c_pre_1 - b_lm * s_mm_c_pre_2 + clm[l,0]
            s_mm_c_pre_2 = np.copy(s_mm_c_pre_1)
            s_mm_c_pre_1 = np.copy(s_mm_c)
        s_m[:,0] = np.copy(s_mm_c)
    #-- return rescaled s_m
    return s_m/SCALE
