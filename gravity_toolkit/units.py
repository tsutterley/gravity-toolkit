#!/usr/bin/env python
u"""
units.py
Written by Tyler Sutterley (05/2024)
Contributions by Hugo Lecomte

Class for converting GRACE/GRACE-FO Level-2 data to specific units

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 09/2023: added property for the approximate mass of the Earth
    Updated 03/2023: include option to not compensate for elastic deformation
        include option to include effects for Earth's oblateness
        added functions for getting unit attributes for known types
        added output for Earth's average angular velocity (omega)
        improve typing for variables in docstrings
    Updated 02/2023: fixed case where maximum spherical harmonic degree is 0
    Updated 01/2023: added function to retrieve named units
    Updated 12/2022: set average Earth's density and radius as class properties
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: made semi-major axis and ellipsoidal flattening arguments
    Updated 07/2020: added class docstring
    Updated 04/2020: include earth parameters as attributes
    Written 03/2020
"""
from __future__ import annotations

import numpy as np

class units(object):
    r"""
    Class for converting spherical harmonic and spatial data to specific units

    Parameters
    ----------
    lmax: int or NoneType, default None
        Upper bound of Spherical Harmonic Degrees
    a_axis: float, default 6.378137e8
        semi-major axis of the Earth's ellipsoid in cm
    flat: float, default 1.0/298.257223563
        flattening of the Earth's ellipsoid

    Attributes
    ----------
    GM: float
        Gravitational Constant of the Earth in cm\ :sup:`3`/s\ :sup:`2`
    g_wmo: float
        standard gravitational acceleration in cm/s\ :sup:`2`
    omega: float
        Angular velocity of the Earth in radians/s
    l: int
        spherical harmonic degree up to ``lmax``
    """
    np.seterr(invalid='ignore')
    def __init__(self,
                 lmax: int | None = None,
                 a_axis: float = 6.378137e8,
                 flat: float = 1.0/298.257223563
        ):
        # Earth Parameters
        # universal gravitational constant [dyn*cm^2/g^2]
        self.G = 6.67430e-8
        # WGS84 Gravitational Constant of the Earth [cm^3/s^2]
        GM_e = 3986004.418e14
        # Gravitational Constant of the Earth's atmosphere
        GM_atm = 3.5e14
        # Gravitational Constant of the Earth (w/o atm)
        self.GM = GM_e - GM_atm
        # Semi-major axis of the ellipsoid [cm]
        self.a_axis = a_axis
        # Flattening of the ellipsoid
        self.flat = flat
        # standard gravitational acceleration (World Meteorological Organization)
        self.g_wmo = 980.665
        # average angular velocity of the Earth [rad/s]
        self.omega = 7292115e-11
        # Degree dependent parameters
        self.norm = None
        self.cmwe = None
        self.mmwe = None
        self.mmGH = None
        self.mmCU = None
        self.mmCH = None
        self.cmVCU = None
        self.mVCU = None
        self.microGal = None
        self.mbar = None
        self.Pa = None
        self.lmax = lmax
        # calculate spherical harmonic degree (0 is falsy)
        self.l = np.arange(self.lmax+1) if (self.lmax is not None) else None

    @property
    def b_axis(self) -> float:
        """semi-minor axis of the Earth's ellipsoid in cm
        """
        return (1.0 - self.flat)*self.a_axis

    @property
    def rad_e(self) -> float:
        """average radius of the Earth in cm with the same volume as the ellipsoid
        """
        return self.a_axis*(1.0 - self.flat)**(1.0/3.0)

    @property
    def rho_e(self) -> float:
        r"""average density of the Earth in g/cm\ :sup:`3`
        """
        return 0.75*self.GM/(self.G*np.pi*self.rad_e**3)

    @property
    def mass(self) -> float:
        """approximate mass of the Earth in g
        """
        return (4.0/3.0)*np.pi*self.rho_e*self.rad_e**3

    def harmonic(self, hl, kl, ll, **kwargs):
        """
        Calculates degree dependent factors for converting spherical
        harmonic units from :cite:p:`Wahr:1998hy`

        Parameters
        ----------
        hl: float
            Vertical Displacement load Love numbers
        kl: float
            Gravitational Potential load Love numbers
        ll: float
            Horizontal Displacement load Love numbers
        include_elastic: bool, default True
            Include compensation for elastic response
        include_ellipsoidal: bool, default False
            Include consideration for Earth's oblateness

        Attributes
        ----------
        norm: float
            fully normalized spherical harmonics
        cmwe: float
            centimeters water equivalent
        mmwe: float
            millimeters water equivalent
        mmGH: float
            millimeters geoid height
        mmCU: float
            millimeters elastic crustal deformation (uplift)
        mmCH: float
            millimeters elastic crustal deformation (horizontal)
        cmVCU: float
            centimeters viscoelastic crustal uplift from :cite:p:`Wahr:2000ek`
        mVCU: float
            meters viscoelastic crustal uplift from :cite:p:`Wahr:2000ek`
        microGal: float
            microGal gravity perturbations
        mbar: float
            millibar equivalent surface pressure
        Pa: float
            pascals equivalent surface pressure
        """
        # set default keyword arguments
        kwargs.setdefault('include_elastic', True)
        kwargs.setdefault('include_ellipsoidal', False)
        fraction = np.ones((self.lmax+1))
        # compensate for elastic deformation within the solid earth
        if kwargs['include_elastic']:
            fraction += kl[self.l]
        # include effects for Earth's oblateness
        if kwargs['include_ellipsoidal']:
            fraction /= (1.0 - self.flat)
        # degree dependent coefficients
        # norm, fully normalized spherical harmonics
        self.norm = np.ones((self.lmax+1))
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe = self.rho_e*self.rad_e*(2.0*self.l+1.0)/fraction/3.0
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe = 10.0*self.rho_e*self.rad_e*(2.0*self.l+1.0)/fraction/3.0
        # mmGH, millimeters geoid height
        self.mmGH = np.ones((self.lmax+1))*(10.0*self.rad_e)
        # mmCU, millimeters elastic crustal deformation (uplift)
        self.mmCU = 10.0*self.rad_e*hl[self.l]/fraction
        # mmCH, millimeters elastic crustal deformation (horizontal)
        self.mmCH = 10.0*self.rad_e*ll[self.l]/fraction
        # cmVCU, centimeters viscoelastic crustal uplift
        self.cmVCU = self.rad_e*(2.0*self.l+1.0)/2.0
        # mVCU, meters viscoelastic crustal uplift
        self.mVCU = self.rad_e*(2.0*self.l+1.0)/200.0
        # microGal, microGal gravity perturbations
        self.microGal = 1.e6*self.GM*(self.l+1.0)/(self.rad_e**2.0)
        # mbar, millibar equivalent surface pressure
        self.mbar = self.g_wmo*self.rho_e*self.rad_e*(2.0*self.l+1.0)/fraction/3e3
        # Pa, pascals equivalent surface pressure
        self.Pa = self.g_wmo*self.rho_e*self.rad_e*(2.0*self.l+1.0)/fraction/30.0
        # return the degree dependent unit conversions
        return self

    def spatial(self, hl, kl, ll, **kwargs):
        """
        Calculates degree dependent factors for converting spatial units
        from :cite:p:`Wahr:1998hy`

        Parameters
        ----------

        hl: float
            Vertical Displacement load Love numbers
        kl: float
            Gravitational Potential load Love numbers
        ll: float
            Horizontal Displacement load Love numbers
        include_elastic: bool, default True
            Include compensation for elastic response

        Attributes
        ----------
        norm: float
            fully normalized spherical harmonics
        cmwe: float
            centimeters water equivalent
        mmwe: float
            millimeters water equivalent
        mmGH: float
            millimeters geoid height
        microGal: float
            microGal gravity perturbations
        """
        # set default keyword arguments
        kwargs.setdefault('include_elastic', True)
        fraction = np.ones((self.lmax+1))
        # compensate for elastic deformation within the solid earth
        if kwargs['include_elastic']:
            fraction += kl[self.l]
        # degree dependent coefficients
        # norm, fully normalized spherical harmonics
        self.norm = np.ones((self.lmax+1))
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe = 3.0*fraction/(1.0+2.0*self.l)/(4.0*np.pi*self.rad_e*self.rho_e)
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe = 3.0*fraction/(1.0+2.0*self.l)/(40.0*np.pi*self.rad_e*self.rho_e)
        # mmGH, millimeters geoid height
        self.mmGH = np.ones((self.lmax+1))/(4.0*np.pi*self.rad_e)
        # microGal, microGal gravity perturbations
        self.microGal = (self.rad_e**2.0)/(4.0*np.pi*1.e6*self.GM)/(self.l+1.0)
        # return the degree dependent unit conversions
        return self

    def get(self, var: str):
        """
        Get the degree dependent factors for a specific unit

        Parameters
        ----------
        var: str
            Unit name for spherical harmonics or spatial fields
        """
        try:
            return getattr(self, var)
        except Exception as exc:
            raise ValueError(f'Unknown units {var}') from exc

    @staticmethod
    def bycode(var: int) -> str:
        """
        Get the name for a units code

        Parameters
        ----------
        var: int
            Named unit code for spherical harmonics or spatial fields
        """
        named_units = [
            'norm', # 0: keep original scale
            'cmwe', # 1: cmwe, centimeters water equivalent
            'mmGH', # 2: mmGH, mm geoid height
            'mmCU', # 3: mmCU, mm elastic crustal deformation
            'microGal', # 4: microGal, microGal gravity perturbations
            'mbar', # 5: mbar, equivalent surface pressure
            'cmVCU' # 6: cmVCU, cm viscoelastic crustal uplift (GIA)
        ]
        try:
            return named_units[var]
        except Exception as exc:
            raise ValueError(f'Unknown units code {var}') from exc

    @staticmethod
    def get_attributes(var: str) -> str:
        """
        Get the variable attributes for a specific unit name

        Parameters
        ----------
        var: tuple
            Units and descriptive longname
        """
        named_attributes = dict(
            norm=('1', 'Fully-Normalized'),
            mmwe=('mm', 'Equivalent_Water_Thickness'),
            cmwe=('cm', 'Equivalent_Water_Thickness'),
            mmGH=('mm', 'Geoid_Height'),
            mmCU=('mm','Elastic_Crustal_Uplift'),
            mmCH=('mm','Horizontal_Elastic_Crustal_Deformation'),
            microGal=(u'\u03BCGal', 'Gravitational_Undulation'),
            mbar=('mbar', 'Equivalent_Surface_Pressure'),
            cmVCU=('cm', 'Viscoelastic_Crustal_Uplift'),
            mVCU=('meters', 'Viscoelastic_Crustal_Uplift')
        )
        try:
            return named_attributes[var]
        except Exception as exc:
            raise ValueError(f'Unknown units {var}') from exc

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
