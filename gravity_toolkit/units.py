#!/usr/bin/env python
u"""
units.py
Written by Tyler Sutterley (08/2020)

Class for converting GRACE/GRACE-FO Level-2 data to specific units

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 08/2020: made semi-major axis and ellipsoidal flattening arguments
    Updated 07/2020: added class docstring
    Updated 04/2020: include earth parameters as attributes
    Written 03/2020
"""
import os
import re
import numpy as np

class units(object):
    """
    Class for converting spherical harmonic and spatial data to specific units
    """
    np.seterr(invalid='ignore')
    def __init__(self, lmax=None, a_axis=6.378137e8, flat=1.0/298.257223563):
        # Earth Parameters
        # Average Density of the Earth [g/cm^3]
        self.rho_e=5.517
        # Semi-major axis of the ellipsoid [cm]
        self.a_axis=a_axis
        # Flattening of the ellipsoid
        self.flat=flat
        # Semi-minor axis of the WGS84 ellipsoid [cm]
        self.b_axis=(1.0-self.flat)*self.a_axis
        # Average Radius of the Earth having the same volume as WGS84 [cm]
        self.rad_e=self.a_axis*(1.0-self.flat)**(1.0/3.0)
        # standard gravitational acceleration (World Meteorological Organization)
        self.g_wmo=980.665
        # Degree dependent parameters
        self.norm=None
        self.cmwe=None
        self.mmwe=None
        self.mmGH=None
        self.mmCU=None
        self.mmCH=None
        self.microGal=None
        self.mbar=None
        self.Pa=None
        self.lmax=lmax
        self.l=np.arange(self.lmax+1) if self.lmax else None

    def harmonic(self, hl, kl, ll):
        """
        Calculates degree dependent factors for harmonic units
        Inputs: hl, kl, ll load love numbers to degree lmax
        """
        # Earth Parameters
        # WGS84 Gravitational Constant of the Earth [cm^3/s^2]
        GM_e=3986004.418e14
        # Gravitational Constant of the Earth's atmosphere
        GM_atm=3.5e14
        # Gravitational Constant of the Earth (w/o atm)
        GM=GM_e-GM_atm
        # degree dependent coefficients
        # norm, geodesy normalized spherical harmonics
        self.norm=np.ones((self.lmax+1))
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe=self.rho_e*self.rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/3.0
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe=10.0*self.rho_e*self.rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/3.0
        # mmGH, millimeters geoid height
        self.mmGH=np.ones((self.lmax+1))*(10.0*self.rad_e)
        # mmCU, millimeters elastic crustal deformation (uplift)
        self.mmCU=10.0*self.rad_e*hl[self.l]/(1.0+kl[self.l])
        # mmCH, millimeters elastic crustal deformation (horizontal)
        self.mmCH=10.0*self.rad_e*ll[self.l]/(1.0+kl[self.l])
        # cmVCU, centimeters viscoelastic crustal uplift
        self.cmVCU=self.rad_e*(2.0*self.l+1.0)/2.0
        # mVCU, meters viscoelastic crustal uplift
        self.mVCU=self.rad_e*(2.0*self.l+1.0)/200.0
        # microGal, microGal gravity perturbations
        self.microGal=1.e6*GM*(self.l+1.0)/(self.rad_e**2.0)
        # mbar, millibar equivalent surface pressure
        self.mbar=self.g_wmo*self.rho_e*self.rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/3e3
        # Pa, pascals equivalent surface pressure
        self.Pa=self.g_wmo*self.rho_e*self.rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/30.0
        # return the degree dependent unit conversions
        return self

    def spatial(self, hl, kl, ll):
        """
        Calculates degree dependent factors for converting spatial units
        Inputs: hl, kl, ll load love numbers to degree lmax
        """
        # degree dependent coefficients
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe=3.0*(1.0+kl[self.l])/(1.0+2.0*self.l)/(4.0*np.pi*self.rad_e*self.rho_e)
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe=3.0*(1.0+kl[self.l])/(1.0+2.0*self.l)/(40.0*np.pi*self.rad_e*self.rho_e)

        # return the degree dependent unit conversions
        return self
