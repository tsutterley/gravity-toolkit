#!/usr/bin/env python
u"""
units.py
Written by Tyler Sutterley (03/2020)

Class for converting GRACE/GRACE-FO Level-2 data to specific units

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
    Written 03/2020
"""
import os
import re
import numpy as np

class units(object):
    np.seterr(invalid='ignore')
    def __init__(self, lmax=None):
        self.norm = None
        self.cmwe = None
        self.mmwe = None
        self.mmGH = None
        self.mmCU = None
        self.mmCH = None
        self.microGal = None
        self.mbar = None
        self.Pa = None
        self.lmax=lmax
        self.l=np.arange(self.lmax+1) if self.lmax else None

    def harmonic(self, hl, kl, ll):
        """
        Calculates degree dependent factors for harmonic units
        Inputs: hl, kl, ll load love numbers to degree lmax
        """
        # Earth Parameters
        # Average Density of the Earth [g/cm^3]
        rho_e = 5.517
        # Average Radius of the Earth [cm]
        rad_e = 6.371e8
        # WGS84 Gravitational Constant of the Earth [cm^3/s^2]
        GM_e = 3986004.418e14
        # Gravitational Constant of the Earth's atmosphere
        GM_atm = 3.5e14
        # Gravitational Constant of the Earth (w/o atm)
        GM = GM_e - GM_atm
        # standard gravitational acceleration (World Meteorological Organization)
        g_wmo = 980.665
        # degree dependent coefficients
        # norm, geodesy normalized spherical harmonics
        self.norm=np.ones((self.lmax+1))
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe=rho_e*rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/3.0
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe=rho_e*rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/30.0
        # mmGH, millimeters geoid height
        self.mmGH=np.ones((self.lmax+1))*(10.0*rad_e)
        # mmCU, millimeters elastic crustal deformation (uplift)
        self.mmCU=10.0*rad_e*hl[self.l]/(1.0+kl[self.l])
        # mmCH, millimeters elastic crustal deformation (horizontal)
        self.mmCH=10.0*rad_e*ll[self.l]/(1.0+kl[self.l])
        # microGal, microGal gravity perturbations
        self.microGal=1.e6*GM*(self.l+1.0)/(rad_e**2.0)
        # mbar, millibar equivalent surface pressure
        self.mbar=g_wmo*rho_e*rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/3e3
        # Pa, pascals equivalent surface pressure
        self.Pa=g_wmo*rho_e*rad_e*(2.0*self.l+1.0)/(1.0+kl[self.l])/30.0
        # return the degree dependent unit conversions
        return self

    def spatial(self, hl, kl, ll):
        """
        Calculates degree dependent factors for converting spatial units
        Inputs: hl, kl, ll load love numbers to degree lmax
        """
        # Earth Parameters
        # Average Density of the Earth [g/cm^3]
        rho_e = 5.517
        # Average Radius of the Earth [cm]
        rad_e = 6.371e8
        # degree dependent coefficients
        # cmwe, centimeters water equivalent [g/cm^2]
        self.cmwe=3.0*(1.0+kl[self.l])/(1.0+2.0*self.l)/(4.0*np.pi*rad_e*rho_e)
        # mmwe, millimeters water equivalent [kg/m^2]
        self.mmwe=3.0*(1.0+kl[self.l])/(1.0+2.0*self.l)/(40.0*np.pi*rad_e*rho_e)

        # return the degree dependent unit conversions
        return self
