#!/usr/bin/env python
u"""
geocenter.py
Written by Tyler Sutterley (06/2019)

Calculates the geocenter variation (in mm) from degree 1 Stokes Coefficients
Calculates the Degree 1 Stokes Coefficients of a geocenter variation (in mm)

CALLING SEQUENCE:
    xyz = geocenter(C10=C10, C11=C11, S11=S11)
    Ylms = geocenter(X=x, Y=y, Z=z, INVERSE=True)

OPTIONS:
    RADIUS: Earth's radius for calculating spherical harmonics from SLR data
    INVERSE: calculates the Stokes Coefficients from geocenter (True/False)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 06/2019: added option RADIUS to manually set the Earth's radius
    Updated 04/2017: changed option from INV to INVERSE and made True/False
    Updated 04/2015: calculate radius of the Earth directly in program
    Updated 02/2014: minor update to if statement
    Updated 03/2013: converted to python
"""
import numpy as np

def geocenter(C10=0,C11=0,S11=0,X=0,Y=0,Z=0,RADIUS=None,INVERSE=False):
    if RADIUS is None:
        #-- WGS84 ellipsoid parameters
        a_axis = 6378137.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.257223563#-- flattening of the ellipsoid
        #-- Mean Earth's Radius in mm having the same volume as WGS84 ellipsoid
        #-- (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1D -f)
        rad_e = 1000.0*a_axis*(1.0 -flat)**(1.0/3.0)
    else:
        rad_e = np.copy(RADIUS)
    #-- check if calculating spherical harmonics or geocenter coefficients
    if INVERSE:
        #-- inverse: geocenter to Stokes Coefficients
        C10 = Z/(rad_e*np.sqrt(3.0))
        C11 = X/(rad_e*np.sqrt(3.0))
        S11 = Y/(rad_e*np.sqrt(3.0))
    else:
        #-- Stokes Coefficients to geocenter
        Z = C10*rad_e*np.sqrt(3.0)
        X = C11*rad_e*np.sqrt(3.0)
        Y = S11*rad_e*np.sqrt(3.0)

    return {'C10':C10, 'C11':C11, 'S11':S11, 'x':X, 'y':Y, 'z':Z}
