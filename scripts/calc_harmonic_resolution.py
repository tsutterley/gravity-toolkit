#!/usr/bin/env python
u"""
calc_harmonic_resolution.py
Written by Tyler Sutterley (04/2022)

Calculates the spatial resolution that can be resolved
    by the spherical harmonics of a certain degree
Default method uses the smallest half-wavelength that can be resolved
    (is equal to approximately 20000/lmax km)
Secondary method calculates the smallest possible bump that can be resolved
    by dividing the area of a sphere by (lmax+1)^2

CALLING SEQUENCE:
    python calc_harmonic_resolution.py --lmax 60 --cap

COMMAND LINE OPTIONS:
    -l X, --lmax X: maximum degree of spherical harmonics
    -R X, --radius X: average radius of the Earth in kilometers
    -C, --cap: calculate the smallest possible bump that can be resolved

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2005)
        http://www.springerlink.com/content/978-3-211-33544-4
    Barthelmes, "Definition of Functionals of the Geopotential and Their
        Calculation from Spherical Harmonic Models", STR09/02 (2009)
        http://icgem.gfz-potsdam.de/str-0902-revised.pdf

UPDATE HISTORY:
    Updated 04/2022: use argparse descriptions within documentation
    Updated 09/2020: using argparse to set parameters
    Updated 10/2019: changing Y/N flags to True/False
    Updated 02/2014: minor update to if statement
    Updated 08/2013: changed SPH_CAP option to (Y/N)
    Written 01/2013
"""
import sys
import argparse
import numpy as np

# PURPOSE: Calculates minimum spatial resolution that can be resolved
# from spherical harmonics of a maximum degree
def calc_harmonic_resolution(LMAX, RADIUS=6371.0008, SPH_CAP=False):
    """
    Calculates minimum spatial resolution that can be resolved from
        spherical harmonics of a maximum degree

    Arguments
    ---------
    LMAX: maximum spherical harmonic degree

    Keyword arguments
    -----------------
    RADIUS: average radius of the Earth in kilometers
    SPH_CAP: calculate the smallest possible bump that can be resolved
    """
    if SPH_CAP:
        # Smallest diameter of a spherical cap that can be resolved by the
        # harmonics.  Size of the smallest bump, half-wavelength, which can
        # be produced by the clm/slm
        psi_min = 4.0*RADIUS*np.arcsin(1.0/(LMAX+1.0))
    else:
        # Shortest half-wavelength that can be resolved by the clm/slm
        # This estimation is based on the number of possible zeros along
        # the equator
        psi_min = np.pi*RADIUS/LMAX
    return psi_min

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lmax','-l', metavar='LMAX',
        type=int, nargs='+',
        help='maximum degree of spherical harmonics')
    parser.add_argument('--radius','-R',
        type=float, default=6371.0008,
        help='Average radius of the Earth in kilometers')
    parser.add_argument('--cap','-C',
        default=False, action='store_true',
        help='Calculate smallest possible bump that can be resolved')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()
    # for each entered spherical harmonic degree
    for LMAX in args.lmax:
        psi_min = calc_harmonic_resolution(LMAX,
            RADIUS=args.radius, SPH_CAP=args.cap)
        print('{0:5d}: {1:0.4f} km'.format(LMAX,psi_min))

# run main program
if __name__ == '__main__':
    main()
