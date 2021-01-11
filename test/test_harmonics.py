#!/usr/bin/env python
u"""
test_harmonics.py (08/2020)
Tests harmonic programs using the Velicogna and Wahr (2013) Greenland synthetic
    1. Converts synthetic spatial distribution to spherical harmonics
    2. Compares output spherical harmonics with validation dataset
    3. Combines harmonics to calculate a truncated and smoothed spatial dataset
    4. Compares output smoothed spatial distribution with validation dataset
"""
import os
import warnings
import pytest
import inspect
import numpy as np
import gravity_toolkit.read_love_numbers
import gravity_toolkit.plm_mohlenkamp
import gravity_toolkit.gen_stokes
import gravity_toolkit.harmonic_summation
import gravity_toolkit.harmonics
import gravity_toolkit.spatial
from gravity_toolkit.utilities import get_data_path

def test_harmonics():
    # path to load Love numbers file
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    love_numbers_file = get_data_path(['data','love_numbers'])
    # read load Love numbers
    hl,kl,ll = gravity_toolkit.read_love_numbers(love_numbers_file)

    # read input spatial distribution file
    distribution_file = 'out.green_ice.grid.0.5.2008.cmh20.gz'
    input_distribution = gravity_toolkit.spatial(spacing=[0.5,0.5], nlat=361,
        nlon=721, extent=[0,360.0,-90,90]).from_ascii(
        os.path.join(filepath,distribution_file),date=False,compression='gzip')

    # spherical harmonic parameters
    # maximum spherical harmonic degree
    LMAX = 60
    # gaussian smoothing radius (km)
    RAD = 250.0

    # calculate colatitudes of input distribution
    theta = (90.0 - input_distribution.lat)*np.pi/180.0
    # use fortran thresholds for colatitude bounds
    theta[theta > np.arccos(-0.9999999)] = np.arccos(-0.9999999)
    theta[theta < np.arccos(0.9999999)] = np.arccos(0.9999999)
    # calculate Legendre polynomials with Martin Mohlenkamp's relation
    PLM = gravity_toolkit.plm_mohlenkamp(LMAX, np.cos(theta))
    # convert to spherical harmonics
    test_Ylms = gravity_toolkit.gen_stokes(input_distribution.data,
        input_distribution.lon, input_distribution.lat, UNITS=1, LMAX=LMAX,
        PLM=PLM, LOVE=(hl,kl,ll))

    # read harmonics from file
    harmonics_file = 'out.geoid.green_ice.0.5.2008.60.gz'
    valid_Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=LMAX).from_ascii(
        os.path.join(filepath,harmonics_file),date=False,compression='gzip')

    # check that harmonic data is equal to machine precision
    difference_Ylms = test_Ylms.copy()
    difference_Ylms.subtract(valid_Ylms)
    harmonic_eps = np.finfo(np.float32).eps
    assert np.all(np.abs(difference_Ylms.clm) < harmonic_eps)
    assert np.all(np.abs(difference_Ylms.slm) < harmonic_eps)

    # cmwe, centimeters water equivalent
    dfactor = gravity_toolkit.units(lmax=LMAX).harmonic(hl,kl,ll)
    wt = 2.0*np.pi*gravity_toolkit.gauss_weights(RAD,LMAX)
    test_Ylms.convolve(dfactor.cmwe*wt)
    # convert harmonics back to spatial domain at same grid spacing
    test_distribution = gravity_toolkit.harmonic_summation(test_Ylms.clm,
        test_Ylms.slm, input_distribution.lon, input_distribution.lat,
        LMAX=LMAX, PLM=PLM).T

    # read input and output spatial distribution files
    distribution_file = 'out.combine.green_ice.0.5.2008.60.gz'
    output_distribution = gravity_toolkit.spatial(spacing=[0.5,0.5], nlat=361,
        nlon=721, extent=[0,360.0,-90,90]).from_ascii(
        os.path.join(filepath,distribution_file),date=False,compression='gzip')

    # check that data is equal to machine precision
    difference_distribution = test_distribution - output_distribution.data
    distribution_eps = np.finfo(np.float16).eps
    assert np.all(np.abs(difference_distribution) < distribution_eps)
