#!/usr/bin/env python
u"""
test_harmonics.py (08/2020)
Tests harmonic programs using the Velicogna and Wahr (2013) Greenland synthetic
    1. Converts synthetic spatial distribution to spherical harmonics
    2. Compares output spherical harmonics with validation dataset
    3. Combines harmonics to calculate a truncated and smoothed spatial dataset
    4. Compares output smoothed spatial distribution with validation dataset
Tests harmonic objects flatten, expansion and iteration routines
"""
import os
import warnings
import pytest
import inspect
import numpy as np
import gravity_toolkit as gravtk
import matplotlib.pyplot as plt

# path to test files
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

# PURPOSE: test harmonic conversion programs
def test_harmonics():
    # path to load Love numbers file
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','love_numbers'])
    # read load Love numbers
    hl,kl,ll = gravtk.read_love_numbers(love_numbers_file)

    # read input spatial distribution file
    distribution_file = 'out.green_ice.grid.0.5.2008.cmh20.gz'
    input_distribution = gravtk.spatial(spacing=[0.5,0.5], nlat=361,
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
    PLM,dPLM = gravtk.associated_legendre(LMAX, np.cos(theta),
        method='mohlenkamp')
    # convert to spherical harmonics
    test_Ylms = gravtk.gen_stokes(input_distribution.data,
        input_distribution.lon, input_distribution.lat, UNITS=1, LMAX=LMAX,
        PLM=PLM, LOVE=(hl,kl,ll))

    # read harmonics from file
    harmonics_file = 'out.geoid.green_ice.0.5.2008.60.gz'
    valid_Ylms = gravtk.harmonics(lmax=LMAX, mmax=LMAX).from_ascii(
        os.path.join(filepath,harmonics_file),date=False,compression='gzip')

    # check that harmonic data is equal to machine precision
    difference_Ylms = test_Ylms.copy()
    difference_Ylms.subtract(valid_Ylms)
    harmonic_eps = np.finfo(np.float32).eps
    assert np.all(np.abs(difference_Ylms.clm) < harmonic_eps)
    assert np.all(np.abs(difference_Ylms.slm) < harmonic_eps)

    # cmwe, centimeters water equivalent
    dfactor = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
    wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
    test_Ylms.convolve(dfactor.cmwe*wt)
    # convert harmonics back to spatial domain at same grid spacing
    test_distribution = gravtk.harmonic_summation(test_Ylms.clm,
        test_Ylms.slm, input_distribution.lon, input_distribution.lat,
        LMAX=LMAX, PLM=PLM).T
    # convert harmonics using fast-fourier transform method
    test_transform = gravtk.harmonic_transform(test_Ylms.clm,
        test_Ylms.slm, input_distribution.lon, input_distribution.lat,
        LMAX=LMAX, PLM=PLM).T

    # read input and output spatial distribution files
    distribution_file = 'out.combine.green_ice.0.5.2008.60.gz'
    output_distribution = gravtk.spatial(spacing=[0.5,0.5], nlat=361,
        nlon=721, extent=[0,360.0,-90,90]).from_ascii(
        os.path.join(filepath,distribution_file),date=False,compression='gzip')

    # check that data is equal to machine precision
    difference_distribution = test_distribution - output_distribution.data
    difference_transform = test_transform - output_distribution.data
    distribution_eps = np.finfo(np.float16).eps
    assert np.all(np.abs(difference_distribution) < distribution_eps)
    assert np.all(np.abs(difference_transform) < distribution_eps)

# PURPOSE: test harmonic objects
def test_iterate():
    # maximum spherical harmonic degree and order
    LMAX = 60
    MMAX = 30
    # number of harmonics
    n_harm = (LMAX**2 + 3*LMAX - (LMAX-MMAX)**2 - (LMAX-MMAX))//2 + 1
    # number of time points for test
    nt = 12

    # create flattened test harmonics
    flat_Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    ll,mm = np.meshgrid(np.arange(LMAX+1), np.arange(MMAX+1))
    ii,jj = np.triu_indices(MMAX+1, k=0, m=LMAX+1)
    flat_Ylms.l = ll[ii,jj].astype(np.int64)
    flat_Ylms.m = mm[ii,jj].astype(np.int64)
    # add date variables
    flat_Ylms.month = 1 + np.arange(nt)
    flat_Ylms.time = 2002.0 + (flat_Ylms.month - 0.5)/12.0
    # create random harmonics
    flat_Ylms.clm = np.random.rand(n_harm,nt)
    flat_Ylms.slm = np.random.rand(n_harm,nt)
    # reshape harmonics object
    valid_Ylms = flat_Ylms.expand(date=True)

    # iterate over harmonics
    for i,test in enumerate(valid_Ylms):
        valid = valid_Ylms.index(i)
        assert np.isclose(valid.l, test.l).all()
        assert np.isclose(valid.m, test.m).all()
        assert np.isclose(valid.month, test.month).all()
        assert np.isclose(valid.time, test.time).all()
        assert np.isclose(valid.clm, test.clm).all()
        assert np.isclose(valid.slm, test.slm).all()

    # flatten harmonics objects
    reshaped_Ylms = valid_Ylms.flatten(date=True)
    assert np.isclose(reshaped_Ylms.l, flat_Ylms.l).all()
    assert np.isclose(reshaped_Ylms.m, flat_Ylms.m).all()
    assert np.isclose(reshaped_Ylms.month, flat_Ylms.month).all()
    assert np.isclose(reshaped_Ylms.time, flat_Ylms.time).all()
    assert np.isclose(reshaped_Ylms.clm, flat_Ylms.clm).all()
    assert np.isclose(reshaped_Ylms.slm, flat_Ylms.slm).all()
