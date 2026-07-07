#!/usr/bin/env python
u"""
test_sea_level.py (07/2026)
"""
import pytest
import inspect
import pathlib
import numpy as np
import gravity_toolkit as gravtk

# path to test files
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: test sea level equation programs
def test_sea_level():
    # path to load Love numbers file
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','love_numbers'])
    # read load Love numbers
    LOVE = gravtk.read_love_numbers(love_numbers_file, FORMAT='class')

    # read land function file
    LANDMASK = filepath.joinpath('land.fcn.1_deg.gz')
    landsea = gravtk.spatial().from_ascii(LANDMASK,
        date=False, spacing=[1.0, 1.0], nlat=180, nlon=360,
        extent=[0.5,359.5,-89.5,89.5], compression='gzip')

    # spherical harmonic parameters
    # maximum spherical harmonic degree
    LMAX = 60
    # read harmonics from file
    harmonics_file = filepath.joinpath('out.geoid.green_ice.0.5.2008.60.gz')
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=LMAX).from_ascii(
        harmonics_file, date=False, compression='gzip')
    # calculate the legendre functions using Martin Mohlenkamp's relation
    th = np.radians(90.0 - landsea.lat)
    PLM, dPLM = gravtk.plm_mohlenkamp(LMAX, np.cos(th))

    # run pseudo-spectral sea level equation solver
    sea_level = gravtk.sea_level_equation(Ylms.clm, Ylms.slm,
        landsea.lon, landsea.lat, landsea.data.T, LMAX=LMAX,
        LOVE=LOVE, BODY_TIDE_LOVE=0,
        FLUID_LOVE=0, DENSITY=1.0, POLAR=True,
        PLM=PLM, ITERATIONS=2, FILL_VALUE=np.nan).T
    
    # check that sea level data is equal to file precision
    valid_file = filepath.joinpath('out.slf.green_ice.1_deg.2008.60.gz')
    validation = gravtk.spatial().from_ascii(valid_file,
        date=False, spacing=[1.0, 1.0], nlat=180, nlon=360,
        extent=[0.5,359.5,-89.5,89.5], compression='gzip')
    # check differences
    difference = validation.data - sea_level
    valid_difference = difference[np.isfinite(difference)]
    assert np.all(np.abs(valid_difference) < 1e-8)
