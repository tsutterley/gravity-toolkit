#!/usr/bin/env python
u"""
test_point_masses.py (01/2023)
UPDATE HISTORY:
    Updated 01/2023: single implicit import of gravity toolkit
    Written 02/2021
"""
import pytest
import numpy as np
import gravity_toolkit as gravtk

# parameterize the number of point masses
@pytest.mark.parametrize("NPTS", np.random.randint(2,2000,size=1))
def test_point_masses(NPTS):
    # create spatial grid
    dlon,dlat = (1.0,1.0)
    lat = np.arange(90.0 - dlat/2.0, -90.0 - dlat/2.0, -dlat)
    lon = np.arange(-180.0 + dlon/2.0, 180.0 + dlon/2.0, dlon)
    gridlon,gridlat = np.meshgrid(lon,lat)
    nlat,nlon = np.shape(gridlon)

    # parameterize point masses
    LAT = lat[0]-dlat*np.random.randint(0,nlat,size=NPTS)
    LON  = lon[0]+dlon*np.random.randint(0,nlon,size=NPTS)
    MASS = 100.0 - 200.0*np.random.randn(NPTS)

    # create test gridded field
    data = np.zeros((nlat,nlon))
    for i in range(NPTS):
        indy,indx = np.nonzero((gridlat == LAT[i]) & (gridlon == LON[i]))
        data[indy,indx] += MASS[i]

    # path to load Love numbers file
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','love_numbers'])
    # read load Love numbers
    hl,kl,ll = gravtk.read_love_numbers(love_numbers_file)
    # calculate harmonics and degree amplitudes for each case
    grid_Ylms = gravtk.gen_stokes(data, lon, lat,
        LMAX=60, UNITS=2, LOVE=(hl,kl,ll))
    point_Ylms = gravtk.gen_point_load(MASS, LON, LAT,
        LMAX=60, UNITS=2, LOVE=(hl,kl,ll))

    # check that harmonic data is equal to machine precision
    difference_Ylms = grid_Ylms.copy()
    difference_Ylms.subtract(point_Ylms)
    harmonic_eps = np.finfo(np.float32).eps
    assert np.all(np.abs(difference_Ylms.clm) < harmonic_eps)
    # verify that the degree amplitudes are within tolerance
    difference_amplitude = grid_Ylms.amplitude - point_Ylms.amplitude
    assert np.all(np.abs(difference_amplitude) < harmonic_eps)
