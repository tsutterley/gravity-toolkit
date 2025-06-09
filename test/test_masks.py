#!/usr/bin/env python
u"""
test_masks.py (06/2025)
Tests that the stored masks can be read
"""
import pytest
import numpy as np
import gravity_toolkit as gravtk

# parameterize land-sea masks
_land_sea_masks = ["landsea_1d.nc", "landsea_hd.nc", "landsea_qd.nc"]
@pytest.mark.parametrize("LANDMASK", _land_sea_masks)
def test_lsmask(LANDMASK):
    u"""
    Test that the land mask can be read and that the
    area of the ocean is close to the expected value
    """
    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    lsmask = gravtk.utilities.get_data_path(['data',LANDMASK])
    landsea = gravtk.spatial().from_netCDF4(lsmask, date=False, varname='LSMASK')
    # degree spacing and grid dimensions
    dlon, dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # colatitude in radians
    gridlon, gridlat = np.meshgrid(landsea.lon, landsea.lat)
    th = (90.0 - gridlat)*np.pi/180.0
    # grid spacing in radians
    dphi = np.pi*np.abs(dlon)/180.0
    dth = np.pi*np.abs(dlat)/180.0
    # create land function
    land_function = np.zeros((nlat, nlon), dtype=np.float64)
    # combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indx,indy] = 1.0
    # calculate the ocean function
    ocean_function = 1.0 - land_function
    # average Radius of the Earth [km]
    rad_e = gravtk.units().rad_e/1e5
    # total area of ocean calculated by integrating the ocean function
    area = np.sum(ocean_function*np.sin(th)*dphi*dth*rad_e**2)
    # assert that the area is close to the expected value
    assert np.isclose(area, 3.62e8, rtol=0.01)

def test_smoothed():
    """
    Test that the smoothed land mask can be read
    """
    # Read Smoothed Ocean and Land Functions
    # will mask out land regions in the final current maps
    LANDMASK = gravtk.utilities.get_data_path(['data','land_fcn_300km.nc'])
    landsea = gravtk.spatial().from_netCDF4(LANDMASK,
        date=False, varname='LSMASK')
    assert landsea.shape == (180, 360)
    assert landsea.spacing == (1.0, 1.0)
