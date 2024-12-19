#!/usr/bin/env python
u"""
mascons.py
Written by Tyler Sutterley (03/2023)
Conversion routines for publicly available GRACE/GRACE-FO mascon solutions

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

REFERENCES:
    grid2mascon.m written by Felix Landerer and David Wiese (JPL)
    mascon2grid.m written by Felix Landerer and David Wiese (JPL)

UPDATE HISTORY:
    Updated 03/2023: improve typing for variables in docstrings
    Updated 11/2022: use lowercase keyword arguments
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: publicly released version
    Updated 09/2019: check number of latitude points for using regional grids
    Updated 08/2018: add GSFC grid and mascon conversion routines
        extract number of mascons and number of variables
        use lat_bound and lon_bound variables as inputs to JPL conversion
    Updated 06/2016: check that mat file exists before attempting conversion
    Updated 03/2016: minor clean up of mscdata calculation, update comments
    Updated 12/2015: added TRANSPOSE option to output spatial routines
    Written 07/2013
"""
import copy
import warnings
import numpy as np

def to_gsfc(gdata, lon, lat, lon_center, lat_center, lon_span, lat_span):
    """
    Converts an input gridded field to an output GSFC mascon array
    :cite:p:`Luthcke:2013ep`

    Parameters
    ----------
    gdata: np.ndarray
        gridded data
    lon: np.ndarray
        column vector of defined longitude points
    lat: np.ndarray
        column vector of defined latitude points
    lon_center: np.ndarray
        mascon longitudinal center points
    lat_center: np.ndarray
        mascon latitudinal center points
    lon_span: np.ndarray
        mascon longitudinal central angles
    lat_span: np.ndarray
        mascon latitudinal central angles

    Returns
    -------
    data: np.ndarray
        row vector of mascons
    lat_center: np.ndarray
        row vector of latitude values for mascon centers
    lon_center: np.ndarray
        row vector of longitude values for mascon centers
    """
    # number of mascons
    nmas = len(lon_center)
    # convert mascon centers to -180:180
    gt180, = np.nonzero(lon_center > 180)
    lon_center[gt180] -= 360.0
    # remove singleton dimensions
    lat = np.squeeze(lat)
    lon = np.squeeze(lon)
    # for mascons centered on 180: use 0:360
    alon = np.copy(lon)
    lt0, = np.nonzero(lon < 0)
    alon[lt0] += 360.0

    # loop over each mascon bin and average gdata with the cos-lat weights
    # for that bin
    mascon_array = {}
    mascon_array['data'] = np.zeros((nmas))
    mascon_array['lon_center'] = np.zeros((nmas))
    mascon_array['lat_center'] = np.zeros((nmas))
    for k in range(0,nmas):
        # create latitudinal and longitudinal bounds for mascon k
        if (lat_center[k] == 90.0) | (lat_center[k] == -90.0):
            # NH and SH polar mascons
            lon_bound = [0.0,360.0]
            lat_bound = lat_center[k] + np.array([-1.0,1.0])*lat_span[k]
        else:
            # convert from mascon centers to mascon bounds
            lon_bound = lon_center[k] + np.array([-0.5,0.5])*lon_span[k]
            lat_bound = lat_center[k] + np.array([-0.5,0.5])*lat_span[k]
        # if mascon is centered on +/-180: use 0:360
        if ((lon_bound[0] <= 180.0) & (lon_bound[1] >= 180.0)):
            ilon = alon.copy()
        elif ((lon_bound[0] <= -180.0) & (lon_bound[1] >= -180.0)):
            lon_bound += 360.0
            ilon = alon.copy()
        else:
            ilon = lon.copy()
        # indices for grid points within the mascon
        I, = np.nonzero((lat >= lat_bound[0]) & (lat < lat_bound[1]))
        J, = np.nonzero((ilon >= lon_bound[0]) & (ilon < lon_bound[1]))
        I,J = (I[np.newaxis,:], J[:,np.newaxis])
        # calculate average data for mascon bin
        mascon_array['data'][k] = np.mean((np.cos(lat[I]*np.pi/180.0) /
            np.mean(np.cos(lat[I]*np.pi/180.0)))*gdata[I,J]/len(I))
        mascon_array['lat_center'][k] = lat_center[k]
        mascon_array['lon_center'][k] = lon_center[k]

    # return python dictionary with the mascon array data, lon and lat
    return mascon_array

def to_jpl(gdata, lon, lat, lon_bound, lat_bound):
    """
    Converts an input gridded field to an output JPL mascon array
    :cite:p:`Watkins:2015jl`

    Parameters
    ----------
    gdata: np.ndarray
        gridded data
    lon: np.ndarray
        column vector of defined longitude points
    lat: np.ndarray
        column vector of defined latitude points
    lon_bound: np.ndarray
        mascon longitudinal bounds from coordinate file
    lat_bound: np.ndarray
        mascon latitudinal bounds from coordinate file

    Returns
    -------
    data: np.ndarray
        row vector of mascons
    mask: np.ndarray
        row vector of mask values showing if mascon has no data
    lat: np.ndarray
        row vector of latitude values for mascons
    lon: np.ndarray
        row vector of longitude values for mascons
    """
    # mascon dimensions
    nmas,nvar = lat_bound.shape
    # remove singleton dimensions
    lat = np.squeeze(lat)
    lon = np.squeeze(lon)

    # loop over each mascon bin and average gdata with the cos-lat weights
    # for that bin
    mascon_array = {}
    mascon_array['data'] = np.zeros((nmas))
    mascon_array['mask'] = np.zeros((nmas),dtype=bool)
    mascon_array['lon'] = np.zeros((nmas))
    mascon_array['lat'] = np.zeros((nmas))
    for k in range(0,nmas):
        # indices for grid points within the mascon
        I, = np.nonzero((lat >= lat_bound[k,1]) & (lat < lat_bound[k,0]))
        J, = np.nonzero((lon >= lon_bound[k,0]) & (lon < lon_bound[k,2]))
        nlt = np.count_nonzero((lat >= lat_bound[k,1]) & (lat < lat_bound[k,0]))
        I,J = (I[np.newaxis,:], J[:,np.newaxis])
        # calculate average data for mascon bin
        mascon_array['data'][k] = np.mean((np.cos(lat[I]*np.pi/180.0) /
            np.mean(np.cos(lat[I]*np.pi/180.0)))*gdata[I,J]/nlt)
        # calculate coordinates of mascon center
        mascon_array['lat'][k] = (lat_bound[k,1]+lat_bound[k,0])/2.0
        mascon_array['lon'][k] = (lon_bound[k,1]+lon_bound[k,2])/2.0
        mascon_array['mask'][k] = bool(nlt == 0)
        # Do a check at the poles to make the lat/lon equal to +/-90/0
        if (np.abs(lat_bound[k,0]) == 90):
            mascon_array['lat'][k] = lat_bound[k,0]
            mascon_array['lon'][k] = 0.0
        if (np.abs(lat_bound[k,1]) == 90):
            mascon_array['lat'][k] = lat_bound[k,1]
            mascon_array['lon'][k] = 0.0
    # replace invalid data with 0
    mascon_array['data'][mascon_array['mask']] = 0.0
    # return python dictionary with the mascon array data, lon and lat
    return mascon_array

def from_gsfc(mscdata, grid_spacing, lon_center, lat_center, lon_span, lat_span,
    **kwargs):
    """
    Converts an input GSFC mascon array to an output gridded field
    :cite:p:`Luthcke:2013ep`

    Parameters
    ----------
    mscdata: np.ndarray
        row vector of mascons
    grid_spacing: np.ndarray
        spacing of the lat/lon grid
    lon_center: float
        mascon np.ndarray center points
    lat_center: np.ndarray
        mascon latitudinal center points
    lon_span: np.ndarray
        mascon longitudinal central angles
    lat_span: np.ndarray
        mascon latitudinal central angles
    transpose: bool, default False
        transpose output matrix (lon/lat)

    Returns
    -------
    mdata: np.ndarray
        distributed mass grid
    """
    # set default keyword arguments
    kwargs.setdefault('transpose', False)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TRANSPOSE='transpose')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # number of mascons
    nmas = len(lon_center)
    # convert mascon centers to -180:180
    gt180, = np.nonzero(lon_center > 180)
    lon_center[gt180] -= 360.0

    # Define output latitude and longitude grids
    lon = np.arange(-180.0+grid_spacing/2.0,180.0+grid_spacing/2.0,grid_spacing)
    lat = np.arange(90.0-grid_spacing/2.0,-90.0-grid_spacing/2.0,-grid_spacing)
    nlon, nlat = (len(lon),len(lat))
    # for mascons centered on 180: use 0:360
    alon = np.copy(lon)
    lt0, = np.nonzero(lon < 0)
    alon[lt0] += 360.0

    # loop over each mascon bin and assign value to grid points inside bin:
    mdata = np.zeros((nlat, nlon))
    for k in range(0, nmas):
        # create latitudinal and longitudinal bounds for mascon k
        if (lat_center[k] == 90.0) | (lat_center[k] == -90.0):
            # NH and SH polar mascons
            lon_bound = [0.0,360.0]
            lat_bound = lat_center[k] + np.array([-1.0,1.0])*lat_span[k]
        else:
            # convert from mascon centers to mascon bounds
            lon_bound = lon_center[k] + np.array([-0.5,0.5])*lon_span[k]
            lat_bound = lat_center[k] + np.array([-0.5,0.5])*lat_span[k]
        # if mascon is centered on +/-180: use 0:360
        if ((lon_bound[0] <= 180.0) & (lon_bound[1] >= 180.0)):
            ilon = alon.copy()
        elif ((lon_bound[0] <= -180.0) & (lon_bound[1] >= -180.0)):
            lon_bound += 360.0
            ilon = alon.copy()
        else:
            ilon = lon.copy()
        # indices for grid points within the mascon
        I, = np.nonzero((lat >= lat_bound[0]) & (lat < lat_bound[1]))
        J, = np.nonzero((ilon >= lon_bound[0]) & (ilon < lon_bound[1]))
        I,J = (I[np.newaxis,:], J[:,np.newaxis])
        mdata[I,J] = mscdata[k]

    # return array
    if kwargs['transpose']:
        return mdata.T
    else:
        return mdata

def from_jpl(mscdata, grid_spacing, lon_bound, lat_bound, **kwargs):
    """
    Converts an input JPL mascon array to an output gridded field
    :cite:p:`Watkins:2015jl`

    Parameters
    ----------
    mscdata: np.ndarray
        row vector of mascons
    grid_spacing: np.ndarray
        spacing of lat/lon grid
    lon_bound: np.ndarray
        mascon longitudinal bounds from coordinate file
    lat_bound: np.ndarray
        mascon latitudinal bounds from coordinate file
    transpose: bool, default False
        transpose output matrix (lon/lat)

    Returns
    -------
    mdata: np.ndarray
        distributed mass grid
    """
    # set default keyword arguments
    kwargs.setdefault('transpose', False)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TRANSPOSE='transpose')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # mascon dimensions
    nmas,nvar = lat_bound.shape

    # Define latitude and longitude grids
    # output lon will not include 360
    # output lat will not include 90
    lon = np.arange(grid_spacing/2.0, 360.0+grid_spacing/2.0, grid_spacing)
    lat = np.arange(-90.0+grid_spacing/2.0, 90.0+grid_spacing/2.0, grid_spacing)
    nlon, nlat = (len(lon),len(lat))

    # loop over each mascon bin and assign value to grid points inside bin:
    mdata = np.zeros((nlat, nlon))
    for k in range(0, nmas):
        I, = np.nonzero((lat >= lat_bound[k,1]) & (lat < lat_bound[k,0]))
        J, = np.nonzero((lon >= lon_bound[k,0]) & (lon < lon_bound[k,2]))
        I,J = (I[np.newaxis,:], J[:,np.newaxis])
        mdata[I,J] = mscdata[k]

    # return array
    if kwargs['transpose']:
        return mdata.T
    else:
        return mdata
