#!/usr/bin/env python
u"""
ocean_stokes.py
Written by Tyler Sutterley (04/2022)

Reads a land-sea mask and converts to a series of spherical harmonics

INPUTS:
    LANDMASK: Mask file to use as input following Sutterley et al. (2020)
        1x1 degree mask distributed from UCAR as part of NCL
            https://www.ncl.ucar.edu/Applications/Data/cdf/landsea.nc
        1,0.5 and 0.25 degree masks distributed from ORNL as part of ISLSCP
            https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html
    LMAX: maximum spherical harmonic degree of the output harmonics

OPTIONS:
    MMAX: maximum spherical harmonic order of the output harmonics
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)
    VARNAME: variable name for mask in netCDF4 file
    SIMPLIFY: simplify land mask by removing isolated points

OUTPUTS:
    clm: Cosine spherical harmonic coefficients (geodesy normalization)
    slm: Sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org

PROGRAM DEPENDENCIES
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    spatial.py: spatial data class for reading, writing and processing data

REFERENCES:
    T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self-Consistent Ice Mass
        Balance and Regional Sea Level From Time-Variable Gravity",
    Earth and Space Science, 7, 2020. https://doi.org/10.1029/2019EA000860

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: added option VARNAME to generalize input masks
    Updated 12/2020: added simplify function to remove isolated points
    Updated 07/2020: added function docstrings
    Updated 06/2020: using spatial data class for input and output operations
    Updated 04/2020 for public release
    Updated 04/2020: added option LOVE for passing load love numbers
    Updated 09/2017: added option MASK for different input masks
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Written 03/2015
"""
import numpy as np
from gravity_toolkit.spatial import spatial
from gravity_toolkit.gen_stokes import gen_stokes

def ocean_stokes(LANDMASK, LMAX, MMAX=None, LOVE=None, VARNAME='LSMASK',
    SIMPLIFY=False):
    """
    Converts data from spherical harmonic coefficients to a spatial field

    Parameters
    ----------
    LANDMASK: str
        netCDF4 land-sea mask file
    LMAX: int
        maximum spherical harmonic degree
    MMAX: int or NoneType, default None
        maximum spherical harmonic order of the output harmonics
    LOVE: tuple
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    VARNAME: str, default 'LSMASK'
        variable name for mask in netCDF4 file
    SIMPLIFY: bool, default False
        simplify land mask by removing isolated points

    Returns
    -------
    clm: float
        cosine spherical harmonic coefficients
    slm: float
        sine spherical harmonic coefficients
    l: int
        spherical harmonic degree to LMAX
    m: int
        spherical harmonic order to MMAX
    """
    # maximum spherical harmonic order
    MMAX = np.copy(LMAX) if MMAX is None else MMAX
    # Read Land-Sea Mask of specified input file
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    landsea = spatial().from_netCDF4(LANDMASK,
        date=False, varname=VARNAME)
    # create land function
    nth,nphi = landsea.shape
    land_function = np.zeros((nth,nphi),dtype=np.float64)
    # combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indx,indy] = 1.0
    # remove isolated points if specified
    if SIMPLIFY:
        land_function -= find_isolated_points(land_function)
    # ocean function reciprocal of land function
    ocean_function = 1.0 - land_function
    # convert to spherical harmonics (1 cm w.e.)
    ocean_Ylms = gen_stokes(ocean_function.T,landsea.lon,landsea.lat,
        UNITS=1,LMIN=0,LMAX=LMAX,MMAX=MMAX,LOVE=LOVE)
    # return the spherical harmonic coefficients
    return ocean_Ylms

def find_isolated_points(mask):
    """
    Simplify mask by removing isolated points

    Parameters
    ----------
    mask: int
        land-sea mask

    Returns
    -------
    isolated: int
        simplified land-sea mask
    """
    nth,_ = mask.shape
    laplacian = -4.0*np.copy(mask)
    laplacian += mask*np.roll(mask,1,axis=1)
    laplacian += mask*np.roll(mask,-1,axis=1)
    temp = np.roll(mask,1,axis=0)
    temp[0,:] = mask[1,:]
    laplacian += mask*temp
    temp = np.roll(mask,-1,axis=0)
    temp[nth-1,:] = mask[nth-2,:]
    laplacian += mask*temp
    # create mask of isolated points
    isolated = np.where(np.abs(laplacian) >= 3, 1, 0)
    return isolated