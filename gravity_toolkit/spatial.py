#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (06/2020)

Data class for reading, writing and processing spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (https://www.h5py.org/)

PROGRAM DEPENDENCIES:
    ncdf_write.py: writes output spatial data to COARDS-compliant netCDF4
    hdf5_write.py: writes output spatial data to HDF5
    ncdf_read.py: reads spatial data from COARDS-compliant netCDF4
    hdf5_read.py: reads spatial data from HDF5

UPDATE HISTORY:
    Written 06/2020
"""
import os
import re
import numpy as np
from gravity_toolkit.ncdf_write import ncdf_write
from gravity_toolkit.hdf5_write import hdf5_write
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.hdf5_read import hdf5_read

class spatial(object):
    np.seterr(invalid='ignore')
    def __init__(self, spacing=[None,None], nlat=None, nlon=None):
        self.data=None
        self.mask=None
        self.lon=None
        self.lat=None
        self.time=None
        self.month=None
        self.fill_value=None
        self.extent=[np.inf,-np.inf,np.inf,-np.inf]
        self.shape=[nlat,nlon,None]
        self.ndim=None
        self.spacing=spacing
        self.filename=None

    def from_ascii(self, filename, verbose=False, date=True,
        columns=['lon','lat','data','time']):
        """
        Read a spatial object from an ascii file
        Inputs: full path of input ascii file
        Options: ascii file contains date information
        """
        self.filename = filename
        print(self.filename) if verbose else None
        #-- read input ascii file (.txt) and split lines
        with open(os.path.expanduser(filename),'r') as f:
            file_contents = f.read().splitlines()
        #-- compile regular expression operator for extracting numerical values
        #-- from input ascii files of spatial data
        regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        #-- output spatial data
        self.lat = np.zeros((self.shape[0]))
        self.lon = np.zeros((self.shape[1]))
        self.data = np.zeros((self.shape[0],self.shape[1]))
        self.mask = np.zeros((self.shape[0],self.shape[1]),dtype=np.bool)
        #-- remove time from list of column names if not date
        columns = [c for c in columns if (c != 'time')]
        #-- extract spatial data array and convert to matrix
        #-- for each line in the file
        for line in file_contents:
            #-- extract columns of interest and assign to dict
            d = {c:r for c,r in zip(columns,rx.findall(line))}
            #-- convert line coordinates to integers
            ilon = np.int(np.float(d['lon'])/self.spacing[0])
            ilat = np.int((90.0-np.float(d['lat']))//self.spacing[1])
            #-- convert fortran exponentials if applicable
            self.data[ilat,ilon] = np.float(d['data'].replace('D','E'))
            self.lon[ilon] = np.float(d['lon'])
            self.lat[ilat] = np.float(d['lat'])
            #-- if the ascii file contains date variables
            if date:
                self.time = np.float(d['time'])
                self.month = np.int(12.0*(self.time - 2002.0)) + 1
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def from_netCDF4(self, filename, date=True, verbose=False,
        varname='z', lonname='lon', latname='lat'):
        """
        Read a spatial object from a netCDF4 file
        Inputs: full path of input netCDF4 file
        Options: netCDF4 file contains date information
        """
        self.filename = filename
        data = ncdf_read(os.path.expanduser(filename), VERBOSE=verbose,
            ATTRIBUTES=False, DATE=date, VARNAME=varname,
            LONNAME=lonname, LATNAME=latname, TIMENAME='time')
        self.data = data['data'].copy()
        self.lon = data['lon'].copy()
        self.lat = data['lat'].copy()
        if date:
            self.time = data['time'].copy()
            self.month = np.int(12.0*(self.time - 2002.0)) + 1
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def from_HDF5(self, filename, date=True, verbose=False,
        varname='z', lonname='lon', latname='lat'):
        """
        Read a spatial object from a HDF5 file
        Inputs: full path of input HDF5 file
        Options: HDF5 file contains date information
        """
        self.filename = filename
        data = hdf5_read(os.path.expanduser(filename), VERBOSE=verbose,
            ATTRIBUTES=False, DATE=date, VARNAME=varname,
            LONNAME=lonname, LATNAME=latname, TIMENAME='time')
        self.data = data['data'].copy()
        self.lon = data['lon'].copy()
        self.m = data['lat'].copy()
        if date:
            self.time = data['time'].copy()
            self.month = np.int(12.0*(self.time - 2002.0)) + 1
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def from_index(self, filename, format=None, date=True, sort=True):
        """
        Read a spatial object from an index of netCDF4 or HDF5 files
        Inputs: full path of index file to be read into a spatial object
        Options:
            format of files in index ( netCDF4 or HDF5)
            netCDF4, or HDF5 contains date information
            sort spatial objects by date information
        """
        self.filename = filename
        #-- Read index file of input spatial data
        with open(os.path.expanduser(filename),'r') as f:
            file_list = f.read().splitlines()
        #-- create a list of spatial objects
        h = []
        #-- for each file in the index
        for i,f in enumerate(file_list):
            if (format == 'ascii'):
                #-- netcdf (.nc)
                h.append(spatial().from_ascii(os.path.expanduser(f),date=date))
            elif (format == 'netCDF4'):
                #-- netcdf (.nc)
                h.append(spatial().from_netCDF4(os.path.expanduser(f),date=date))
            elif (format == 'HDF5'):
                #-- HDF5 (.H5)
                h.append(spatial().from_HDF5(os.path.expanduser(f),date=date))
        #-- create a single spatial object from the list
        return self.from_list(h,date=date,sort=sort)

    def from_list(self, object_list, date=True, sort=True):
        """
        Build a sorted spatial object from a list of other spatial objects
        Inputs: list of spatial object to be merged
        Options:
            spatial objects contain date information
            sort spatial objects by date information
        """
        #-- number of spatial objects in list
        n = len(object_list)
        #-- indices to sort data objects if spatial list contain dates
        if date and sort:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
        else:
            list_sort = np.arange(n)
        #-- truncate to maximum degree and order
        self.spacing = object_list[0].spacing
        self.extents = object_list[0].extents
        self.shape = object_list[0].shape
        #-- create output spatial
        self.data = np.zeros((self.shape[0],self.shape[1],n))
        #-- create list of files
        self.filename = []
        #-- output dates
        if date:
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int)
        #-- for each indice
        for t,i in enumerate(list_sort):
            self.data[:,:,t] = object_list[i].data[:,:].copy()
            if date:
                self.time[t] = object_list[i].time[:].copy()
                self.month[t] = object_list[i].month[:].copy()
            #-- append filename to list
            if getattr(object_list[i], 'filename'):
                self.filename.append(object_list[i].filename)
        #-- update the dimensions
        self.update_dimensions()
        self.update_mask()
        #-- return the single spatial object
        return self

    def from_dict(self, d):
        """
        Convert a dict object to a spatial object
        Inputs: dictionary object to be converted
        """
        #-- find valid keys within dictionary
        k = [k for k in ['lon','lat','data','time','month'] if k in d.keys()]
        #-- assign variables to self
        for key in k:
            setattr(self, key, d[key].copy())
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def to_ascii(self, filename, date=True, verbose=False):
        """
        Write a spatial object to ascii file
        Inputs: full path of output ascii file
        Options: spatial objects contain date information
        """
        self.filename = filename
        print(self.filename) if verbose else None
        #-- open the output file
        fid = open(os.path.expanduser(filename), 'w')
        if date:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f} {3:10.4f}'
        else:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f}'
        #-- write to file for each valid latitude and longitude
        ii,jj = np.nonzero(self.data != self.fill_value)
        for ln,lt,dt in zip(self.lon[jj],self.lat[ii],self.data[ii,jj]):
            print(file_format.format(ln,lt,dt,self.time), file=fid)
        #-- close the output file
        fid.close()

    def to_netCDF4(self, filename, date=True, varname='z',
        units=None, longname=None, title=None, verbose=False):
        """
        Write a spatial object to netCDF4 file
        Inputs: full path of output netCDF4 file
        Options: spatial objects contain date information
        """
        self.filename = filename
        ncdf_write(self.data, self.lon, self.lat, self.time,
            FILENAME=os.path.expanduser(filename), UNITS=units, LONGNAME=longname,
            TIME_UNITS='years', TIME_LONGNAME='Date_in_Decimal_Years',
            TITLE=title, FILL_VALUE=self.fill_value, DATE=date, VERBOSE=verbose)

    def to_HDF5(self, filename, date=True, varname='z',
        units=None, longname=None, title=None, verbose=False):
        """
        Write a spatial object to HDF5 file
        Inputs: full path of output HDF5 file
        Options: spatial objects contain date information
        """
        self.filename = filename
        hdf5_write(self.data, self.lon, self.lat, self.time,
            FILENAME=os.path.expanduser(filename), UNITS=units, LONGNAME=longname,
            TIME_UNITS='years', TIME_LONGNAME='Date_in_Decimal_Years',
            TITLE=title, FILL_VALUE=self.fill_value, DATE=date, VERBOSE=verbose)

    def update_spacing(self):
        """
        Calculate the step size of spatial object
        """
        #-- calculate degree spacing
        dlat = np.abs(self.lat[1] - self.lat[0])
        dlon = np.abs(self.lon[1] - self.lon[0])
        self.spacing = (dlon,dlat)
        return self

    def update_extents(self):
        """
        Calculate the bounds of spatial object
        """
        self.extent[0] = np.min(self.lon)
        self.extent[1] = np.max(self.lon)
        self.extent[2] = np.min(self.lat)
        self.extent[3] = np.max(self.lat)

    def update_dimensions(self):
        """
        Update the dimensions of the spatial object with new extents
        """
        self.shape = np.shape(self.data)
        self.ndim = np.ndim(self.data)
        return self

    def update_mask(self):
        """
        Update the mask of the spatial object
        """
        self.mask = (self.data == self.fill_value)
        return self

    def copy(self):
        """
        Copy a spatial object to a new spatial object
        """
        temp = spatial()
        #-- assign variables to self
        for key in ['lon','lat','data','time','month']:
            val = getattr(self, key)
            setattr(temp, key, np.copy(val))
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a spatial object if non-existent
        """
        #-- change time dimensions to be iterable
        if (np.ndim(self.time) == 0):
            self.time = np.array([self.time])
            self.month = np.array([self.month])
        #-- output spatial with a third dimension
        if (np.ndim(self.data) == 2):
            self.data = self.data[:,:,None]
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        return self

    def squeeze(self):
        """
        Remove singleton dimensions from a spatial object
        """
        #-- squeeze singleton dimensions
        self.time = np.squeeze(self.time)
        self.month = np.squeeze(self.month)
        self.data = np.squeeze(self.data)
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        return self

    def index(self, indice, date=True):
        """
        Subset a spatial object to specific index
        Inputs: indice in matrix to subset
        Options: spatial objects contain date information
        """
        #-- output spatial object
        temp = spatial()
        #-- subset output spatial field
        temp.data = self.data[:,:,indice].copy()
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        #-- subset output dates
        if date:
            temp.time = self.time[indice].copy()
            temp.month = self.month[indice].copy()
        #-- subset filenames
        if getattr(self, 'filename'):
            temp.filename = self.filename[indice]
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        return temp

    def subset(self, months):
        """
        Subset a spatial object to specific GRACE/GRACE-FO months
        Inputs: GRACE/GRACE-FO months
        """
        #-- check if months is an array or a single value
        if (np.ndim(months) == 0):
            months = np.array([months])
        #-- number of months
        n = len(months)
        #-- check that all months are available
        months_check = list(set(months) - set(self.month))
        if months_check:
            m = ','.join(['{0:03d}'.format(m) for m in months_check])
            raise IOError('GRACE/GRACE-FO months {0} not Found'.format(m))
        #-- indices to sort data objects
        months_list = [i for i,m in enumerate(self.month) if m in months]
        #-- output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1])
        #-- create output spatial object
        temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
        temp.lon = np.copy(self.lon)
        temp.lat = np.copy(self.lat)
        temp.time = np.zeros((n))
        temp.month = np.zeros((n),dtype=np.int)
        temp.filename = []
        #-- for each indice
        for t,i in enumerate(months_list):
            temp.data[:,:,t] = self.data[:,:,i].copy()
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
            if getattr(self, 'filename'):
                temp.filename.append(self.filename[i])
        #-- remove singleton dimensions if importing a single value
        return temp.squeeze()

    def mean(self, apply=False):
        """
        Compute mean spatial field and remove from data if specified
        Option: apply to remove the mean field from the input data
        """
        #-- output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1])
        #-- create output mean spatial object
        temp.data = np.mean((self.data),axis=2)
        #-- calculate the spatial anomalies by removing the mean field
        if apply:
            for i,t in enumerate(self.time):
                self.data[:,:,i] -= temp.data[:,:]
        #-- return the mean field
        return temp
