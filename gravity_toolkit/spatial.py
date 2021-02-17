#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (02/2021)

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
    Updated 02/2021: added replace_masked to replace masked values in data
    Updated 01/2021: added scaling factor and scaling factor error function
        from Lander and Swenson (2012) https://doi.org/10.1029/2011WR011453
    Updated 12/2020: added transpose function, can calculate mean over indices
        output attributes dictionary from netCDF4 and HDF5 files
        can create a spatial object from an open file-like object
    Updated 09/2020: added header option to skip rows in ascii files
    Updated 08/2020: added compression options for ascii, netCDF4 and HDF5 files
    Updated 07/2020: added class docstring and using kwargs for output to file
        added case_insensitive_filename function to search directories
    Updated 06/2020: added zeros_like() for creating an empty spatial object
    Written 06/2020
"""
import os
import re
import io
import copy
import gzip
import zipfile
import numpy as np
from gravity_toolkit.ncdf_write import ncdf_write
from gravity_toolkit.hdf5_write import hdf5_write
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.hdf5_read import hdf5_read

class spatial(object):
    """
    Data class for reading, writing and processing spatial data
    """
    np.seterr(invalid='ignore')
    def __init__(self, spacing=[None,None], nlat=None, nlon=None,
        extent=[None]*4, fill_value=None):
        self.data=None
        self.mask=None
        self.lon=None
        self.lat=None
        self.time=None
        self.month=None
        self.fill_value=fill_value
        self.attributes=dict()
        self.extent=extent
        self.spacing=spacing
        self.shape=[nlat,nlon,None]
        self.ndim=None
        self.filename=None

    def case_insensitive_filename(self,filename):
        """
        Searches a directory for a filename without case dependence
        """
        #-- check if filename is open file object
        if isinstance(filename, io.IOBase):
            self.filename = copy.copy(filename)
        else:
            #-- tilde-expand input filename
            self.filename = os.path.expanduser(filename)
            #-- check if file presently exists with input case
            if not os.access(self.filename,os.F_OK):
                #-- search for filename without case dependence
                basename = os.path.basename(filename)
                directory = os.path.dirname(os.path.expanduser(filename))
                f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
                if not f:
                    raise IOError('{0} not found in file system'.format(filename))
                self.filename = os.path.join(directory,f.pop())
        return self

    def from_ascii(self, filename, date=True, compression=None, verbose=False,
        columns=['lon','lat','data','time'], header=0):
        """
        Read a spatial object from an ascii file
        Inputs: full path of input ascii file
        Options:
            ascii file contains date information
            ascii file is compressed or streamed from memory
            verbose output of file information
            column names of ascii file
            rows of header lines to skip
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        print(self.filename) if verbose else None
        #-- open the ascii file and extract contents
        if (compression == 'gzip'):
            #-- read input ascii data from gzip compressed file and split lines
            with gzip.open(self.filename,'r') as f:
                file_contents = f.read().decode('ISO-8859-1').splitlines()
        elif (compression == 'zip'):
            #-- read input ascii data from zipped file and split lines
            base,extension = os.path.splitext(self.filename)
            with zipfile.ZipFile(self.filename) as z:
                file_contents = z.read(base).decode('ISO-8859-1').splitlines()
        elif (compression == 'bytes'):
            #-- read input file object and split lines
            file_contents = self.filename.read().splitlines()
        else:
            #-- read input ascii file (.txt, .asc) and split lines
            with open(self.filename,'r') as f:
                file_contents = f.read().splitlines()
        #-- compile regular expression operator for extracting numerical values
        #-- from input ascii files of spatial data
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        #-- output spatial dimensions
        if (None not in self.extent):
            self.lat = np.linspace(self.extent[3],self.extent[2],self.shape[0])
            self.lon = np.linspace(self.extent[0],self.extent[1],self.shape[1])
        else:
            self.lat = np.zeros((self.shape[0]))
            self.lon = np.zeros((self.shape[1]))
        #-- output spatial data
        self.data = np.zeros((self.shape[0],self.shape[1]))
        self.mask = np.zeros((self.shape[0],self.shape[1]),dtype=np.bool)
        #-- remove time from list of column names if not date
        columns = [c for c in columns if (c != 'time')]
        #-- extract spatial data array and convert to matrix
        #-- for each line in the file
        for line in file_contents[header:]:
            #-- extract columns of interest and assign to dict
            #-- convert fortran exponentials if applicable
            d = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
            #-- convert line coordinates to integers
            ilon = np.int(np.float(d['lon'])/self.spacing[0])
            ilat = np.int((90.0-np.float(d['lat']))//self.spacing[1])
            self.data[ilat,ilon] = np.float(d['data'])
            self.mask[ilat,ilon] = False
            self.lon[ilon] = np.float(d['lon'])
            self.lat[ilat] = np.float(d['lat'])
            #-- if the ascii file contains date variables
            if date:
                self.time = np.array(d['time'],dtype='f')
                self.month = np.array(12.0*(self.time-2002.0)+1,dtype='i')
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def from_netCDF4(self, filename, date=True, compression=None, verbose=False,
        varname='z', lonname='lon', latname='lat', timename='time'):
        """
        Read a spatial object from a netCDF4 file
        Inputs: full path of input netCDF4 file
        Options:
            netCDF4 file contains date information
            netCDF4 file is compressed or streamed from memory
            verbose output of file information
            netCDF4 variable names of data, longitude, latitude and time
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- read data from netCDF4 file
        data = ncdf_read(self.filename, COMPRESSION=compression,
            VERBOSE=verbose, DATE=date, VARNAME=varname,
            LONNAME=lonname, LATNAME=latname, TIMENAME=timename)
        self.data = data['data'].copy()
        if '_FillValue' in data['attributes']['data'].keys():
            self.fill_value = data['attributes']['data']['_FillValue']
        self.mask = np.zeros(self.data.shape, dtype=np.bool)
        self.lon = data['lon'].copy()
        self.lat = data['lat'].copy()
        if date:
            self.time = data['time'].copy()
            self.month = np.array(12.0*(self.time-2002.0)+1,dtype='i')
        #-- update attributes
        self.attributes.update(data['attributes'])
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def from_HDF5(self, filename, date=True, compression=None, verbose=False,
        varname='z', lonname='lon', latname='lat', timename='time'):
        """
        Read a spatial object from a HDF5 file
        Inputs: full path of input HDF5 file
        Options:
            HDF5 file contains date information
            HDF5 file is compressed or streamed from memory
            verbose output of file information
            HDF5 variable names of data, longitude, latitude and time
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- read data from HDF5 file
        data = hdf5_read(self.filename, COMPRESSION=compression,
            VERBOSE=verbose, DATE=date, VARNAME=varname,
            LONNAME=lonname, LATNAME=latname, TIMENAME=timename)
        self.data = data['data'].copy()
        if '_FillValue' in data['attributes']['data'].keys():
            self.fill_value = data['attributes']['_FillValue']
        self.mask = np.zeros(self.data.shape, dtype=np.bool)
        self.lon = data['lon'].copy()
        self.lat = data['lat'].copy()
        if date:
            self.time = data['time'].copy()
            self.month = np.array(12.0*(self.time-2002.0)+1,dtype='i')
        #-- update attributes
        self.attributes.update(data['attributes'])
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
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- Read index file of input spatial data
        with open(self.filename,'r') as f:
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

    def from_list(self, object_list, date=True, sort=True, clear=False):
        """
        Build a sorted spatial object from a list of other spatial objects
        Inputs: list of spatial object to be merged
        Options:
            spatial objects contain date information
            sort spatial objects by date information
            clear the spatial list from memory
        """
        #-- number of spatial objects in list
        n = len(object_list)
        #-- indices to sort data objects if spatial list contain dates
        if date and sort:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
        else:
            list_sort = np.arange(n)
        #-- extract dimensions and grid spacing
        self.spacing = object_list[0].spacing
        self.extent = object_list[0].extent
        self.shape = object_list[0].shape
        #-- create output spatial grid and mask
        self.data = np.zeros((self.shape[0],self.shape[1],n))
        self.mask = np.zeros((self.shape[0],self.shape[1],n),dtype=np.bool)
        self.fill_value = object_list[0].fill_value
        self.lon = object_list[0].lon.copy()
        self.lat = object_list[0].lat.copy()
        #-- create list of files and attributes
        self.filename = []
        self.attributes = []
        #-- output dates
        if date:
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int)
        #-- for each indice
        for t,i in enumerate(list_sort):
            self.data[:,:,t] = object_list[i].data[:,:].copy()
            self.mask[:,:,t] |= object_list[i].mask[:,:]
            if date:
                self.time[t] = np.atleast_1d(object_list[i].time)
                self.month[t] = np.atleast_1d(object_list[i].month)
            #-- append filename to list
            if getattr(object_list[i], 'filename'):
                self.filename.append(object_list[i].filename)
            #-- append attributes to list
            if getattr(object_list[i], 'attributes'):
                self.attributes.append(object_list[i].attributes)
        #-- update the dimensions
        self.update_dimensions()
        self.update_mask()
        #-- clear the input list to free memory
        if clear:
            object_list = None
        #-- return the single spatial object
        return self

    def from_dict(self, d):
        """
        Convert a dict object to a spatial object
        Inputs: dictionary object to be converted
        """
        #-- assign variables to self
        for key in ['lon','lat','data','error','time','month']:
            try:
                setattr(self, key, d[key].copy())
            except (AttributeError, KeyError):
                pass
        #-- create output mask for data
        self.mask = np.zeros_like(self.data,dtype=np.bool)
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
        self.filename = os.path.expanduser(filename)
        print(self.filename) if verbose else None
        #-- open the output file
        fid = open(self.filename, 'w')
        if date:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f} {3:10.4f}'
        else:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f}'
        #-- write to file for each valid latitude and longitude
        ii,jj = np.nonzero((self.data != self.fill_value) & (~self.mask))
        for ln,lt,dt in zip(self.lon[jj],self.lat[ii],self.data[ii,jj]):
            print(file_format.format(ln,lt,dt,self.time), file=fid)
        #-- close the output file
        fid.close()

    def to_netCDF4(self, filename, date=True, **kwargs):
        """
        Write a spatial object to netCDF4 file
        Inputs: full path of output netCDF4 file
        Options: spatial objects contain date information
        **kwargs: keyword arguments for ncdf_write
        """
        self.filename = os.path.expanduser(filename)
        KWARGS = {}
        for key,val in kwargs.items():
            KWARGS[key.upper()] = val
        if 'TIME_UNITS' not in KWARGS.keys():
            KWARGS['TIME_UNITS'] = 'years'
        if 'TIME_LONGNAME' not in KWARGS.keys():
            KWARGS['TIME_LONGNAME'] = 'Date_in_Decimal_Years'
        ncdf_write(self.data, self.lon, self.lat, self.time,
            FILENAME=self.filename, DATE=date,
            FILL_VALUE=self.fill_value, **KWARGS)

    def to_HDF5(self, filename, date=True, **kwargs):
        """
        Write a spatial object to HDF5 file
        Inputs: full path of output HDF5 file
        Options: spatial objects contain date information
        **kwargs: keyword arguments for hdf5_write
        """
        self.filename = os.path.expanduser(filename)
        KWARGS = {}
        for key,val in kwargs.items():
            KWARGS[key.upper()] = val
        if 'TIME_UNITS' not in KWARGS.keys():
            KWARGS['TIME_UNITS'] = 'years'
        if 'TIME_LONGNAME' not in KWARGS.keys():
            KWARGS['TIME_LONGNAME'] = 'Date_in_Decimal_Years'
        hdf5_write(self.data, self.lon, self.lat, self.time,
            FILENAME=self.filename, DATE=date,
            FILL_VALUE=self.fill_value, **KWARGS)

    def to_masked_array(self):
        """
        Convert a spatial object to a masked numpy array
        """
        return np.ma.array(self.data, mask=self.mask,
            fill_value=self.fill_value)

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
        Update the dimensions of the spatial object
        """
        self.shape = np.shape(self.data)
        self.ndim = np.ndim(self.data)
        return self

    def update_mask(self):
        """
        Update the mask of the spatial object
        """
        if self.fill_value is not None:
            self.mask |= (self.data == self.fill_value)
            self.mask |= np.isnan(self.data)
            self.data[self.mask] = self.fill_value
        return self

    def copy(self):
        """
        Copy a spatial object to a new spatial object
        """
        temp = spatial(fill_value=self.fill_value)
        #-- copy attributes or update attributes dictionary
        if isinstance(self.attributes,list):
            setattr(temp,'attributes',self.attributes)
        elif isinstance(self.attributes,dict):
            temp.attributes.update(self.attributes)
        #-- assign variables to self
        var = ['lon','lat','data','mask','error','time','month']
        for key in var:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.copy(val))
            except AttributeError:
                pass
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        temp.replace_masked()
        return temp

    def zeros_like(self):
        """
        Create a spatial object using the dimensions of another
        """
        temp = spatial(fill_value=self.fill_value)
        #-- assign variables to self
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        var = ['data','mask','error','time','month']
        for key in var:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.zeros_like(val))
            except AttributeError:
                pass
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        temp.replace_masked()
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a spatial object if non-existent
        """
        #-- change time dimensions to be iterable
        self.time = np.atleast_1d(self.time)
        self.month = np.atleast_1d(self.month)
        #-- output spatial with a third dimension
        if (np.ndim(self.data) == 2):
            self.data = self.data[:,:,None]
            self.mask = self.mask[:,:,None]
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def squeeze(self):
        """
        Remove singleton dimensions from a spatial object
        """
        #-- squeeze singleton dimensions
        self.time = np.squeeze(self.time)
        self.month = np.squeeze(self.month)
        self.data = np.squeeze(self.data)
        self.mask = np.squeeze(self.mask)
        #-- get spacing and dimensions
        self.update_spacing()
        self.update_extents()
        self.update_dimensions()
        self.update_mask()
        return self

    def index(self, indice, date=True):
        """
        Subset a spatial object to specific index
        Inputs: indice in matrix to subset
        Options: spatial objects contain date information
        """
        #-- output spatial object
        temp = spatial(fill_value=self.fill_value)
        #-- subset output spatial field
        temp.data = self.data[:,:,indice].copy()
        temp.mask = self.mask[:,:,indice].copy()
        #-- subset output spatial error
        try:
            temp.error = self.error[:,:,indice].copy()
        except AttributeError:
            pass
        #-- copy dimensions
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
        months = np.atleast_1d(months)
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
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        #-- create output spatial object
        temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
        temp.mask = np.zeros((temp.shape[0],temp.shape[1],n))
        #-- create output spatial error
        try:
            getattr(self, 'error')
            temp.error = np.zeros((temp.shape[0],temp.shape[1],n))
        except AttributeError:
            pass
        #-- copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        temp.time = np.zeros((n))
        temp.month = np.zeros((n),dtype=np.int)
        temp.filename = []
        #-- for each indice
        for t,i in enumerate(months_list):
            temp.data[:,:,t] = self.data[:,:,i].copy()
            temp.mask[:,:,t] = self.mask[:,:,i].copy()
            try:
                temp.error[:,:,t] = self.error[:,:,i].copy()
            except AttributeError:
                pass
            #-- copy time dimensions
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
            #-- subset filenmaes
            if getattr(self, 'filename'):
                temp.filename.append(self.filename[i])
        #-- remove singleton dimensions if importing a single value
        return temp.squeeze()

    def offset(self, var):
        """
        Offset a spatial object by a constant
        Inputs: scalar value to which the spatial object will be offset
        """
        temp = self.copy()
        #-- offset by a single constant or a time-variable scalar
        if (np.ndim(var) == 0):
            temp.data = self.data + var
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            n = len(var)
            temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
            temp.mask = np.zeros((temp.shape[0],temp.shape[1],n),dtype=np.bool)
            for i,v in enumerate(var):
                temp.data[:,:,i] = self.data[:,:] + v
                temp.mask[:,:,i] = np.copy(self.mask[:,:])
        elif (np.ndim(var) == 1) and (self.ndim == 3):
            for i,v in enumerate(var):
                temp.data[:,:,i] = self.data[:,:,i] + v
        elif (np.ndim(var) == 2) and (self.ndim == 2):
            temp.data = self.data + var
        elif (np.ndim(var) == 2) and (self.ndim == 3):
            for i,t in enumerate(self.time):
                temp.data[:,:,i] = self.data[:,:,i] + var
        elif (np.ndim(var) == 3) and (self.ndim == 3):
            for i,t in enumerate(self.time):
                temp.data[:,:,i] = self.data[:,:,i] + var[:,:,i]
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def scale(self, var):
        """
        Multiply a spatial object by a constant
        Inputs: scalar value to which the spatial object will be multiplied
        """
        temp = self.copy()
        #-- multiply by a single constant or a time-variable scalar
        if (np.ndim(var) == 0):
            temp.data = var*self.data
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            n = len(var)
            temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
            temp.mask = np.zeros((temp.shape[0],temp.shape[1],n),dtype=np.bool)
            for i,v in enumerate(var):
                temp.data[:,:,i] = v*self.data[:,:]
                temp.mask[:,:,i] = np.copy(self.mask[:,:])
        elif (np.ndim(var) == 1) and (self.ndim == 3):
            for i,v in enumerate(var):
                temp.data[:,:,i] = v*self.data[:,:,i]
        elif (np.ndim(var) == 2) and (self.ndim == 2):
            temp.data = var*self.data
        elif (np.ndim(var) == 2) and (self.ndim == 3):
            for i,t in enumerate(self.time):
                temp.data[:,:,i] = var*self.data[:,:,i]
        elif (np.ndim(var) == 3) and (self.ndim == 3):
            for i,t in enumerate(self.time):
                temp.data[:,:,i] = var[:,:,i]*self.data[:,:,i]
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def kfactor(self, var):
        """
        Calculate the scaling factor and scaling factor errors
            from two spatial objects following
            Landerer and Swenson (2012)
        Inputs: spatial object to used for scaling
        Returns: scaling factor and scaling factor error
        """
        #-- copy to not modify original inputs
        temp1 = self.copy()
        temp2 = var.copy()
        #-- expand dimensions and replace invalid values with 0
        temp1.expand_dims().replace_invalid(0.0)
        temp2.expand_dims().replace_invalid(0.0)
        #-- dimensions of input spatial object
        nlat,nlon,nt = temp1.shape
        #-- allocate for scaling factor and scaling factor error
        temp = spatial(nlat=nlat, nlon=nlon, fill_value=0.0)
        temp.data = np.zeros((nlat, nlon))
        temp.error = np.zeros((nlat, nlon))
        #-- copy latitude and longitude variables
        temp.lon = np.copy(temp1.lon)
        temp.lat = np.copy(temp1.lat)
        #-- find valid data points and set mask
        temp.mask = np.any(temp1.mask | temp2.mask, axis=2)
        indy,indx = np.nonzero(np.logical_not(temp.mask))
        #-- calculate point-based scaling factors as centroids
        val1 = np.sum(temp1.data[indy,indx,:]*temp2.data[indy,indx,:],axis=1)
        val2 = np.sum(temp1.data[indy,indx,:]**2,axis=1)
        temp.data[indy,indx] = val1/val2
        #-- calculate difference between scaled and original
        variance = temp1.scale(temp.data).offset(-temp2.data)
        #-- calculate scaling factor errors as RMS of variance
        temp.error = np.sqrt((variance.sum(power=2).data)/nt)
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        #-- return the scaling factors and scaling factor errors
        return temp

    def mean(self, apply=False, indices=Ellipsis):
        """
        Compute mean spatial field and remove from data if specified
        Option:
            apply to remove the mean field from the input data
            indices of spatial object to compute mean
        """
        #-- output spatial object
        temp = spatial(nlon=self.data.shape[0],nlat=self.data.shape[1],
            fill_value=self.fill_value)
        #-- copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        #-- create output mean spatial object
        temp.data = np.mean(self.data[:,:,indices],axis=2)
        temp.mask = np.any(self.mask[:,:,indices],axis=2)
        #-- calculate the mean time
        try:
            val = getattr(self, 'time')
            temp.time = np.mean(val[indices])
        except (AttributeError,TypeError):
            pass
        #-- calculate the spatial anomalies by removing the mean field
        if apply:
            for i,t in enumerate(self.time):
                self.data[:,:,i] -= temp.data[:,:]
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def reverse(self, axis=0):
        """
        Reverse the order of data and dimensions along an axis
        Option: axis to reorder
        """
        #-- output spatial object
        temp = self.copy()
        temp.expand_dims()
        #-- copy dimensions and reverse order
        if (axis == 0):
            temp.lat = temp.lat[::-1].copy()
            temp.data = temp.data[::-1,:,:].copy()
            temp.mask = temp.mask[::-1,:,:].copy()
        elif (axis == 1):
            temp.lon = temp.lon[::-1].copy()
            temp.data = temp.data[:,::-1,:].copy()
            temp.mask = temp.mask[:,::-1,:].copy()
        #-- squeeze output spatial object
        #-- get spacing and dimensions
        #-- update mask
        temp.squeeze()
        return temp

    def transpose(self, axes=None):
        """
        Reverse or permute the axes of a spatial object
        Option: order of the output axes
        """
        #-- output spatial object
        temp = self.copy()
        #-- copy dimensions and reverse order
        temp.data = np.transpose(temp.data, axes=axes)
        temp.mask = np.transpose(temp.mask, axes=axes)
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def sum(self, power=1):
        """
        Compute summation of spatial field
        Option: apply a power before calculating summation
        """
        #-- output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        #-- copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        #-- create output summation spatial object
        temp.data = np.sum(np.power(self.data,power),axis=2)
        temp.mask = np.any(self.mask,axis=2)
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def power(self, power):
        """
        Raise a spatial object to a power
        Inputs: power to which the spatial object will be raised
        """
        temp = self.copy()
        temp.data = np.power(self.data,power)
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        return temp

    def max(self):
        """
        Compute maximum value of spatial field
        """
        #-- output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        #-- copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        #-- create output maximum spatial object
        temp.data = np.max(self.data,axis=2)
        temp.mask = np.any(self.mask,axis=2)
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def min(self):
        """
        Compute minimum value of spatial field
        """
        #-- output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        #-- copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        #-- create output minimum spatial object
        temp.data = np.min(self.data,axis=2)
        temp.mask = np.any(self.mask,axis=2)
        #-- get spacing and dimensions
        temp.update_spacing()
        temp.update_extents()
        temp.update_dimensions()
        #-- update mask
        temp.update_mask()
        return temp

    def replace_invalid(self, fill_value, mask=None):
        """
        Replace the masked values with a new fill_value
        """
        #-- validate current mask
        self.update_mask()
        #-- update the mask if specified
        if mask is not None:
            if (np.shape(mask) == self.shape):
                self.mask |= mask
            elif (np.ndim(mask) == 2) & (self.ndim == 3):
                temp = np.broadcast_to(mask, self.shape)
                self.mask |= temp
        #-- update the fill value
        self.fill_value = fill_value
        #-- replace invalid values with new fill value
        self.data[self.mask] = self.fill_value
        return self

    def replace_masked(self):
        """
        Replace the masked values with fill_value
        """
        if self.fill_value is not None:
            self.data[self.mask] = self.fill_value
        return self