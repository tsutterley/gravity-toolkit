#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (10/2024)

Data class for reading, writing and processing spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 10/2024: allow 2D and 3D arrays in output netCDF4 files
    Updated 06/2024: use wrapper to importlib for optional dependencies
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 10/2023: place time and month variables in try/except block
    Updated 09/2023: prevent double printing of filenames when using debug
    Updated 08/2023: add string representation of the spatial object
    Updated 05/2023: use pathlib to define and operate on paths
        more operatations on spatial error if in possible data keys
        rename reverse function to flip to match numpy nomenclature
    Updated 03/2023: customizable file-level attributes to netCDF4 and HDF5
        add attributes fetching to from_dict function
        retrieve all root attributes from HDF5 and netCDF4 datasets
        fix indexing of filenames in single string case
        add indexing of filenames to spatial object iterator
        use copy.copy and not numpy.copy in copy spatial object function
        fix mask and shape of subsetted spatial grid objects
        add extend_matrix function and add error output to from_list
        convert spacing, extent, shape and ndim to spatial class properties
        improve typing for variables in docstrings
        set case insensitive filename to None if filename is empty
    Updated 02/2023: use monospaced text to note spatial objects in docstrings
    Updated 12/2022: add software information to output HDF5 and netCDF4
        make spatial objects iterable and with length
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: fix output latitude HDF5 and netCDF4 attributes
        place index filename within try/except statement
    Updated 04/2022: updated docstrings to numpy documentation format
        using internal netCDF4 and HDF5 readers and writers
        include utf-8 encoding in reads to be windows compliant
        include filename attribute when copying spatial objects
    Updated 12/2021: logging case_insensitive_filename output for debugging
    Updated 11/2021: fix kwargs to index and hdf5 read functions
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: use functions for converting to and from GRACE months
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
        use file not found exceptions in case insensitive filename
    Updated 02/2021: added replace_masked to replace masked values in data
        use adjust_months function to fix special cases of GRACE/GRACE-FO months
        added generic reader, generic writer and write to list functions
        generalize ascii, netCDF4 and HDF5 readers and writers
        replaced numpy bool to prevent deprecation warning
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
import re
import io
import copy
import gzip
import time
import uuid
import logging
import pathlib
import zipfile
import numpy as np
import gravity_toolkit.version
from gravity_toolkit.time import adjust_months, calendar_to_grace
from gravity_toolkit.utilities import import_dependency

# attempt imports
h5py = import_dependency('h5py')
netCDF4 = import_dependency('netCDF4')

class spatial(object):
    """
    Data class for reading, writing and processing spatial data

    Attributes
    ----------
    data: np.ndarray
        spatial grid data
    mask: np.ndarray
        spatial grid mask
    lon: np.ndarray
        grid longitudes
    lat: np.ndarray
        grid latitudes
    time: np.ndarray
        time variable of the spatial data
    month: np.ndarray
        GRACE/GRACE-FO months variable of the spatial data
    fill_value: float or NoneType, default None
        invalid value for spatial grid data
    attributes: dict
        attributes of ``spatial`` variables
    filename: str
        input or output filename

    """
    np.seterr(invalid='ignore')
    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault('fill_value',None)
        # set default class attributes
        self.data=None
        self.mask=None
        self.lon=None
        self.lat=None
        self.time=None
        self.month=None
        self.fill_value=kwargs['fill_value']
        self.attributes=dict()
        self.filename=None
        # iterator
        self.__index__ = 0

    def case_insensitive_filename(self, filename):
        """
        Searches a directory for a filename without case dependence

        Parameters
        ----------
        filename: str, io.IOBase, pathlib.Path or None
            input filename
        """
        # check if filename is open file object
        if isinstance(filename, io.IOBase):
            self.filename = copy.copy(filename)
        elif isinstance(filename, type(None)) or not bool(filename):
            self.filename = None
        else:
            # tilde-expand input filename
            self.filename = pathlib.Path(filename).expanduser().absolute()
            # check if file presently exists with input case
            if not self.filename.exists():
                # search for filename without case dependence
                f = [f.name for f in self.filename.parent.iterdir() if
                    re.match(self.filename.name, f.name, re.I)]
                if not f:
                    msg = f'{filename} not found in file system'
                    raise FileNotFoundError(msg)
                self.filename = self.filename.with_name(f.pop())
        # return the filename
        return self

    def compressuser(self, filename=None):
        """
        Tilde-compresses a file to be relative to the home directory

        Parameters
        ----------
        filename: str or None, default None
            output filename
        """
        if filename is None:
            filename = self.filename
        else:
            filename = pathlib.Path(filename).expanduser().absolute()
        # attempt to compress filename relative to home directory
        try:
            relative_to = filename.relative_to(pathlib.Path().home())
        except (ValueError, AttributeError) as exc:
            return filename
        else:
            return pathlib.Path('~').joinpath(relative_to)

    def from_ascii(self, filename, date=True, **kwargs):
        """
        Read a ``spatial`` object from an ascii file

        Parameters
        ----------
        filename: str
            full path of input ascii file
        date: bool, default True
            ascii file has date information
        compression: str or NoneType, default None
            file compression type

                - ``'gzip'``
                - ``'zip'``
                - ``'bytes'``
        spacing: list, default [None,None]
            grid step size ``[longitude,latitude]``
        extent: list, default [None,None,None,None]
            spatial grid bounds
            ``[minimum longitude, maximum longitude,
            minimum latitude, maximum latitude]``
        nlat: int or NoneType, default None
            length of latitude dimension
        nlon: int or NoneType, default None
            length of longitude dimension
        columns: list, default ['lon','lat','data','time']
            variable names for each column
        header: int, default 0
            Number of rows of header lines to skip
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('verbose',False)
        kwargs.setdefault('compression',None)
        kwargs.setdefault('spacing',[None,None])
        kwargs.setdefault('nlat',None)
        kwargs.setdefault('nlon',None)
        kwargs.setdefault('extent',[None]*4)
        kwargs.setdefault('columns',['lon','lat','data','time'])
        kwargs.setdefault('header',0)
        # open the ascii file and extract contents
        logging.info(str(self.filename))
        if (kwargs['compression'] == 'gzip'):
            # read input ascii data from gzip compressed file and split lines
            with gzip.open(self.filename, mode='r') as f:
                file_contents = f.read().decode('ISO-8859-1').splitlines()
        elif (kwargs['compression'] == 'zip'):
            # read input ascii data from zipped file and split lines
            stem = self.filename.stem
            with zipfile.ZipFile(self.filename) as z:
                file_contents = z.read(stem).decode('ISO-8859-1').splitlines()
        elif (kwargs['compression'] == 'bytes'):
            # read input file object and split lines
            file_contents = self.filename.read().splitlines()
        else:
            # read input ascii file (.txt, .asc) and split lines
            with self.filename.open(mode='r', encoding='utf8') as f:
                file_contents = f.read().splitlines()
        # compile regular expression operator for extracting numerical values
        # from input ascii files of spatial data
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        # output spatial dimensions
        if (None not in kwargs['extent']) and kwargs['nlat'] and kwargs['nlon']:
            extent = kwargs.get('extent')
            self.lat = np.linspace(extent[3], extent[2], kwargs['nlat'])
            self.lon = np.linspace(extent[0], extent[1], kwargs['nlon'])
            dlon = np.abs(self.lon[1] - self.lon[0])
            dlat = np.abs(self.lat[1] - self.lat[0])
        elif (None not in kwargs['extent']) and (None not in kwargs['spacing']):
            extent = kwargs.get('extent')
            dlon, dlat = kwargs.get('spacing')
            self.lat = np.arange(extent[3], extent[2] - dlat, dlat)
            self.lon = np.arange(extent[0], extent[1] + dlon, dlon)
        elif kwargs['nlat'] and kwargs['nlon'] and (None not in kwargs['spacing']):
            dlon, dlat = kwargs.get('spacing')
            self.lat = np.zeros((kwargs['nlat']))
            self.lon = np.zeros((kwargs['nlon']))
        else:
            raise ValueError('Unknown dimensions for input ``spatial`` object')
        # get spatial dimensions
        nlat = len(self.lat)
        nlon = len(self.lon)
        # output spatial data
        self.data = np.zeros((nlat, nlon))
        self.mask = np.zeros((nlat, nlon), dtype=bool)
        # remove time from list of column names if not date
        columns = [c for c in kwargs['columns'] if (c != 'time')]
        # extract spatial data array and convert to matrix
        # for each line in the file
        header = kwargs['header']
        for line in file_contents[header:]:
            # extract columns of interest and assign to dict
            # convert fortran exponentials if applicable
            d = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
            # convert line coordinates to integers
            ilon = np.int64(np.float64(d['lon'])/dlon)
            ilat = np.int64((90.0 - np.float64(d['lat']))//dlat)
            self.data[ilat, ilon] = np.float64(d['data'])
            self.mask[ilat, ilon] = False
            self.lon[ilon] = np.float64(d['lon'])
            self.lat[ilat] = np.float64(d['lat'])
            # if the ascii file contains date variables
            if date:
                self.time = np.array(d['time'],dtype='f')
                self.month = calendar_to_grace(self.time)
        # if the ascii file contains date variables
        if date:
            # adjust months to fix special cases if necessary
            self.month = adjust_months(self.month)
        # update mask
        self.update_mask()
        return self

    def from_netCDF4(self, filename, **kwargs):
        """
        Read a ``spatial`` object from a netCDF4 file

        Parameters
        ----------
        filename: str
            full path of input netCDF4 file
        date: bool, default True
            netCDF4 file has date information
        compression: str or NoneType, default None
            file compression type

                - ``'gzip'``
                - ``'zip'``
                - ``'bytes'``
        varname: str, default 'data'
            name for data variable
        lonname: str, default 'lon'
            name for longitude variable
        latname: str, default 'lat'
            name for latitude variable
        timename: str, default 'time'
            name for time-dimension variable
        field_mapping: dict, default {}
            mapping between output variables and input netCDF4
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date',True)
        kwargs.setdefault('compression',None)
        kwargs.setdefault('varname','z')
        kwargs.setdefault('lonname','lon')
        kwargs.setdefault('latname','lat')
        kwargs.setdefault('timename','time')
        kwargs.setdefault('field_mapping',{})
        kwargs.setdefault('verbose',False)
        # Open the NetCDF4 file for reading
        if (kwargs['compression'] == 'gzip'):
            # read as in-memory (diskless) netCDF4 dataset
            with gzip.open(self.filename, mode='r') as f:
                fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
        elif (kwargs['compression'] == 'zip'):
            # read zipped file and extract file into in-memory file object
            stem = self.filename.stem
            with zipfile.ZipFile(self.filename) as z:
                # first try finding a netCDF4 file with same base filename
                # if none found simply try searching for a netCDF4 file
                try:
                    f,=[f for f in z.namelist() if re.match(stem,f,re.I)]
                except:
                    f,=[f for f in z.namelist() if re.search(r'\.nc(4)?$',f)]
                # read bytes from zipfile as in-memory (diskless) netCDF4 dataset
                fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=z.read(f))
        elif (kwargs['compression'] == 'bytes'):
            # read as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=filename.read())
        else:
            # read netCDF4 dataset
            fileID = netCDF4.Dataset(self.filename, 'r')
        # Output NetCDF file information
        logging.info(fileID.filepath())
        logging.info(list(fileID.variables.keys()))
        # set automasking
        fileID.set_auto_mask(False)
        # list of variable attributes
        attributes_list = ['description','units','long_name','calendar',
            'standard_name','_FillValue','missing_value']
        # mapping between output keys and netCDF4 variable names
        if not kwargs['field_mapping']:
            fields = [kwargs['lonname'],kwargs['latname'],kwargs['varname']]
            if kwargs['date']:
                fields.append(kwargs['timename'])
            kwargs['field_mapping'] = self.default_field_mapping(fields)
        # for each variable
        for field,key in kwargs['field_mapping'].items():
            # Getting the data from each NetCDF variable
            # remove singleton dimensions
            setattr(self, field, np.squeeze(fileID.variables[key][:]))
            # Getting attributes of included variables
            self.attributes[field] = {}
            for attr in attributes_list:
                # try getting the attribute
                try:
                    self.attributes[field][attr] = \
                        fileID.variables[key].getncattr(attr)
                except (KeyError,ValueError,AttributeError):
                    pass
        # get global netCDF4 attributes
        self.attributes['ROOT'] = {}
        for att_name in fileID.ncattrs():
            self.attributes['ROOT'][att_name] = fileID.getncattr(att_name)
        # Closing the NetCDF file
        fileID.close()
        # switching data array to lat/lon if lon/lat
        sz = self.data.shape
        if (self.data.ndim == 2) and (len(self.lon) == sz[0]):
            self.data = self.data.T
        # set fill value and mask
        if '_FillValue' in self.attributes['data'].keys():
            self.fill_value = self.attributes['data']['_FillValue']
            self.mask = (self.data == self.fill_value)
        else:
            self.mask = np.zeros(self.data.shape, dtype=bool)
        # set GRACE/GRACE-FO month if file has date variables
        if kwargs['date']:
            self.month = calendar_to_grace(self.time)
            # adjust months to fix special cases if necessary
            self.month = adjust_months(self.month)
        # update mask
        self.update_mask()
        return self

    def from_HDF5(self, filename, **kwargs):
        """
        Read a ``spatial`` object from a HDF5 file

        Parameters
        ----------
        filename: str
            full path of input HDF5 file
        date: bool, default True
            HDF5 file has date information
        compression: str or NoneType, default None
            file compression type

                - ``'gzip'``
                - ``'zip'``
                - ``'bytes'``
        varname: str, default 'data'
            name for data variable
        lonname: str, default 'lon'
            name for longitude variable
        latname: str, default 'lat'
            name for latitude variable
        timename: str, default 'time'
            name for time-dimension variable
        field_mapping: dict, default {}
            mapping between output variables and input HDF5
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date',True)
        kwargs.setdefault('compression',None)
        kwargs.setdefault('varname','z')
        kwargs.setdefault('lonname','lon')
        kwargs.setdefault('latname','lat')
        kwargs.setdefault('timename','time')
        kwargs.setdefault('field_mapping',{})
        kwargs.setdefault('verbose',False)
        # Open the HDF5 file for reading
        if (kwargs['compression'] == 'gzip'):
            # read gzip compressed file and extract into in-memory file object
            with gzip.open(self.filename, mode='r') as f:
                fid = io.BytesIO(f.read())
            # set filename of BytesIO object
            fid.filename = self.filename.name
            # rewind to start of file
            fid.seek(0)
            # read as in-memory (diskless) HDF5 dataset from BytesIO object
            fileID = h5py.File(fid, 'r')
        elif (kwargs['compression'] == 'zip'):
            # read zipped file and extract file into in-memory file object
            stem = self.filename.stem
            with zipfile.ZipFile(self.filename) as z:
                # first try finding a HDF5 file with same base filename
                # if none found simply try searching for a HDF5 file
                try:
                    f,=[f for f in z.namelist() if re.match(stem,f,re.I)]
                except:
                    f,=[f for f in z.namelist() if re.search(r'\.H(DF)?5$',f,re.I)]
                # read bytes from zipfile into in-memory BytesIO object
                fid = io.BytesIO(z.read(f))
            # set filename of BytesIO object
            fid.filename = self.filename.name
            # rewind to start of file
            fid.seek(0)
            # read as in-memory (diskless) HDF5 dataset from BytesIO object
            fileID = h5py.File(fid, mode='r')
        elif (kwargs['compression'] == 'bytes'):
            # read as in-memory (diskless) HDF5 dataset
            fileID = h5py.File(filename, mode='r')
        else:
            # read HDF5 dataset
            fileID = h5py.File(self.filename, 'r')
        # Output HDF5 file information
        logging.info(fileID.filename)
        logging.info(list(fileID.keys()))
        # list of variable attributes
        attributes_list = ['description','units','long_name','calendar',
            'standard_name','_FillValue','missing_value']
        # mapping between output keys and HDF5 variable names
        if not kwargs['field_mapping']:
            fields = [kwargs['lonname'],kwargs['latname'],kwargs['varname']]
            if kwargs['date']:
                fields.append(kwargs['timename'])
            kwargs['field_mapping'] = self.default_field_mapping(fields)
        # for each variable
        for field,key in kwargs['field_mapping'].items():
            # Getting the data from each HDF5 variable
            # remove singleton dimensions
            setattr(self, field, np.squeeze(fileID[key][:]))
            # Getting attributes of included variables
            self.attributes[field] = {}
            for attr in attributes_list:
                try:
                    self.attributes[field][attr] = fileID[key].attrs[attr]
                except (KeyError, AttributeError):
                    pass
        # get global HDF5 attributes
        self.attributes['ROOT'] = {}
        for att_name,att_val in fileID.attrs.items():
            self.attributes['ROOT'][att_name] = att_val
        # Closing the HDF5 file
        fileID.close()
        # switching data array to lat/lon if lon/lat
        sz = self.data.shape
        if (self.data.ndim == 2) and (len(self.lon) == sz[0]):
            self.data = self.data.T
        # set fill value and mask
        if '_FillValue' in self.attributes['data'].keys():
            self.fill_value = self.attributes['data']['_FillValue']
            self.mask = (self.data == self.fill_value)
        else:
            self.mask = np.zeros(self.data.shape, dtype=bool)
        # set GRACE/GRACE-FO month if file has date variables
        if kwargs['date']:
            self.month = calendar_to_grace(self.time)
            # adjust months to fix special cases if necessary
            self.month = adjust_months(self.month)
        # update mask
        self.update_mask()
        return self

    def from_index(self, filename, **kwargs):
        """
        Read a ``spatial`` object from an index of
        ascii, netCDF4 or HDF5 files

        Parameters
        ----------
        filename: str
            full path of index file
        format: str or NoneType, default None
            format of individual files within index

                - ``'ascii'``
                - ``'netCDF4'``
                - ``'HDF5'``
        date: bool, default True
            files contains date information
        sort: bool, default True
            sort ``spatial`` objects by date information
        **kwargs: dict
            keyword arguments for input readers
        """
        # set default keyword arguments
        kwargs.setdefault('format',None)
        kwargs.setdefault('date',True)
        kwargs.setdefault('sort',True)
        # set filename
        self.case_insensitive_filename(filename)
        # file parser for reading index files
        # removes commented lines (can comment out files in the index)
        # removes empty lines (if there are extra empty lines)
        parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)
        # Read index file of input spatial data
        with self.filename.open(mode='r', encoding='utf8') as f:
            file_list = [l for l in f.read().splitlines() if parser.match(l)]
        # create a list of spatial objects
        s = []
        # for each file in the index
        for i,f in enumerate(file_list):
            if (kwargs['format'] == 'ascii'):
                # netcdf (.nc)
                s.append(spatial().from_ascii(f, **kwargs))
            elif (kwargs['format'] == 'netCDF4'):
                # netcdf (.nc)
                s.append(spatial().from_netCDF4(f, **kwargs))
            elif (kwargs['format'] == 'HDF5'):
                # HDF5 (.H5)
                s.append(spatial().from_HDF5(f, **kwargs))
        # create a single spatial object from the list
        return self.from_list(s,date=kwargs['date'],sort=kwargs['sort'])

    def from_list(self, object_list, **kwargs):
        """
        Build a sorted ``spatial`` object from a list of
        other ``spatial`` objects

        Parameters
        ----------
        object_list: list
            list of ``spatial`` objects to be merged
        date: bool, default True
            files contains date information
        sort: bool, default True
            sort ``spatial`` objects by date information
        clear: bool, default True
            clear the list of ``spatial`` objects from memory
        """
        # set default keyword arguments
        kwargs.setdefault('date',True)
        kwargs.setdefault('sort',True)
        kwargs.setdefault('clear',False)
        # number of spatial objects in list
        n = len(object_list)
        # indices to sort data objects if spatial list contain dates
        if kwargs['date'] and kwargs['sort']:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
        else:
            list_sort = np.arange(n)
        # extract grid spacing
        shape = object_list[0].shape
        # create output spatial grid and mask
        self.data = np.zeros((shape[0], shape[1], n))
        self.mask = np.zeros((shape[0], shape[1], n),dtype=bool)
        # add error if in original list attributes
        if hasattr(object_list[0], 'error'):
            self.error = np.zeros((shape[0], shape[1], n))
        self.fill_value = object_list[0].fill_value
        self.lon = object_list[0].lon.copy()
        self.lat = object_list[0].lat.copy()
        # create list of files and attributes
        self.filename = []
        self.attributes = []
        # output dates
        if kwargs['date']:
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int64)
        # for each indice
        for t,i in enumerate(list_sort):
            self.data[:,:,t] = object_list[i].data[:,:].copy()
            self.mask[:,:,t] |= object_list[i].mask[:,:]
            if hasattr(object_list[i], 'error'):
                self.error[:,:,t] = object_list[i].error[:,:].copy()
            if kwargs['date']:
                self.time[t] = np.atleast_1d(object_list[i].time)
                self.month[t] = np.atleast_1d(object_list[i].month)
            # append filename to list
            if getattr(object_list[i], 'filename'):
                self.filename.append(object_list[i].filename)
            # append attributes to list
            if getattr(object_list[i], 'attributes'):
                self.attributes.append(object_list[i].attributes)
        # adjust months to fix special cases if necessary
        if kwargs['date']:
            self.month = adjust_months(self.month)
        # update mask
        self.update_mask()
        # clear the input list to free memory
        if kwargs['clear']:
            object_list = None
        # return the single spatial object
        return self

    def from_file(self, filename, format=None, date=True, **kwargs):
        """
        Read a ``spatial`` object from a specified format

        Parameters
        ----------
        filename: str
            full path of input file
        format: str or NoneType, default None
            file format

                - ``'ascii'``
                - ``'netCDF4'``
                - ``'HDF5'``
        date: bool, default True
            file contains date information
        verbose: bool, default False
            print file and variable information
        **kwargs: dict
            keyword arguments for input readers
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default verbosity
        kwargs.setdefault('verbose',False)
        # read from file
        if (format == 'ascii'):
            # ascii (.txt)
            return spatial().from_ascii(filename, date=date, **kwargs)
        elif (format == 'netCDF4'):
            # netcdf (.nc)
            return spatial().from_netCDF4(filename, date=date, **kwargs)
        elif (format == 'HDF5'):
            # HDF5 (.H5)
            return spatial().from_HDF5(filename, date=date, **kwargs)

    def from_dict(self, d, **kwargs):
        """
        Convert a ``dict`` object to a ``spatial`` object

        Parameters
        ----------
        d: dict
            dictionary object to be converted
        """
        # assign variables to self
        for key in ['lon','lat','data','error','time','month','directory']:
            try:
                setattr(self, key, d[key].copy())
            except (AttributeError, KeyError):
                pass
        # create output mask for data
        self.mask = np.zeros_like(self.data, dtype=bool)
        # add attributes to root if in dictionary
        self.attributes['ROOT'] = d.get('attributes')
        # update mask
        self.update_mask()
        return self

    def to_ascii(self, filename, **kwargs):
        """
        Write a ``spatial`` object to ascii file

        Parameters
        ----------
        filename: str
            full path of output ascii file
        date: bool, default True
            ``spatial`` objects contain date information
        verbose: bool, default False
            Output file and variable information
        """
        self.filename = pathlib.Path(filename).expanduser().absolute()
        # set default verbosity and parameters
        kwargs.setdefault('date',True)
        kwargs.setdefault('verbose',False)
        logging.info(str(self.filename))
        # open the output file
        fid = self.filename.open(mode='w', encoding='utf8')
        if hasattr(self, 'error') and kwargs['date']:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f} {3:12.4f} {4:10.4f}'
        elif hasattr(self, 'error'):
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f} {3:12.4f}'
        elif kwargs['date']:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f} {4:10.4f}'
        else:
            file_format = '{0:10.4f} {1:10.4f} {2:12.4f}'
        # write to file for each valid latitude and longitude
        ii,jj = np.nonzero((self.data != self.fill_value) & (~self.mask))
        for i,j in zip(ii,jj):
            ln = self.lon[j]
            lt = self.lat[i]
            data = self.data[i,j]
            error = self.error[i,j] if hasattr(self, 'error') else 0.0
            print(file_format.format(ln,lt,data,error,self.time), file=fid)
        # close the output file
        fid.close()

    def to_netCDF4(self, filename, **kwargs):
        """
        Write a ``spatial`` object to netCDF4 file

        Parameters
        ----------
        filename: str
            full path of output netCDF4 file
        varname: str, default 'z'
            data variable name in netCDF4 file
        lonname: str
            longitude variable name in netCDF4 file
        latname: str
            latitude variable name in netCDF4 file
        field_mapping: dict, default {}
            mapping between input variables and output netCDF4
        attributes: dict, default {}
            output netCDF4 variable and file-level attributes
        units: str or NoneType, default: None
            data variable units
        longname: str or NoneType, default: None
            data variable unit description
        time_units: str, default 'years'
            time variable units
        time_longname: str, default 'Date_in_Decimal_Years'
            time variable unit description
        title: str or NoneType, default None
            title attribute of dataset
        source: str or NoneType, default None
            source attribute of dataset
        reference: str or NoneType, default None
            reference attribute of dataset
        date: bool, default True
            ``spatial`` objects contain date information
        clobber: bool, default True
            Overwrite an existing netCDF4 file
        verbose: bool, default False
            Output file and variable information
        """
        # set default verbosity and parameters
        kwargs.setdefault('verbose',False)
        kwargs.setdefault('varname','z')
        kwargs.setdefault('lonname','lon')
        kwargs.setdefault('latname','lat')
        kwargs.setdefault('timename','time')
        kwargs.setdefault('field_mapping',{})
        attributes = self.attributes.get('ROOT') or {}
        kwargs.setdefault('attributes',dict(ROOT=attributes))
        kwargs.setdefault('units',None)
        kwargs.setdefault('longname',None)
        kwargs.setdefault('time_units','years')
        kwargs.setdefault('time_longname','Date_in_Decimal_Years')
        kwargs.setdefault('title',None)
        kwargs.setdefault('source',None)
        kwargs.setdefault('reference',None)
        kwargs.setdefault('date',True)
        kwargs.setdefault('clobber',True)
        kwargs.setdefault('verbose',False)
        # setting NetCDF clobber attribute
        clobber = 'w' if kwargs['clobber'] else 'a'
        # opening NetCDF file for writing
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fileID = netCDF4.Dataset(self.filename, clobber, format="NETCDF4")
        # mapping between output keys and netCDF4 variable names
        if not kwargs['field_mapping']:
            fields = [kwargs['lonname'],kwargs['latname'],kwargs['varname']]
            if kwargs['date']:
                fields.append(kwargs['timename'])
            kwargs['field_mapping'] = self.default_field_mapping(fields)
        # create attributes dictionary for output variables
        if not all(key in kwargs['attributes'] for key in kwargs['field_mapping'].values()):
            # Defining attributes for longitude and latitude
            kwargs['attributes'][kwargs['field_mapping']['lon']] = {}
            kwargs['attributes'][kwargs['field_mapping']['lon']]['long_name'] = 'longitude'
            kwargs['attributes'][kwargs['field_mapping']['lon']]['units'] = 'degrees_east'
            kwargs['attributes'][kwargs['field_mapping']['lat']] = {}
            kwargs['attributes'][kwargs['field_mapping']['lat']]['long_name'] = 'latitude'
            kwargs['attributes'][kwargs['field_mapping']['lat']]['units'] = 'degrees_north'
            # Defining attributes for dataset
            kwargs['attributes'][kwargs['field_mapping']['data']] = {}
            kwargs['attributes'][kwargs['field_mapping']['data']]['long_name'] = kwargs['longname']
            kwargs['attributes'][kwargs['field_mapping']['data']]['units'] = kwargs['units']
            # Defining attributes for date if applicable
            if kwargs['date']:
                kwargs['attributes'][kwargs['field_mapping']['time']] = {}
                kwargs['attributes'][kwargs['field_mapping']['time']]['long_name'] = kwargs['time_longname']
                kwargs['attributes'][kwargs['field_mapping']['time']]['units'] = kwargs['time_units']
        # add default global (file-level) attributes
        if kwargs['title']:
            kwargs['attributes']['ROOT']['title'] = kwargs['title']
        if kwargs['source']:
            kwargs['attributes']['ROOT']['source'] = kwargs['source']
        if kwargs['reference']:
            kwargs['attributes']['ROOT']['reference'] = kwargs['reference']
        # netCDF4 dimension variables
        dimensions = []
        dimensions.append('lat')
        dimensions.append('lon')
        # expand dimensions if containing date variables
        if kwargs['date']:
            self.expand_dims()
            dimensions.append('time')
        dims = tuple(kwargs['field_mapping'][key] for key in dimensions)
        # defining the NetCDF dimensions and variables
        nc = {}
        # NetCDF dimensions
        for i,field in enumerate(dimensions):
            temp = getattr(self,field)
            key = kwargs['field_mapping'][field]
            fileID.createDimension(key, len(temp))
            nc[key] = fileID.createVariable(key, temp.dtype, (key,))
        # NetCDF spatial data
        variables = set(kwargs['field_mapping'].keys()) - set(dimensions)
        for field in sorted(variables):
            temp = getattr(self,field)
            ndim = temp.ndim
            key = kwargs['field_mapping'][field]
            nc[key] = fileID.createVariable(key, temp.dtype, dims[:ndim],
                fill_value=self.fill_value, zlib=True)
        # filling NetCDF variables
        for field,key in kwargs['field_mapping'].items():
            nc[key][:] = getattr(self,field)
            # filling netCDF dataset attributes
            for att_name,att_val in kwargs['attributes'][key].items():
                # skip variable attribute if None
                if not att_val:
                    continue
                # skip variable attributes if in list
                if att_name not in ('DIMENSION_LIST','CLASS','NAME','_FillValue'):
                    nc[key].setncattr(att_name, att_val)
        # global attributes of NetCDF4 file
        for att_name,att_val in kwargs['attributes']['ROOT'].items():
            fileID.setncattr(att_name, att_val)
        # add software information
        fileID.software_reference = gravity_toolkit.version.project_name
        fileID.software_version = gravity_toolkit.version.full_version
        # date created
        fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
        # Output NetCDF structure information
        logging.info(str(self.filename))
        logging.info(list(fileID.variables.keys()))
        # Closing the NetCDF file
        fileID.close()

    def to_HDF5(self, filename, **kwargs):
        """
        Write a ``spatial`` object to HDF5 file

        Parameters
        ----------
        filename: str
            full path of output HDF5 file
        varname: str, default 'z'
            data variable name in HDF5 file
        lonname: str
            longitude variable name in HDF5 file
        latname: str
            latitude variable name in HDF5 file
        field_mapping: dict, default {}
            mapping between input variables and output HDF5
        attributes: dict, default {}
            output HDF5 variable and file-level attributes
        units: str or NoneType, default: None
            data variable units
        longname: str or NoneType, default: None
            data variable unit description
        time_units: str, default 'years'
            time variable units
        time_longname: str, default 'Date_in_Decimal_Years'
            time variable unit description
        title: str or NoneType, default None
            description attribute of dataset
        source: str or NoneType, default None
            source attribute of dataset
        reference: str or NoneType, default None
            reference attribute of dataset
        date: bool, default True
            ``spatial`` objects contain date information
        clobber: bool, default True
            Overwrite an existing HDF5 file
        verbose: bool, default False
            Output file and variable information
        """
        # set default verbosity and parameters
        kwargs.setdefault('verbose',False)
        kwargs.setdefault('varname','z')
        kwargs.setdefault('lonname','lon')
        kwargs.setdefault('latname','lat')
        kwargs.setdefault('timename','time')
        kwargs.setdefault('field_mapping',{})
        attributes = self.attributes.get('ROOT') or {}
        kwargs.setdefault('attributes',dict(ROOT=attributes))
        kwargs.setdefault('units',None)
        kwargs.setdefault('longname',None)
        kwargs.setdefault('time_units','years')
        kwargs.setdefault('time_longname','Date_in_Decimal_Years')
        kwargs.setdefault('title',None)
        kwargs.setdefault('source',None)
        kwargs.setdefault('reference',None)
        kwargs.setdefault('date',True)
        kwargs.setdefault('clobber',True)
        kwargs.setdefault('verbose',False)
        # setting NetCDF clobber attribute
        clobber = 'w' if kwargs['clobber'] else 'w-'
        # opening NetCDF file for writing
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fileID = h5py.File(self.filename, clobber)
        # mapping between output keys and HDF5 variable names
        if not kwargs['field_mapping']:
            fields = [kwargs['lonname'],kwargs['latname'],kwargs['varname']]
            if kwargs['date']:
                fields.append(kwargs['timename'])
            kwargs['field_mapping'] = self.default_field_mapping(fields)
        # create attributes dictionary for output variables
        if not all(key in kwargs['attributes'] for key in kwargs['field_mapping'].values()):
            # Defining attributes for longitude and latitude
            kwargs['attributes'][kwargs['field_mapping']['lon']] = {}
            kwargs['attributes'][kwargs['field_mapping']['lon']]['long_name'] = 'longitude'
            kwargs['attributes'][kwargs['field_mapping']['lon']]['units'] = 'degrees_east'
            kwargs['attributes'][kwargs['field_mapping']['lat']] = {}
            kwargs['attributes'][kwargs['field_mapping']['lat']]['long_name'] = 'latitude'
            kwargs['attributes'][kwargs['field_mapping']['lat']]['units'] = 'degrees_north'
            # Defining attributes for dataset
            kwargs['attributes'][kwargs['field_mapping']['data']] = {}
            kwargs['attributes'][kwargs['field_mapping']['data']]['long_name'] = kwargs['longname']
            kwargs['attributes'][kwargs['field_mapping']['data']]['units'] = kwargs['units']
            # Defining attributes for date if applicable
            if kwargs['date']:
                kwargs['attributes'][kwargs['field_mapping']['time']] = {}
                kwargs['attributes'][kwargs['field_mapping']['time']]['long_name'] = kwargs['time_longname']
                kwargs['attributes'][kwargs['field_mapping']['time']]['units'] = kwargs['time_units']
        # add default global (file-level) attributes
        if kwargs['title']:
            kwargs['attributes']['ROOT']['title'] = kwargs['title']
        if kwargs['source']:
            kwargs['attributes']['ROOT']['source'] = kwargs['source']
        if kwargs['reference']:
            kwargs['attributes']['ROOT']['reference'] = kwargs['reference']
        # HDF5 dimension variables
        dimensions = []
        dimensions.append('lat')
        dimensions.append('lon')
        # expand dimensions if containing date variables
        if kwargs['date']:
            self.expand_dims()
            dimensions.append('time')
        dims = tuple(kwargs['field_mapping'][key] for key in dimensions)
        # Defining the HDF5 dataset variables
        h5 = {}
        for field,key in kwargs['field_mapping'].items():
            temp = getattr(self,field)
            key = kwargs['field_mapping'][field]
            h5[key] = fileID.create_dataset(key, temp.shape,
                data=temp, dtype=temp.dtype, compression='gzip')
            # filling HDF5 dataset attributes
            for att_name,att_val in kwargs['attributes'][key].items():
                # skip variable attribute if None
                if not att_val:
                    continue
                # skip variable attributes if in list
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    h5[key].attrs[att_name] = att_val
        # add dimensions
        variables = set(kwargs['field_mapping'].keys()) - set(dimensions)
        for field in sorted(variables):
            key = kwargs['field_mapping'][field]
            for i,dim in enumerate(dims):
                h5[key].dims[i].label = dim
                h5[key].dims[i].attach_scale(h5[dim])
            # Dataset contains missing values
            if (self.fill_value is not None):
                h5[key].attrs['_FillValue'] = self.fill_value
        # global attributes of HDF5 file
        for att_name,att_val in kwargs['attributes']['ROOT'].items():
            fileID.attrs[att_name] = att_val
        # add software information
        fileID.attrs['software_reference'] = gravity_toolkit.version.project_name
        fileID.attrs['software_version'] = gravity_toolkit.version.full_version
        # date created
        fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
        # Output HDF5 structure information
        logging.info(str(self.filename))
        logging.info(list(fileID.keys()))
        # Closing the NetCDF file
        fileID.close()

    def to_index(self, filename, file_list, format=None, date=True, **kwargs):
        """
        Write a ``spatial`` object to index of ascii, netCDF4 or HDF5 files

        Parameters
        ----------
        filename: str
            full path of index file to be written
        file_list: list
            list of filenames for each output file
        format: str or NoneType, default None
            format of files in index

                - ``'ascii'``
                - ``'netCDF4'``
                - ``'HDF5'``
        date: bool, default True
            ``spatial`` object contains date information
        verbose: bool, default False
            print file and variable information
        kwargs: dict
            keyword arguments for output writers
        """
        # Write index file of output spatial files
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fid = self.filename.open(mode='w', encoding='utf8')
        # set default verbosity
        kwargs.setdefault('verbose',False)
        # for each file to be in the index
        for i,f in enumerate(file_list):
            # print filename to index
            print(self.compressuser(f), file=fid)
            # index spatial object at i
            s = self.index(i, date=date)
            # write to file
            if (format == 'ascii'):
                # ascii (.txt)
                s.to_ascii(f, date=date, **kwargs)
            elif (format == 'netCDF4'):
                # netcdf (.nc)
                s.to_netCDF4(f, date=date, **kwargs)
            elif (format == 'HDF5'):
                # HDF5 (.H5)
                s.to_HDF5(f, date=date, **kwargs)
        # close the index file
        fid.close()

    def to_file(self, filename, format=None, date=True, **kwargs):
        """
        Write a ``spatial`` object to a specified format

        Parameters
        ----------
        filename: str
            full path of output file
        format: str or NoneType, default None
            file format

                - ``'ascii'``
                - ``'netCDF4'``
                - ``'HDF5'``
        date: bool, default True
            ``spatial`` object contains date information
        verbose: bool, default False
            print file and variable information
        kwargs: dict
            keyword arguments for output writers
        """
        # set default verbosity
        kwargs.setdefault('verbose',False)
        # write to file
        if (format == 'ascii'):
            # ascii (.txt)
            self.to_ascii(filename, date=date, **kwargs)
        elif (format == 'netCDF4'):
            # netcdf (.nc)
            self.to_netCDF4(filename, date=date, **kwargs)
        elif (format == 'HDF5'):
            # HDF5 (.H5)
            self.to_HDF5(filename, date=date, **kwargs)

    def default_field_mapping(self, variables):
        """
        Builds field mappings from a variable list


        Parameters
        ----------
        variables: list
            netCDF4/HDF5 variables names to be mapped

                - ``lonname``
                - ``latname``
                - ``varname``
                - ``timename``

        Returns
        -------
        field_mapping: dict
            Field mappings for netCDF4/HDF5 read and write functions
        """
        # get each variable name and add to field mapping dictionary
        field_mapping = {}
        for i, var in enumerate(['lon', 'lat', 'data', 'time']):
            try:
                field_mapping[var] = copy.copy(variables[i])
            except IndexError as exc:
                pass
        # return the field mapping
        return field_mapping

    def to_masked_array(self):
        """
        Convert a ``spatial`` object to a masked numpy array
        """
        return np.ma.array(self.data, mask=self.mask,
            fill_value=self.fill_value)

    def update_mask(self):
        """
        Update the mask of the ``spatial`` object
        """
        if self.fill_value is not None:
            self.mask |= (self.data == self.fill_value)
            self.mask |= np.isnan(self.data)
            self.data[self.mask] = self.fill_value
            if hasattr(self, 'error'):
                self.error[self.mask] = self.fill_value
        return self

    def copy(self):
        """
        Copy a ``spatial`` object to a new ``spatial`` object
        """
        temp = spatial(fill_value=self.fill_value)
        # copy attributes or update attributes dictionary
        if isinstance(self.attributes, list):
            setattr(temp,'attributes',self.attributes)
        elif isinstance(self.attributes, dict):
            temp.attributes.update(self.attributes)
        # assign variables to self
        var = ['lon','lat','data','mask','error','time','month','filename']
        for key in var:
            try:
                val = getattr(self, key)
                setattr(temp, key, copy.copy(val))
            except AttributeError:
                pass
        # update mask
        temp.replace_masked()
        return temp

    def zeros_like(self):
        """
        Create a ``spatial`` object using the dimensions of another
        """
        temp = spatial(fill_value=self.fill_value)
        # assign variables to self
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        var = ['data','mask','error','time','month']
        for key in var:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.zeros_like(val))
            except AttributeError:
                pass
        # update mask
        temp.replace_masked()
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a ``spatial`` object if non-existent
        """
        # change time dimensions to be iterable
        self.time = np.atleast_1d(self.time)
        self.month = np.atleast_1d(self.month)
        # output spatial with a third dimension
        if (np.ndim(self.data) == 2):
            self.data = self.data[:,:,None]
            # try expanding mask variable
            try:
                self.mask = self.mask[:,:,None]
            except Exception as exc:
                pass
            # try expanding spatial error
            try:
                self.error = self.error[:,:,None]
            except AttributeError as exc:
                pass
        # update mask
        self.update_mask()
        return self

    # PURPOSE: Extend a global matrix
    def extend_matrix(self):
        """
        Extends a global matrix to wrap along longitudes

        Returns
        -------
        temp: float
            extended matrix
        """
        temp = self.copy()
        # shape of the original data object
        ny, nx, *nt = self.shape
        # extended longitude array [x-1,x0,...,xN,xN+1]
        temp.lon = np.zeros((nx+2), dtype=self.lon.dtype)
        temp.lon[0] = self.lon[0] - self.spacing[0]
        temp.lon[1:-1] = self.lon[:]
        temp.lon[-1] = self.lon[-1] + self.spacing[1]
        # attempt to extend possible data variables
        for key in ['data','mask','error']:
            try:
                # get the original data variable
                var = getattr(self, key)
                # extended data matrices along longitude axis
                if (self.ndim == 2):
                    tmp = np.zeros((ny, nx+2), dtype=var.dtype)
                    tmp[:,0] = var[:,-1]
                    tmp[:,1:-1] = var[:,:]
                    tmp[:,-1] = var[:,0]
                elif (self.ndim == 3):
                    var = getattr(self, key)
                    tmp = np.zeros((ny, nx+2, nt[0]), dtype=var.dtype)
                    tmp[:,0,:] = var[:,-1,:]
                    tmp[:,1:-1,:] = var[:,:,:]
                    tmp[:,-1,:] = var[:,0,:]
                # set the output extended data variable
                setattr(temp, key, tmp)
            except Exception as exc:
                pass
        # update mask
        temp.update_mask()
        # return the extended spatial object
        return temp

    def squeeze(self):
        """
        Remove singleton dimensions from a ``spatial`` object
        """
        # squeeze singleton dimensions
        self.time = np.squeeze(self.time)
        self.month = np.squeeze(self.month)
        # attempt to squeeze possible data variables
        for key in ['data','mask','error']:
            try:
                setattr(self, key, np.squeeze(getattr(self, key)))
            except Exception as exc:
                pass
        # update mask
        self.update_mask()
        return self

    def index(self, indice, date=True):
        """
        Subset a ``spatial`` object to specific index

        Parameters
        ----------
        indice: int
            index in matrix for subsetting
        date: bool, default True
            ``spatial`` objects contain date information
        """
        # output spatial object
        temp = spatial(fill_value=self.fill_value)
        # attempt to subset possible data variables
        for key in ['data','mask','error']:
            try:
                tmp = getattr(self, key)
                setattr(temp, key, tmp[:,:,indice].copy())
            except Exception as exc:
                pass
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # subset output dates
        if date:
            temp.time = self.time[indice].copy()
            temp.month = self.month[indice].copy()
        # subset filenames if applicable
        if getattr(self, 'filename'):
            if isinstance(self.filename, (list, tuple, np.ndarray)):
                temp.filename = str(self.filename[indice])
            elif isinstance(self.filename, str):
                temp.filename = copy.copy(self.filename)
        # remove singleton dimensions if importing a single value
        return temp.squeeze()

    def subset(self, months):
        """
        Subset a ``spatial`` object to specific GRACE/GRACE-FO months

        Parameters
        ----------
        months: int
            GRACE/GRACE-FO to subset
        """
        # check if months is an array or a single value
        months = np.atleast_1d(months)
        # number of months
        n = len(months)
        # check that all months are available
        months_check = list(set(months) - set(self.month))
        if months_check:
            m = ','.join([f'{m:03d}' for m in months_check])
            raise IOError(f'GRACE/GRACE-FO months {m} not Found')
        # indices to sort data objects
        months_list = [i for i,m in enumerate(self.month) if m in months]
        # output spatial object
        temp = self.zeros_like()
        # create output spatial object
        temp.data = np.zeros((self.shape[0],self.shape[1],n))
        temp.mask = np.zeros((self.shape[0],self.shape[1],n), dtype=bool)
        # create output spatial error
        try:
            getattr(self, 'error')
            temp.error = np.zeros((self.shape[0],self.shape[1],n))
        except AttributeError:
            pass
        # allocate for output dates
        temp.time = np.zeros((n))
        temp.month = np.zeros((n),dtype=np.int64)
        temp.filename = []
        # for each indice
        for t,i in enumerate(months_list):
            temp.data[:,:,t] = self.data[:,:,i].copy()
            temp.mask[:,:,t] = self.mask[:,:,i].copy()
            try:
                temp.error[:,:,t] = self.error[:,:,i].copy()
            except AttributeError:
                pass
            # copy time dimensions
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
            # subset filename
            if getattr(self, 'filename'):
                if isinstance(self.filename, (list, tuple, np.ndarray)):
                    temp.filename.append(str(self.filename[i]))
                elif isinstance(self.filename, str):
                    temp.filename.append(self.filename)
        # remove singleton dimensions if importing a single value
        return temp.squeeze()

    def offset(self, var):
        """
        Offset a ``spatial`` object by a constant

        Parameters
        ----------
        var: float
            scalar value to which the ``spatial`` object will be offset
        """
        temp = self.copy()
        # offset by a single constant or a time-variable scalar
        if (np.ndim(var) == 0):
            temp.data = self.data + var
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            n = len(var)
            temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
            temp.mask = np.zeros((temp.shape[0],temp.shape[1],n),dtype=bool)
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
        # update mask
        temp.update_mask()
        return temp

    def scale(self, var):
        """
        Multiply a ``spatial`` object by a constant

        Parameters
        ----------
        var: float
            scalar value to which the ``spatial`` object will be multiplied
        """
        temp = self.copy()
        # multiply by a single constant or a time-variable scalar
        if (np.ndim(var) == 0):
            temp.data = var*self.data
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            n = len(var)
            temp.data = np.zeros((temp.shape[0],temp.shape[1],n))
            temp.mask = np.zeros((temp.shape[0],temp.shape[1],n),dtype=bool)
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
        # update mask
        temp.update_mask()
        return temp

    def mean(self, apply=False, indices=Ellipsis):
        """
        Compute mean spatial field and remove from data if specified

        Parameters
        ----------
        apply: bool, default False
            remove the mean field from the input ``spatial`` object
        indices: int, default Ellipsis
            indices of input ``spatial`` object to compute mean
        """
        # output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # create output mean spatial object
        temp.data = np.mean(self.data[:,:,indices],axis=2)
        temp.mask = np.any(self.mask[:,:,indices],axis=2)
        # calculate the mean time
        try:
            val = getattr(self, 'time')
            temp.time = np.mean(val[indices])
        except (AttributeError,TypeError):
            pass
        # calculate the spatial anomalies by removing the mean field
        if apply:
            for i,t in enumerate(self.time):
                self.data[:,:,i] -= temp.data[:,:]
        # update mask
        temp.update_mask()
        return temp

    def flip(self, axis=0):
        """
        Reverse the order of data and dimensions along an axis

        Parameters
        ----------
        axis: int, default 0
            axis to reorder
        """
        # output spatial object
        temp = self.copy()
        # copy dimensions and reverse order
        if (axis == 0):
            temp.lat = temp.lat[::-1].copy()
        elif (axis == 1):
            temp.lon = temp.lon[::-1].copy()
        elif (axis == 2):
            temp.time = temp.time[::-1].copy()
        # attempt to reverse possible data variables
        for key in ['data','mask','error']:
            try:
                setattr(temp, key, np.flip(getattr(self, key), axis=axis))
            except Exception as exc:
                pass
        # update mask
        temp.update_mask()
        return temp

    def transpose(self, axes=None):
        """
        Transpose or permute the axes of a ``spatial`` object

        Parameters
        ----------
        axis: int or NoneType, default None
            order of the output axes
        """
        # output spatial object
        temp = self.copy()
        # attempt to transpose possible data variables
        for key in ['data','mask','error']:
            try:
                setattr(temp, key, np.transpose(getattr(self, key), axes=axes))
            except Exception as exc:
                pass
        # update mask
        temp.update_mask()
        return temp

    def sum(self, power=1):
        """
        Compute summation of a ``spatial`` object

        Parameters
        ----------
        power: int, default 1
            apply a power before calculating summation
        """
        # output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # create output summation spatial object
        temp.data = np.sum(np.power(self.data,power),axis=2)
        temp.mask = np.any(self.mask,axis=2)
        # update mask
        temp.update_mask()
        return temp

    def power(self, power):
        """
        Raise a ``spatial`` object to a power

        Parameters
        ----------
        power: int
            power to which the ``spatial`` object will be raised
        """
        temp = self.copy()
        temp.data = np.power(self.data,power)
        return temp

    def max(self):
        """
        Compute maximum value of a ``spatial`` object
        """
        # output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # create output maximum spatial object
        temp.data = np.max(self.data,axis=2)
        temp.mask = np.any(self.mask,axis=2)
        # update mask
        temp.update_mask()
        return temp

    def min(self):
        """
        Compute minimum value of a ``spatial`` object
        """
        # output spatial object
        temp = spatial(nlon=self.shape[0],nlat=self.shape[1],
            fill_value=self.fill_value)
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # create output minimum spatial object
        temp.data = np.min(self.data,axis=2)
        temp.mask = np.any(self.mask,axis=2)
        # update mask
        temp.update_mask()
        return temp

    def replace_invalid(self, fill_value, mask=None):
        """
        Replace the masked values with a new ``fill_value``

        Parameters
        ----------
        fill_value: float
            Replacement invalid value
        mask: bool or NoneType, default None
            Update the current mask
        """
        # validate current mask
        self.update_mask()
        # update the mask if specified
        if mask is not None:
            if (np.shape(mask) == self.shape):
                self.mask |= mask
            elif (np.ndim(mask) == 2) & (self.ndim == 3):
                # broadcast mask over third dimension
                temp = np.repeat(mask[:,:,np.newaxis],self.shape[2],axis=2)
                self.mask |= temp
        # update the fill value
        self.fill_value = fill_value
        # replace invalid values with new fill value
        self.data[self.mask] = self.fill_value
        if hasattr(self, 'error'):
            self.error[self.mask] = self.fill_value
        return self

    def replace_masked(self):
        """
        Replace the masked values with ``fill_value``
        """
        if (self.fill_value is not None):
            self.data[self.mask] = self.fill_value
        if (self.fill_value is not None) and hasattr(self, 'error'):
            self.error[self.mask] = self.fill_value
        return self

    @property
    def dtype(self):
        """Main data type of ``spatial`` object"""
        return self.data.dtype

    @property
    def spacing(self):
        """Step size of ``spatial`` object ``[longitude,latitude]``
        """
        dlat = np.abs(self.lat[1] - self.lat[0])
        dlon = np.abs(self.lon[1] - self.lon[0])
        return (dlon,dlat)

    @property
    def extent(self):
        """Bounds of ``spatial`` object
        ``[minimum longitude, maximum longitude,
        minimum latitude, maximum latitude]``
        """
        lonmin = np.min(self.lon)
        lonmax = np.max(self.lon)
        latmin = np.min(self.lat)
        latmax = np.max(self.lat)
        return [lonmin, lonmax, latmin, latmax]

    @property
    def shape(self):
        """Dimensions of ``spatial`` object
        """
        return np.shape(self.data)

    @property
    def ndim(self):
        """Number of dimensions in ``spatial`` object
        """
        return np.ndim(self.data)

    def __str__(self):
        """String representation of the ``spatial`` object
        """
        properties = ['gravity_toolkit.spatial']
        extent = ', '.join(map(str, self.extent))
        properties.append(f"    extent: {extent}")
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        if self.month:
            properties.append(f"    start_month: {min(self.month)}")
            properties.append(f"    end_month: {max(self.month)}")
        return '\n'.join(properties)

    def __len__(self):
        """Number of months
        """
        return len(self.month) if np.any(self.month) else 0

    def __iter__(self):
        """Iterate over GRACE/GRACE-FO months
        """
        self.__index__ = 0
        return self

    def __next__(self):
        """Get the next month of data
        """
        # output spatial object
        temp = spatial(fill_value=self.fill_value)
        # subset output spatial field and dates
        try:
            temp.data = self.data[:,:,self.__index__].copy()
            temp.mask = self.mask[:,:,self.__index__].copy()
        except IndexError as exc:
            raise StopIteration from exc
        # subset output spatial time and month
        try:
            temp.time = self.time[self.__index__].copy()
            temp.month = self.month[self.__index__].copy()
        except AttributeError as exc:
            pass
        # subset output spatial error
        try:
            temp.error = self.error[:,:,self.__index__].copy()
        except AttributeError as exc:
            pass
        # subset filename
        if getattr(self, 'filename'):
            if isinstance(self.filename, (list, tuple, np.ndarray)):
                temp.filename = str(self.filename[self.__index__])
            elif isinstance(self.filename, str):
                temp.filename = copy.copy(self.filename)
        # copy dimensions
        temp.lon = self.lon.copy()
        temp.lat = self.lat.copy()
        # add to index
        self.__index__ += 1
        return temp

# PURPOSE: additional routines for the spatial module
# for outputting scaling factor data
class scaling_factors(spatial):
    """
    Inheritance of ``spatial`` class for outputting scaling factors

    Attributes
    ----------
    data: float
        spatial scaling factor data
    error: float
        spatial scaling factor errors
    magnitude: float
        spatial magnitude of the original data
    mask: bool
        spatial grid mask
    x: float
        x-coordinate array
    y: float
        y-coordinate array
    lon: float
        grid longitudes
    lat: float
        grid latitudes
    fill_value: float or NoneType, default None
        invalid value for spatial grid data
    attributes: dict
        attributes of spatial variables
    extent: list, default [None,None,None,None]
        spatial grid bounds
        ``[minimum x, maximum x, minimum y, maximum y]``
    spacing: list, default [None,None]
        grid step size ``[x, y]``
    shape: tuple
        dimensions of spatial object
    ndim: int
        number of dimensions of spatial object
    filename: str
        input or output filename
    """
    np.seterr(invalid='ignore')
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.error = None
        self.magnitude = None

    def from_ascii(self, filename, **kwargs):
        """
        Read a ``scaling_factors`` object from an ascii file

        Parameters
        ----------
        filename: str
            full path of input ascii file
        compression: str or NoneType, default None
            file compression type

                - ``'gzip'``
                - ``'zip'``
                - ``'bytes'``
        spacing: list, default [None,None]
            grid step size ``[longitude,latitude]``
        extent: list, default [None,None,None,None]
            spatial grid bounds
            ``[minimum longitude, maximum longitude,
            minimum latitude, maximum latitude]``
        nlat: int or NoneType, default None
            length of latitude dimension
        nlon: int or NoneType, default None
            length of longitude dimension
        columns: list, default ['lon','lat','kfactor','error','magnitude']
            variable names for each column
        header: int, default 0
            Number of rows of header lines to skip
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('verbose',False)
        kwargs.setdefault('compression',None)
        kwargs.setdefault('spacing',[None,None])
        kwargs.setdefault('nlat',None)
        kwargs.setdefault('nlon',None)
        kwargs.setdefault('extent',[None]*4)
        default_columns = ['lon','lat','kfactor','error','magnitude']
        kwargs.setdefault('columns',default_columns)
        kwargs.setdefault('header',0)
        # open the ascii file and extract contents
        logging.info(str(self.filename))
        if (kwargs['compression'] == 'gzip'):
            # read input ascii data from gzip compressed file and split lines
            with gzip.open(self.filename, mode='r') as f:
                file_contents = f.read().decode('ISO-8859-1').splitlines()
        elif (kwargs['compression'] == 'zip'):
            # read input ascii data from zipped file and split lines
            stem = self.filename.stem
            with zipfile.ZipFile(self.filename) as z:
                file_contents = z.read(stem).decode('ISO-8859-1').splitlines()
        elif (kwargs['compression'] == 'bytes'):
            # read input file object and split lines
            file_contents = self.filename.read().splitlines()
        else:
            # read input ascii file (.txt, .asc) and split lines
            with self.filename.open(mode='r', encoding='utf8') as f:
                file_contents = f.read().splitlines()
        # compile regular expression operator for extracting numerical values
        # from input ascii files of spatial data
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        # output spatial dimensions
        if (None not in kwargs['extent']) and kwargs['nlat'] and kwargs['nlon']:
            extent = kwargs.get('extent')
            self.lat = np.linspace(extent[3], extent[2], kwargs['nlat'])
            self.lon = np.linspace(extent[0], extent[1], kwargs['nlon'])
            dlon = np.abs(self.lon[1] - self.lon[0])
            dlat = np.abs(self.lat[1] - self.lat[0])
        elif (None not in kwargs['extent']) and (None not in kwargs['spacing']):
            extent = kwargs.get('extent')
            dlon, dlat = kwargs.get('spacing')
            self.lat = np.arange(extent[3], extent[2] - dlat, dlat)
            self.lon = np.arange(extent[0], extent[1] + dlon, dlon)
        elif kwargs['nlat'] and kwargs['nlon'] and (None not in kwargs['spacing']):
            dlon, dlat = kwargs.get('spacing')
            self.lat = np.zeros((kwargs['nlat']))
            self.lon = np.zeros((kwargs['nlon']))
        else:
            raise ValueError('Unknown dimensions for input ``spatial`` object')
        # get spatial dimensions
        nlat = len(self.lat)
        nlon = len(self.lon)
        # output spatial data
        self.data = np.zeros((nlat, nlon))
        self.mask = np.zeros((nlat, nlon), dtype=bool)
        # remove time from list of column names if not date
        columns = [c for c in kwargs['columns'] if (c != 'time')]
        # extract spatial data array and convert to matrix
        # for each line in the file
        header = kwargs['header']
        for line in file_contents[header:]:
            # extract columns of interest and assign to dict
            # convert fortran exponentials if applicable
            d = {c:r.replace('D','E') for c,r in zip(columns,rx.findall(line))}
            # convert line coordinates to integers
            ilon = np.int64(np.float64(d['lon'])/dlon)
            ilat = np.int64((90.0-np.float64(d['lat']))//dlat)
            # get scaling factor, error and magnitude
            self.data[ilat,ilon] = np.float64(d['data'])
            self.error[ilat,ilon] = np.float64(d['error'])
            self.magnitude[ilat,ilon] = np.float64(d['magnitude'])
            # set mask
            self.mask[ilat,ilon] = False
            # set latitude and longitude
            self.lon[ilon] = np.float64(d['lon'])
            self.lat[ilat] = np.float64(d['lat'])
        # update mask
        self.update_mask()
        return self

    def to_ascii(self, filename, **kwargs):
        """
        Write a ``scaling_factors`` object to ascii file

        Parameters
        ----------
        filename: str
            full path of output ascii file
        verbose: bool, default False
            Output file and variable information
        """
        self.filename = pathlib.Path(filename).expanduser().absolute()
        # set default verbosity and parameters
        kwargs.setdefault('verbose',False)
        logging.info(str(self.filename))
        # open the output file
        fid = self.filename.open(mode='w', encoding='utf8')
        # write to file for each valid latitude and longitude
        ii,jj = np.nonzero((self.data != self.fill_value) & (~self.mask))
        for i,j in zip(ii,jj):
            print((f'{self.lon[j]:10.4f} {self.lat[i]:10.4f} '
                  f'{self.data[i,j]:12.4f} {self.error[i,j]:12.4f} '
                  f'{self.magnitude[i,j]:12.4f}'), file=fid)
        # close the output file
        fid.close()

    def kfactor(self, var):
        """
        Calculate the scaling factor and scaling factor errors
        from two ``spatial`` or ``scaling_factors`` objects
        following :cite:t:`Landerer:2012kf` and :cite:p:`Hsu:2017hd`

        Parameters
        ----------
        var: obj
            ``spatial`` object to used for scaling

        Returns
        -------
        temp: obj
            scaling factor, scaling error and magnitude
        """
        # copy to not modify original inputs
        temp1 = self.copy()
        temp2 = var.copy()
        # expand dimensions and replace invalid values with 0
        temp1.expand_dims().replace_invalid(0.0)
        temp2.expand_dims().replace_invalid(0.0)
        # dimensions of input spatial object
        nlat, nlon, nt = temp1.shape
        # allocate for scaling factor and scaling factor error
        temp = scaling_factors(nlat=nlat, nlon=nlon, fill_value=0.0)
        temp.data = np.zeros((nlat, nlon))
        temp.error = np.zeros((nlat, nlon))
        # copy latitude and longitude variables
        temp.lon = np.copy(temp1.lon)
        temp.lat = np.copy(temp1.lat)
        # find valid data points and set mask
        temp.mask = np.any(temp1.mask | temp2.mask, axis=2)
        indy,indx = np.nonzero(np.logical_not(temp.mask))
        # calculate point-based scaling factors as centroids
        val1 = np.sum(temp1.data[indy,indx,:]*temp2.data[indy,indx,:],axis=1)
        val2 = np.sum(temp1.data[indy,indx,:]**2,axis=1)
        temp.data[indy,indx] = val1/val2
        # calculate difference between scaled and original
        variance = temp1.scale(temp.data).offset(-temp2.data)
        # calculate scaling factor errors as RMS of variance
        temp.error = np.sqrt((variance.sum(power=2).data)/nt)
        # calculate magnitude of original data
        temp.magnitude = temp2.sum(power=2.0).power(0.5).data[:]
        # update mask
        temp.update_mask()
        # return the scaling factors and scaling factor errors
        return temp

    def update_mask(self):
        """
        Update the mask of the ``scaling_factors`` object
        """
        if self.fill_value is not None:
            self.mask |= (self.data == self.fill_value)
            self.mask |= np.isnan(self.data)
            self.data[self.mask] = self.fill_value
            # replace fill values within scaling factor errors
            if getattr(self, 'error') is not None:
                self.error[self.mask] = self.fill_value
            # replace fill values within scaling factor magnitudes
            if getattr(self, 'magnitude') is not None:
                self.magnitude[self.mask] = self.fill_value
        return self

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
