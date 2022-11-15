#!/usr/bin/env python
u"""
harmonics.py
Written by Tyler Sutterley (06/2023)
Contributions by Hugo Lecomte

Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data

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
    read_GRACE_harmonics.py: read spherical harmonic data from SHM files
    read_gfc_harmonics.py: reads spherical harmonic data from gfc files
    read_ICGEM_harmonics.py: reads gravity model coefficients from GFZ ICGEM
    destripe_harmonics.py: filters spherical harmonics for correlated errors

UPDATE HISTORY:
    Updated 06/2023: fix GRACE/GRACE-FO months in drift function
    Updated 05/2023: use reify decorators for complex form and amplitude
        use pathlib to define and operate on paths
    Updated 03/2023: customizable file-level attributes to netCDF4 and HDF5
        add attributes fetching to the from_dict and to_dict functions
        retrieve all root attributes from HDF5 and netCDF4 datasets
        only attempt to squeeze from final dimension in harmonics objects
        add indexing of filenames to harmonics object iterator
        use copy.copy and not numpy.copy in copy harmonics object function
        convert shape and ndim to harmonic class properties
        improve typing for variables in docstrings
        set case insensitive filename to None if filename is empty
    Updated 02/2023: fix expand case where data is a single degree
        fixed case where maximum spherical harmonic degree is 0
        use monospaced text for harmonics objects in docstrings
    Updated 01/2023: made amplitude a property of the harmonics class
        added property for the complex form of the spherical harmonics
    Updated 12/2022: add software information to output HDF5 and netCDF4
        moved GIA model reader to be an inherited class of harmonics
        make harmonics objects iterable and with length
        added function for creating sized empty harmonic objects
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
        using internal netCDF4 and HDF5 readers and writers
        added function for converting to a python dictionary
        include utf-8 encoding in reads to be windows compliant
        added GIA model reader and drift functions
        include filename attribute when modifying harmonic objects
        allow input ascii files to have additional columns
    Updated 12/2021: logging case_insensitive_filename output for debugging
    Updated 11/2021: kwargs to index, netCDF4 and HDF5 read functions
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: added time-variable gravity data from gfc format
        use functions for converting to and from GRACE months
    Updated 08/2021: added from spherical harmonic model (SHM) format
    Updated 07/2021: fixed gfc format in from file wrapper
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
        use file not found exceptions in case insensitive filename
    Updated 02/2021: added degree amplitude function
        use adjust_months function to fix special cases of GRACE/GRACE-FO months
        added generic reader, generic writer and write to list functions
        generalize ascii, netCDF4 and HDF5 readers and writers
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added verbose option for gfc files
        can calculate spherical harmonic mean over a range of time indices
        will also calculate the mean time and month of a harmonics object
        can create a harmonics object from an open file-like object
    Updated 11/2020: added plotting functions for visualization
    Updated 08/2020: added compression options for ascii, netCDF4 and HDF5 files
    Updated 07/2020: added class docstring and using kwargs for output to file
        added case_insensitive_filename function to search directories
    Updated 06/2020: output list of filenames with from_list()
        zeros_like() creates a new harmonics object with dimensions of another
        add ndim and shape attributes of harmonics objects
    Updated 04/2020: added from_gfc to read static gravity model coefficients
        add to_ascii and iterate over temporal fields in convolve and destripe
        make date optional for harmonic read functions.  add more math functions
        add option to sort if reading from an index or merging a list
        add options to flatten and expand harmonics matrices or arrays
    Written 03/2020
"""
from __future__ import print_function, division

import re
import io
import copy
import gzip
import time
import uuid
import logging
import pathlib
import zipfile
import warnings
import matplotlib
import numpy as np
import gravity_toolkit.version
from gravity_toolkit.time import adjust_months,calendar_to_grace
import scipy as sc
import matplotlib.pyplot as plt
import gravity_toolkit.wavelets as wv
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.hdf5_stokes import hdf5_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes
from gravity_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from gravity_toolkit.destripe_harmonics import destripe_harmonics
from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
from gravity_toolkit.utilities import reify

# attempt imports
try:
    from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
except (ImportError, ModuleNotFoundError) as exc:
    warnings.filterwarnings("module")
    warnings.warn("geoid_toolkit not available", ImportWarning)
try:
    import h5py
except (ImportError, ModuleNotFoundError) as exc:
    warnings.filterwarnings("module")
    warnings.warn("h5py not available", ImportWarning)
try:
    import netCDF4
except (ImportError, ModuleNotFoundError) as exc:
    warnings.filterwarnings("module")
    warnings.warn("netCDF4 not available", ImportWarning)
# ignore warnings
warnings.filterwarnings("ignore")

class harmonics(object):
    """
    Data class for reading, writing and processing spherical harmonic data

    Attributes
    ----------
    lmax: int
        maximum degree of the spherical harmonic field
    mmax: int
        maximum order of the spherical harmonic field
    clm: np.ndarray
        cosine spherical harmonics
    slm: np.ndarray
        sine spherical harmonics
    time: np.ndarray
        time variable of the spherical harmonics
    month: np.ndarray
        GRACE/GRACE-FO months variable of the spherical harmonics
    attributes: dict
        attributes of ``harmonics`` variables
    filename: str
        input or output filename
    flattened: bool
        ``harmonics`` object is compressed into arrays
    """
    np.seterr(invalid='ignore')
    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault('lmax',None)
        kwargs.setdefault('mmax',None)
        # set default class attributes
        self.clm=None
        self.slm=None
        self.time=None
        self.month=None
        self.lmax=kwargs['lmax']
        self.mmax=kwargs['mmax']
        # calculate spherical harmonic degree and order (0 is falsy)
        self.l=np.arange(self.lmax+1) if (self.lmax is not None) else None
        self.m=np.arange(self.mmax+1) if (self.mmax is not None) else None
        self.attributes=dict()
        self.filename=None
        self.flattened=False
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
                    errmsg = f'{filename} not found in file system'
                    raise FileNotFoundError(errmsg)
                self.filename = self.filename.with_name(f.pop())
        # print filename
        logging.debug(self.filename)
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

    def from_ascii(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from an ascii file

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
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date', True)
        kwargs.setdefault('verbose', False)
        kwargs.setdefault('compression', None)
        # open the ascii file and extract contents
        logging.info(self.filename)
        if kwargs['compression'] == 'gzip':
            # read input ascii data from gzip compressed file and split lines
            with gzip.open(self.filename, mode='r') as f:
                file_contents = f.read().decode('ISO-8859-1').splitlines()
        elif kwargs['compression'] == 'zip':
            # read input ascii data from zipped file and split lines
            stem = self.filename.stem
            with zipfile.ZipFile(self.filename) as z:
                file_contents = z.read(stem).decode('ISO-8859-1').splitlines()
        elif kwargs['compression'] == 'bytes':
            # read input file object and split lines
            file_contents = self.filename.read().splitlines()
        else:
            # read input ascii file (.txt, .asc) and split lines
            with open(self.filename, mode='r', encoding='utf8') as f:
                file_contents = f.read().splitlines()
        # compile regular expression operator for extracting numerical values
        # from input ascii files of spherical harmonics
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        # find maximum degree and order of harmonics
        self.lmax = 0
        self.mmax = 0
        # for each line in the file
        for line in file_contents:
            l1,m1,clm1,slm1,*aux = rx.findall(line)
            # convert line degree and order to integers
            l1,m1 = np.array([l1,m1],dtype=np.int64)
            self.lmax = np.copy(l1) if (l1 > self.lmax) else self.lmax
            self.mmax = np.copy(m1) if (m1 > self.mmax) else self.mmax
        # output spherical harmonics data
        self.clm = np.zeros((self.lmax+1,self.mmax+1))
        self.slm = np.zeros((self.lmax+1,self.mmax+1))
        # if the ascii file contains date variables
        if kwargs['date']:
            self.time = np.float64(aux[0])
            self.month = np.int64(calendar_to_grace(self.time))
            # adjust months to fix special cases if necessary
            self.month = adjust_months(self.month)
        # extract harmonics and convert to matrix
        # for each line in the file
        for line in file_contents:
            l1,m1,clm1,slm1,*aux = rx.findall(line)
            # convert line degree and order to integers
            ll,mm = np.array([l1,m1],dtype=np.int64)
            # convert fortran exponentials if applicable
            self.clm[ll,mm] = np.float64(clm1.replace('D','E'))
            self.slm[ll,mm] = np.float64(slm1.replace('D','E'))
        # assign degree and order fields
        self.update_dimensions()
        return self

    def from_netCDF4(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from a netCDF4 file

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
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date', True)
        kwargs.setdefault('verbose', False)
        kwargs.setdefault('compression', None)
        # Open the NetCDF4 file for reading
        if kwargs['compression'] == 'gzip':
            # read as in-memory (diskless) netCDF4 dataset
            with gzip.open(self.filename, mode='r') as f:
                fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
        elif kwargs['compression'] == 'zip':
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
        elif kwargs['compression'] == 'bytes':
            # read as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=filename.read())
        else:
            # read netCDF4 dataset
            fileID = netCDF4.Dataset(self.filename, mode='r')
        # Output NetCDF file information
        logging.info(fileID.filepath())
        logging.info(list(fileID.variables.keys()))
        # read flattened spherical harmonics
        temp = harmonics()
        temp.filename = copy.copy(self.filename)
        # create list of variables to retrieve
        fields = ['l','m','clm','slm']
        # retrieve date variables if specified
        if kwargs['date']:
            fields.extend(['time','month'])
        # Getting the data from each NetCDF variable
        for field in fields:
            setattr(temp, field, fileID.variables[field][:].copy())
        # calculate maximum degree and order
        temp.lmax = np.max(temp.l)
        temp.mmax = np.max(temp.m)
        # expand the spherical harmonics to dimensions
        self = temp.expand(date=kwargs['date'])
        # adjust months to fix special cases if necessary
        self.month = adjust_months(self.month)
        # attributes of clm/slm and included variables
        for key in fields:
            # attempt to get attribute for variable
            try:
                self.attributes[key] = [
                    fileID.variables[key].units,
                    fileID.variables[key].long_name
                    ]
            except (KeyError,ValueError,AttributeError):
                pass
        # get global netCDF4 attributes
        self.attributes['ROOT'] = {}
        for att_name in fileID.ncattrs():
            self.attributes['ROOT'][att_name] = fileID.getncattr(att_name)
        # Closing the NetCDF file
        fileID.close()
        # remove singleton dimensions and
        # assign degree and order fields
        self.squeeze(update_dimensions=True)
        return self

    def from_HDF5(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from a HDF5 file

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
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date', True)
        kwargs.setdefault('verbose', False)
        kwargs.setdefault('compression', None)
        # Open the HDF5 file for reading
        if kwargs['compression'] == 'gzip':
            # read gzip compressed file and extract into in-memory file object
            with gzip.open(self.filename, mode='r') as f:
                fid = io.BytesIO(f.read())
            # set filename of BytesIO object
            fid.filename = self.filename.name
            # rewind to start of file
            fid.seek(0)
            # read as in-memory (diskless) HDF5 dataset from BytesIO object
            fileID = h5py.File(fid, mode='r')
        elif kwargs['compression'] == 'zip':
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
        elif kwargs['compression'] == 'bytes':
            # read as in-memory (diskless) HDF5 dataset
            fileID = h5py.File(self.filename, mode='r')
        else:
            # read HDF5 dataset
            fileID = h5py.File(self.filename, mode='r')
        # Output HDF5 file information
        logging.info(fileID.filename)
        logging.info(list(fileID.keys()))
        # read flattened spherical harmonics
        temp = harmonics()
        temp.filename = copy.copy(self.filename)
        # create list of variables to retrieve
        fields = ['l','m','clm','slm']
        # retrieve date variables if specified
        if kwargs['date']:
            fields.extend(['time','month'])
        # Getting the data from each HDF5 variable
        for field in fields:
            setattr(temp, field, fileID[field][:].copy())
        # calculate maximum degree and order
        temp.lmax = np.max(temp.l)
        temp.mmax = np.max(temp.m)
        # expand the spherical harmonics to dimensions
        self = temp.expand(date=kwargs['date'])
        # adjust months to fix special cases if necessary
        self.month = adjust_months(self.month)
        # attributes of clm/slm and included variables
        for key in fields:
            # attempt to get attribute for variable
            try:
                self.attributes[key] = [
                    fileID[key].attrs['units'],
                    fileID[key].attrs['long_name']
                    ]
            except (KeyError, AttributeError):
                pass
        # get global HDF5 attributes
        self.attributes['ROOT'] = {}
        for att_name,att_val in fileID.attrs.items():
            self.attributes['ROOT'][att_name] = att_val
        # Closing the HDF5 file
        fileID.close()
        # remove singleton dimensions and
        # assign degree and order fields
        self.squeeze(update_dimensions=True)
        return self

    def from_gfc(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from a gfc gravity model file

        Parameters
        ----------
        filename: str
            full path of input gfc file
        date: bool, default True
            gfc file has date information
        tide: str or NoneType, default None
            permanent tide system of output gravity fields

                - ``'tide_free'``: no permanent direct and indirect tidal potentials
                - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
                - ``'zero_tide'``: permanent direct tidal potential removed
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default parameters
        kwargs.setdefault('date',False)
        kwargs.setdefault('tide',None)
        kwargs.setdefault('verbose',False)
        # read data from gfc file
        if kwargs['date']:
            Ylms = read_gfc_harmonics(self.filename,
                TIDE=kwargs['tide'])
        else:
            Ylms = read_ICGEM_harmonics(self.filename,
                TIDE=kwargs['tide'])
        # Output file information
        logging.info(self.filename)
        logging.info(list(Ylms.keys()))
        # copy variables for gravity model
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.lmax = np.int64(Ylms['max_degree'])
        self.mmax = np.int64(Ylms['max_degree'])
        # copy date variables
        if kwargs['date']:
            self.time = Ylms['time'].copy()
            # adjust months to fix special cases if necessary
            self.month = adjust_months(Ylms['month'])
        # geophysical parameters of gravity model
        self.GM = np.float64(Ylms['earth_gravity_constant'])
        self.R = np.float64(Ylms['radius'])
        self.tide = Ylms['tide_system']
        # assign degree and order fields
        self.update_dimensions()
        return self

    def from_SHM(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from a spherical harmonic model file

        Parameters
        ----------
        filename: str
            full path of input spherical harmonic model file
        MMAX: int or NoneType, default None
            Maximum order of spherical harmonics
        POLE_TIDE: bool, default False
            correct GSM data for pole tide drift
        verbose: bool, default False
            print file and variable information
        """
        # set filename
        self.case_insensitive_filename(filename)
        # set default keyword arguments
        kwargs.setdefault('verbose',False)
        # read data from SHM file
        Ylms = read_GRACE_harmonics(self.filename, self.lmax, **kwargs)
        # Output file information
        logging.info(self.filename)
        logging.info(list(Ylms.keys()))
        # copy variables for gravity model
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.time = Ylms['time'].copy()
        self.month = np.int64(calendar_to_grace(self.time))
        # copy header information for gravity model
        self.header = Ylms['header']
        # assign degree and order fields
        self.update_dimensions()
        return self

    def from_index(self, filename, **kwargs):
        """
        Read a ``harmonics`` object from an index of ascii, netCDF4 or HDF5 files

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
            sort ``harmonics`` objects by date information
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
        # Read index file of input spherical harmonics
        with open(self.filename, mode='r', encoding='utf8') as f:
            file_list = [l for l in f.read().splitlines() if parser.match(l)]
        # create a list of harmonic objects
        h = []
        # for each file in the index
        for i,f in enumerate(file_list):
            if kwargs['format'] == 'ascii':
                # ascii (.txt)
                h.append(harmonics().from_ascii(f, date=kwargs['date']))
            elif kwargs['format'] == 'netCDF4':
                # netcdf (.nc)
                h.append(harmonics().from_netCDF4(f, date=kwargs['date']))
            elif kwargs['format'] == 'HDF5':
                # HDF5 (.H5)
                h.append(harmonics().from_HDF5(f, date=kwargs['date']))
        # create a single harmonic object from the list
        return self.from_list(h,date=kwargs['date'],sort=kwargs['sort'])

    def from_list(self, object_list, **kwargs):
        """
        Build a sorted ``harmonics`` object from a list of
        other ``harmonics`` objects

        Parameters
        ----------
        object_list: list
            list of ``harmonics`` objects to be merged
        date: bool, default True
            files contains date information
        sort: bool, default True
            sort ``harmonics`` objects by date information
        clear: bool, default True
            clear the list of ``harmonics`` objects from memory
        """
        # set default keyword arguments
        kwargs.setdefault('date',True)
        kwargs.setdefault('sort',True)
        kwargs.setdefault('clear',False)
        # number of harmonic objects in list
        n = len(object_list)
        # indices to sort data objects if harmonics list contain dates
        if kwargs['date'] and kwargs['sort']:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
        else:
            list_sort = np.arange(n)
        # truncate to maximum degree and order
        self.lmax = np.min([d.lmax for d in object_list])
        self.mmax = np.min([d.mmax for d in object_list])
        # create output harmonics
        self.clm = np.zeros((self.lmax+1,self.mmax+1,n))
        self.slm = np.zeros((self.lmax+1,self.mmax+1,n))
        # create list of files
        self.filename = []
        # output dates
        if kwargs['date']:
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int64)
        # for each indice
        for t,i in enumerate(list_sort):
            self.clm[:,:,t] = object_list[i].clm[:self.lmax+1,:self.mmax+1]
            self.slm[:,:,t] = object_list[i].slm[:self.lmax+1,:self.mmax+1]
            if kwargs['date']:
                self.time[t] = np.atleast_1d(object_list[i].time)
                self.month[t] = np.atleast_1d(object_list[i].month)
            # append filename to list
            if getattr(object_list[i], 'filename'):
                self.filename.append(object_list[i].filename)
        # adjust months to fix special cases if necessary
        if kwargs['date']:
            self.month = adjust_months(self.month)
        # assign degree and order fields
        self.update_dimensions()
        # clear the input list to free memory
        if kwargs['clear']:
            object_list = None
        # return the single harmonic object
        return self

    def from_file(self, filename, format=None, date=True, **kwargs):
        """
        Read a ``harmonics`` object from a specified format

        Parameters
        ----------
        filename: str
            full path of input file
        format: str or NoneType, default None
            file format

                - ``'ascii'``
                - ``'netCDF4'``
                - ``'HDF5'``
                - ``'gfc'``
                - ``'SHM'``
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
        kwargs.setdefault('verbose', False)
        # read from file
        if format == 'ascii':
            # ascii (.txt)
            return harmonics().from_ascii(filename, date=date, **kwargs)
        elif format == 'netCDF4':
            # netcdf (.nc)
            return harmonics().from_netCDF4(filename, date=date, **kwargs)
        elif format == 'HDF5':
            # HDF5 (.H5)
            return harmonics().from_HDF5(filename, date=date, **kwargs)
        elif format == 'gfc':
            # ICGEM gravity model (.gfc)
            return harmonics().from_gfc(filename, **kwargs)
        elif format == 'SHM':
            # spherical harmonic model
            return harmonics().from_SHM(filename, self.lmax, **kwargs)

    def from_dict(self, d, **kwargs):
        """
        Convert a ``dict`` object to a ``harmonics`` object

        Parameters
        ----------
        d: dict
            dictionary object to be converted
        """
        # assign dictionary variables to self
        for key in ['l', 'm', 'clm', 'slm', 'time', 'month']:
            try:
                setattr(self, key, d[key].copy())
            except (AttributeError, KeyError):
                pass
        # maximum degree and order
        self.lmax = np.max(d['l'])
        self.mmax = np.max(d['m'])
        # add attributes to root if in dictionary
        self.attributes['ROOT'] = d.get('attributes')
        # assign degree and order fields
        self.update_dimensions()
        return self

    def to_ascii(self, filename, date=True, **kwargs):
        """
        Write a ``harmonics`` object to ascii file

        Parameters
        ----------
        filename: str
            full path of output ascii file
        date: bool, default True
            ``harmonics`` objects contain date information
        verbose: bool, default False
            Output file and variable information
        """
        self.filename = pathlib.Path(filename).expanduser().absolute()
        # set default verbosity
        kwargs.setdefault('verbose',False)
        logging.info(self.filename)
        # open the output file
        fid = open(self.filename, mode='w', encoding='utf8')
        if date:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e} {4:10.4f}'
        else:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e}'
        # write to file for each spherical harmonic degree and order
        for m in range(0, self.mmax+1):
            for l in range(m, self.lmax+1):
                args = (l, m, self.clm[l,m], self.slm[l,m], self.time)
                print(file_format.format(*args), file=fid)
        # close the output file
        fid.close()

    def to_netCDF4(self, filename, **kwargs):
        """
        Write a ``harmonics`` object to netCDF4 file

        Parameters
        ----------
        filename: str
            full path of output netCDF4 file
        units: str, default: 'Geodesy_Normalization'
            spherical harmonic units
        time_units: str, default 'years'
            time variable units
        time_longname: str, default 'Date_in_Decimal_Years'
            time variable description
        months_name: str, default 'month'
            name of months variable
        months_units: str, default 'number'
            months variable units
        months_longname: str, default 'GRACE_month'
            months variable description
        field_mapping: dict, default {}
            mapping between input variables and output netCDF4
        attributes: dict, default {}
            output netCDF4 variable and file-level attributes
        title: str or NoneType, default None
            title attribute of dataset
        source: str or NoneType, default None
            source attribute of dataset
        reference: str or NoneType, default None
            reference attribute of dataset
        date: bool, default True
            ``harmonics`` objects contain date information
        clobber: bool, default True
            Overwrite an existing netCDF4 file
        verbose: bool, default False
            Output file and variable information
        """
        # set default keyword arguments
        kwargs.setdefault('units','Geodesy_Normalization')
        kwargs.setdefault('time_units','years')
        kwargs.setdefault('time_longname','Date_in_Decimal_Years')
        kwargs.setdefault('months_name','month')
        kwargs.setdefault('months_units','number')
        kwargs.setdefault('months_longname','GRACE_month')
        kwargs.setdefault('field_mapping',{})
        attributes = self.attributes.get('ROOT') or {}
        kwargs.setdefault('attributes',dict(ROOT=attributes))
        kwargs.setdefault('title',None)
        kwargs.setdefault('source',None)
        kwargs.setdefault('reference',None)
        kwargs.setdefault('date',True)
        kwargs.setdefault('clobber',True)
        kwargs.setdefault('verbose',False)
        # setting NetCDF clobber attribute
        clobber = 'w' if kwargs['clobber'] else 'a'
        # opening netCDF file for writing
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fileID = netCDF4.Dataset(self.filename, clobber, format="NETCDF4")
        # flatten harmonics
        temp = self.flatten(date=kwargs['date'])
        # mapping between output keys and netCDF4 variable names
        if not kwargs['field_mapping']:
            kwargs['field_mapping']['l'] = 'l'
            kwargs['field_mapping']['m'] = 'm'
            kwargs['field_mapping']['clm'] = 'clm'
            kwargs['field_mapping']['slm'] = 'slm'
            if kwargs['date']:
                kwargs['field_mapping']['time'] = 'time'
                kwargs['field_mapping']['month'] = kwargs['months_name']
        # create attributes dictionary for output variables
        if not all(key in kwargs['attributes'] for key in kwargs['field_mapping']):
            # Defining attributes for degree and order
            kwargs['attributes'][kwargs['field_mapping']['l']] = {}
            kwargs['attributes'][kwargs['field_mapping']['l']]['long_name'] = 'spherical_harmonic_degree'
            kwargs['attributes'][kwargs['field_mapping']['l']]['units'] = 'Wavenumber'
            kwargs['attributes'][kwargs['field_mapping']['m']] = {}
            kwargs['attributes'][kwargs['field_mapping']['m']]['long_name'] = 'spherical_harmonic_order'
            kwargs['attributes'][kwargs['field_mapping']['m']]['units'] = 'Wavenumber'
            # Defining attributes for dataset
            kwargs['attributes'][kwargs['field_mapping']['clm']] = {}
            kwargs['attributes'][kwargs['field_mapping']['clm']]['long_name'] = 'cosine_spherical_harmonics'
            kwargs['attributes'][kwargs['field_mapping']['clm']]['units'] = kwargs['units']
            kwargs['attributes'][kwargs['field_mapping']['slm']] = {}
            kwargs['attributes'][kwargs['field_mapping']['slm']]['long_name'] = 'sine_spherical_harmonics'
            kwargs['attributes'][kwargs['field_mapping']['slm']]['units'] = kwargs['units']
            # Defining attributes for date if applicable
            if kwargs['date']:
                # attributes for date and month (or integer date)
                kwargs['attributes'][kwargs['field_mapping']['time']] = {}
                kwargs['attributes'][kwargs['field_mapping']['time']]['long_name'] = kwargs['time_longname']
                kwargs['attributes'][kwargs['field_mapping']['time']]['units'] = kwargs['time_units']
                kwargs['attributes'][kwargs['field_mapping']['month']] = {}
                kwargs['attributes'][kwargs['field_mapping']['month']]['long_name'] = kwargs['months_longname']
                kwargs['attributes'][kwargs['field_mapping']['month']]['units'] = kwargs['months_units']
        # add default global (file-level) attributes
        if kwargs['title']:
            kwargs['attributes']['ROOT']['title'] = kwargs['title']
        if kwargs['source']:
            kwargs['attributes']['ROOT']['source'] = kwargs['source']
        if kwargs['reference']:
            kwargs['attributes']['ROOT']['reference'] = kwargs['reference']
        # netCDF4 dimension variables
        dimensions = []
        dimensions.append('lm')
        fileID.createDimension('lm', len(temp.l))
        # defining netCDF4 temporal dimension
        if kwargs['date']:
            temp.expand_dims(update_dimensions=False)
            dimensions.append('time')
            fileID.createDimension('time', len(temp.time))
        # defining and filling the netCDF variables
        nc = {}
        for field,key in kwargs['field_mapping'].items():
            val = getattr(temp, field)
            if field in ('l','m'):
                dims = (dimensions[0],)
            elif field in ('time','month'):
                dims = (dimensions[1],)
            else:
                dims = tuple(dimensions)
            # create netCDF4 variable
            nc[key] = fileID.createVariable(key, val.dtype, dims)
            nc[key][:] = val[:]
            # filling netCDF dataset attributes
            for att_name,att_val in kwargs['attributes'][key].items():
                # skip variable attribute if None
                if not att_val:
                    continue
                # skip variable attributes if in list
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    nc[key].setncattr(att_name, att_val)
        # global attributes of NetCDF4 file
        for att_name,att_val in kwargs['attributes']['ROOT'].items():
            fileID.setncattr(att_name, att_val)
        # add software information
        fileID.software_reference = gravity_toolkit.version.project_name
        fileID.software_version = gravity_toolkit.version.full_version
        # date created
        fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
        # Output netCDF structure information
        logging.info(self.filename)
        logging.info(list(fileID.variables.keys()))
        # Closing the netCDF file
        fileID.close()

    def to_HDF5(self, filename, **kwargs):
        """
        Write a ``harmonics`` object to HDF5 file

        Parameters
        ----------
        filename: str
            full path of output HDF5 file
        units: str, default: 'Geodesy_Normalization'
            spherical harmonic units
        time_units: str, default 'years'
            time variable units
        time_longname: str, default 'Date_in_Decimal_Years'
            time variable description
        months_name: str, default 'month'
            name of months variable
        months_units: str, default 'number'
            months variable units
        months_longname: str, default 'GRACE_month'
            months variable description
        field_mapping: dict, default {}
            mapping between input variables and output HDF5
        attributes: dict, default {}
            output HDF5 variable and file-level attributes
        title: str or NoneType, default None
            description attribute of dataset
        source: str or NoneType, default None
            source attribute of dataset
        reference: str or NoneType, default None
            reference attribute of dataset
        date: bool, default True
            ``harmonics`` objects contain date information
        clobber: bool, default True
            Overwrite an existing HDF5 file
        verbose: bool, default False
            Output file and variable information
        """
        # set default keyword arguments
        kwargs.setdefault('units','Geodesy_Normalization')
        kwargs.setdefault('time_units','years')
        kwargs.setdefault('time_longname','Date_in_Decimal_Years')
        kwargs.setdefault('months_name','month')
        kwargs.setdefault('months_units','number')
        kwargs.setdefault('months_longname','GRACE_month')
        kwargs.setdefault('field_mapping',{})
        attributes = self.attributes.get('ROOT') or {}
        kwargs.setdefault('attributes',dict(ROOT=attributes))
        kwargs.setdefault('title',None)
        kwargs.setdefault('source',None)
        kwargs.setdefault('reference',None)
        kwargs.setdefault('date',True)
        kwargs.setdefault('clobber',True)
        kwargs.setdefault('verbose',False)
        # setting HDF5 clobber attribute
        clobber = 'w' if kwargs['clobber'] else 'w-'
        # opening HDF5 file for writing
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fileID = h5py.File(self.filename, clobber)
        # flatten harmonics
        temp = self.flatten(date=kwargs['date'])
        if kwargs['date']:
            temp.expand_dims(update_dimensions=False)
        # mapping between output keys and HDF5 variable names
        if not kwargs['field_mapping']:
            kwargs['field_mapping']['l'] = 'l'
            kwargs['field_mapping']['m'] = 'm'
            kwargs['field_mapping']['clm'] = 'clm'
            kwargs['field_mapping']['slm'] = 'slm'
            if kwargs['date']:
                kwargs['field_mapping']['time'] = 'time'
                kwargs['field_mapping']['month'] = kwargs['months_name']
        # create attributes dictionary for output variables
        if not all(key in kwargs['attributes'] for key in kwargs['field_mapping']):
            # Defining attributes for degree and order
            kwargs['attributes'][kwargs['field_mapping']['l']] = {}
            kwargs['attributes'][kwargs['field_mapping']['l']]['long_name'] = 'spherical_harmonic_degree'
            kwargs['attributes'][kwargs['field_mapping']['l']]['units'] = 'Wavenumber'
            kwargs['attributes'][kwargs['field_mapping']['m']] = {}
            kwargs['attributes'][kwargs['field_mapping']['m']]['long_name'] = 'spherical_harmonic_order'
            kwargs['attributes'][kwargs['field_mapping']['m']]['units'] = 'Wavenumber'
            # Defining attributes for dataset
            kwargs['attributes'][kwargs['field_mapping']['clm']] = {}
            kwargs['attributes'][kwargs['field_mapping']['clm']]['long_name'] = 'cosine_spherical_harmonics'
            kwargs['attributes'][kwargs['field_mapping']['clm']]['units'] = kwargs['units']
            kwargs['attributes'][kwargs['field_mapping']['slm']] = {}
            kwargs['attributes'][kwargs['field_mapping']['slm']]['long_name'] = 'sine_spherical_harmonics'
            kwargs['attributes'][kwargs['field_mapping']['slm']]['units'] = kwargs['units']
            # Defining attributes for date if applicable
            if kwargs['date']:
                # attributes for date and month (or integer date)
                kwargs['attributes'][kwargs['field_mapping']['time']] = {}
                kwargs['attributes'][kwargs['field_mapping']['time']]['long_name'] = kwargs['time_longname']
                kwargs['attributes'][kwargs['field_mapping']['time']]['units'] = kwargs['time_units']
                kwargs['attributes'][kwargs['field_mapping']['month']] = {}
                kwargs['attributes'][kwargs['field_mapping']['month']]['long_name'] = kwargs['months_longname']
                kwargs['attributes'][kwargs['field_mapping']['month']]['units'] = kwargs['months_units']
        # add default global (file-level) attributes
        if kwargs['title']:
            kwargs['attributes']['ROOT']['title'] = kwargs['title']
        if kwargs['source']:
            kwargs['attributes']['ROOT']['source'] = kwargs['source']
        if kwargs['reference']:
            kwargs['attributes']['ROOT']['reference'] = kwargs['reference']
        # Defining the HDF5 dataset variables
        h5 = {}
        for field,key in kwargs['field_mapping'].items():
            val = getattr(temp, field)
            h5[key] = fileID.create_dataset(key, val.shape,
                data=val, dtype=val.dtype, compression='gzip')
            # filling HDF5 dataset attributes
            for att_name,att_val in kwargs['attributes'][key].items():
                # skip variable attribute if None
                if not att_val:
                    continue
                # skip variable attributes if in list
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    h5[key].attrs[att_name] = att_val
        # global attributes of HDF5 file
        for att_name,att_val in kwargs['attributes']['ROOT'].items():
            fileID.attrs[att_name] = att_val
        # add software information
        fileID.attrs['software_reference'] = gravity_toolkit.version.project_name
        fileID.attrs['software_version'] = gravity_toolkit.version.full_version
        # date created
        fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
        # Output HDF5 structure information
        logging.info(self.filename)
        logging.info(list(fileID.keys()))
        # Closing the HDF5 file
        fileID.close()

    def to_index(self, filename, file_list, format=None, date=True, **kwargs):
        """
        Write a ``harmonics`` object to index of ascii, netCDF4 or HDF5 files

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
            ``harmonics`` object contains date information
        verbose: bool, default False
            print file and variable information
        kwargs: dict
            keyword arguments for output writers
        """
        # Write index file of output spherical harmonics
        self.filename = pathlib.Path(filename).expanduser().absolute()
        fid = open(self.filename, mode='w', encoding='utf8')
        # set default verbosity
        kwargs.setdefault('verbose',False)
        # for each file to be in the index
        for i,f in enumerate(file_list):
            # print filename to index
            print(self.compressuser(f), file=fid)
            # index harmonics object at i
            h = self.index(i, date=date)
            # write to file
            if format == 'ascii':
                # ascii (.txt)
                h.to_ascii(f, date=date, **kwargs)
            elif format == 'netCDF4':
                # netcdf (.nc)
                h.to_netCDF4(f, date=date, **kwargs)
            elif format == 'HDF5':
                # HDF5 (.H5)
                h.to_HDF5(f, date=date, **kwargs)
        # close the index file
        fid.close()

    def to_file(self, filename, format=None, date=True, **kwargs):
        """
        Write a ``harmonics`` object to a specified format

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
            ``harmonics`` object contains date information
        verbose: bool, default False
            print file and variable information
        kwargs: dict
            keyword arguments for output writers
        """
        # set default verbosity
        kwargs.setdefault('verbose',False)
        # write to file
        if format == 'ascii':
            # ascii (.txt)
            self.to_ascii(filename, date=date, **kwargs)
        elif format == 'netCDF4':
            # netcdf (.nc)
            self.to_netCDF4(filename, date=date, **kwargs)
        elif format == 'HDF5':
            # HDF5 (.H5)
            self.to_HDF5(filename, date=date, **kwargs)

    def to_dict(self):
        """
        Convert a ``harmonics`` object to a ``dict`` object

        Returns
        -------
        d: dict
            converted dictionary object
        """
        # assign dictionary variables from self
        d = {}
        for key in ['l','m','clm','slm','time','month','attributes']:
            try:
                d[key] = getattr(self, key)
            except (AttributeError, KeyError):
                pass
        # return the dictionary object
        return d

    def to_masked_array(self):
        """
        Convert a ``harmonics`` object to a masked numpy array
        """
        # assign degree and order fields
        self.update_dimensions()
        # verify dimensions and get shape
        ndim_prev = np.copy(self.ndim)
        self.expand_dims()
        l1, m1, nt = self.shape
        # create single triangular matrices with harmonics
        Ylms = np.ma.zeros((self.lmax + 1, 2*self.lmax + 1, nt))
        Ylms.mask = np.ones((self.lmax + 1, 2*self.lmax + 1, nt),dtype=bool)
        for m in range(-self.mmax, self.mmax + 1):
            mm = np.abs(m)
            for l in range(mm, self.lmax + 1):
                if m < 0:
                    Ylms.data[l, self.lmax+m, :] = self.slm[l, mm, :]
                    Ylms.mask[l, self.lmax+m, :] = False
                else:
                    Ylms.data[l, self.lmax+m, :] = self.clm[l,mm,:]
                    Ylms.mask[l, self.lmax+m, :] = False
        # reshape to previous
        if self.ndim != ndim_prev:
            self.squeeze()
        # return the triangular matrix
        return Ylms

    def update_dimensions(self):
        """
        Update the dimension variables of the ``harmonics`` object
        """
        # calculate spherical harmonic degree and order (0 is falsy)
        self.l = np.arange(self.lmax + 1) if (self.lmax is not None) else None
        self.m = np.arange(self.mmax + 1) if (self.mmax is not None) else None
        return self

    def add(self, temp):
        """
        Add two ``harmonics`` objects

        Parameters
        ----------
        temp: obj
            harmonic object to be added
        """
        # assign degree and order fields
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if self.ndim == 2:
            self.clm[:l1, :m1] += temp.clm[:l1, :m1]
            self.slm[:l1, :m1] += temp.slm[:l1, :m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1, :m1, i] += temp.clm[:l1, :m1]
                self.slm[:l1, :m1, i] += temp.slm[:l1, :m1]
        else:
            old_month = self.month
            exclude1 = set(self.month) - set(temp.month)

            self.month = np.array(list(sorted(set(self.month) - exclude1)))
            self.time = np.array([self.time[i] for i in range(len(self.time)) if not (old_month[i] in exclude1)])

            for i in range(len(old_month)):
                for j in range(len(temp.month)):
                    if old_month[i] == temp.month[j]:
                        self.clm[:l1, :m1, i] += temp.clm[:l1, :m1, j]
                        self.slm[:l1, :m1, i] += temp.slm[:l1, :m1, j]

            to_keep = []
            for i in range(len(old_month)):
                if not(old_month[i] in exclude1):
                    to_keep.append(i)
            self.clm = self.clm[:, :, to_keep]
            self.slm = self.slm[:, :, to_keep]
        return self

    def subtract(self, temp):
        """
        Subtract one ``harmonics`` object from another

        Parameters
        ----------
        temp: obj
            harmonic object to be subtracted
        """
        # assign degree and order fields
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if self.ndim == 2:
            self.clm[:l1, :m1] -= temp.clm[:l1, :m1]
            self.slm[:l1, :m1] -= temp.slm[:l1, :m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1, :m1, i] -= temp.clm[:l1, :m1]
                self.slm[:l1, :m1, i] -= temp.slm[:l1, :m1]
        else:
            old_month = self.month
            exclude1 = set(self.month) - set(temp.month)

            self.month = np.array(list(sorted(set(self.month) - exclude1)))
            self.time = np.array([self.time[i] for i in range(len(self.time)) if not (old_month[i] in exclude1)])

            for i in range(len(old_month)):
                for j in range(len(temp.month)):
                    if old_month[i] == temp.month[j]:
                        self.clm[:l1, :m1, i] -= temp.clm[:l1, :m1, j]
                        self.slm[:l1, :m1, i] -= temp.slm[:l1, :m1, j]

            to_keep = []
            for i in range(len(old_month)):
                if not(old_month[i] in exclude1):
                    to_keep.append(i)
            self.clm = self.clm[:, :, to_keep]
            self.slm = self.slm[:, :, to_keep]
        return self

    def multiply(self, temp):
        """
        Multiply two ``harmonics`` objects

        Parameters
        ----------
        temp: obj
            harmonic object to be multiplied
        """
        # assign degree and order fields
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax + 1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax + 1 if (temp.mmax > self.mmax) else temp.mmax+1
        if self.ndim == 2:
            self.clm[:l1, :m1] *= temp.clm[:l1, :m1]
            self.slm[:l1, :m1] *= temp.slm[:l1, :m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1, :m1, i] *= temp.clm[:l1, :m1]
                self.slm[:l1, :m1, i] *= temp.slm[:l1, :m1]
        else:
            old_month = self.month
            exclude1 = set(self.month) - set(temp.month)

            self.month = np.array(list(sorted(set(self.month) - exclude1)))
            self.time = np.array([self.time[i] for i in range(len(self.time)) if not (old_month[i] in exclude1)])

            for i in range(len(old_month)):
                for j in range(len(temp.month)):
                    if old_month[i] == temp.month[j]:
                        self.clm[:l1, :m1, i] *= temp.clm[:l1, :m1, j]
                        self.slm[:l1, :m1, i] *= temp.slm[:l1, :m1, j]

            to_keep = []
            for i in range(len(old_month)):
                if not (old_month[i] in exclude1):
                    to_keep.append(i)
            self.clm = self.clm[:, :, to_keep]
            self.slm = self.slm[:, :, to_keep]
        return self

    def divide(self, temp):
        """
        Divide one ``harmonics`` object from another

        Parameters
        ----------
        temp: obj
            harmonic object to be divided
        """
        # assign degree and order fields
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax + 1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax + 1 if (temp.mmax > self.mmax) else temp.mmax+1
        # indices for cosine spherical harmonics (including zonals)
        lc,mc = np.tril_indices(l1, m=m1)
        # indices for sine spherical harmonics (excluding zonals)
        m0 = np.nonzero(mc != 0)
        ls,ms = (lc[m0], mc[m0])
        if self.ndim == 2:
            self.clm[lc, mc] /= temp.clm[lc, mc]
            self.slm[ls, ms] /= temp.slm[ls, ms]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[lc, mc, i] /= temp.clm[lc, mc]
                self.slm[ls, ms, i] /= temp.slm[ls, ms]
        else:
            old_month = self.month
            exclude1 = set(self.month) - set(temp.month)

            self.month = np.array(list(sorted(set(self.month) - exclude1)))
            self.time = np.array([self.time[i] for i in range(len(self.time)) if not (old_month[i] in exclude1)])

            for i in range(len(old_month)):
                for j in range(len(temp.month)):
                    if old_month[i] == temp.month[j]:
                        self.clm[:l1, :m1, i] /= temp.clm[:l1, :m1, j]
                        self.slm[:l1, :m1, i] /= temp.slm[:l1, :m1, j]

            to_keep = []
            for i in range(len(old_month)):
                if not (old_month[i] in exclude1):
                    to_keep.append(i)
            self.clm = self.clm[:, :, to_keep]
            self.slm = self.slm[:, :, to_keep]
        return self

    def copy(self):
        """
        Copy a ``harmonics`` object to a new ``harmonics`` object
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        # copy attributes or update attributes dictionary
        if isinstance(self.attributes, list):
            setattr(temp,'attributes',self.attributes)
        elif isinstance(self.attributes, dict):
            temp.attributes.update(self.attributes)
        # try to assign variables to self
        for key in ['clm','slm','time','month','filename']:
            try:
                val = getattr(self, key)
                setattr(temp, key, copy.copy(val))
            except AttributeError:
                pass
        # assign degree and order fields
        temp.update_dimensions()
        return temp

    def zeros(self, lmax=None, mmax=None, nt=None):
        """
        Create an empty ``harmonics`` object
        """
        # assign maximum degree and order
        if lmax is not None:
            self.lmax = np.int64(lmax)
        if mmax is not None:
            self.mmax = np.int64(mmax)
        # assign variables to self
        if nt is not None:
            self.clm = np.zeros((self.lmax+1, self.mmax+1, nt))
            self.slm = np.zeros((self.lmax+1, self.mmax+1, nt))
            self.time = np.zeros((nt))
            self.month = np.zeros((nt), dtype=int)
        else:
            self.clm = np.zeros((self.lmax+1, self.mmax+1))
            self.slm = np.zeros((self.lmax+1, self.mmax+1))
        # assign degree and order fields
        self.update_dimensions()
        return self

    def zeros_like(self):
        """
        Create a ``harmonics`` object using the dimensions of another
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        # assign variables to temp
        for key in ['clm','slm','time','month']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.zeros_like(val))
            except AttributeError:
                pass
        # assign degree and order fields
        temp.update_dimensions()
        return temp

    def expand_dims(self, update_dimensions=True):
        """
        Add a singleton dimension to a ``harmonics`` object if non-existent

        Parameters
        ----------
        update_dimensions: bool, default True
            Update the degree and order dimensions
        """
        # change time dimensions to be iterable
        self.time = np.atleast_1d(self.time)
        self.month = np.atleast_1d(self.month)
        # output harmonics with a third dimension
        if (self.ndim == 2) and not self.flattened:
            self.clm = self.clm[:, :, None]
            self.slm = self.slm[:, :, None]
        elif (self.ndim == 1) and self.flattened:
            self.clm = self.clm[:, None]
            self.slm = self.slm[:, None]
        # assign degree and order fields
        if update_dimensions:
            self.update_dimensions()
        # return the expanded harmonics object
        return self

    def squeeze(self, update_dimensions=True):
        """
        Remove singleton dimensions from a ``harmonics`` object

        Parameters
        ----------
        update_dimensions: bool, default True
            Update the degree and order dimensions
        """
        # squeeze singleton dimensions
        try:
            self.clm = np.squeeze(self.clm, axis=-1)
            self.slm = np.squeeze(self.slm, axis=-1)
            self.time = np.squeeze(self.time)
            self.month = np.squeeze(self.month)
        except ValueError as exc:
            pass
        # assign degree and order fields
        if update_dimensions:
            self.update_dimensions()
        else:
            self.ndim = self.clm.ndim
            self.shape = self.clm.shape
        return self

    def flatten(self, date=True):
        """
        Flatten ``harmonics`` matrices into arrays

        Parameters
        ----------
        date: bool, default True
            ``harmonics`` objects contain date information
        """
        n_harm = (self.lmax**2 + 3*self.lmax - (self.lmax-self.mmax)**2 -
            (self.lmax - self.mmax))//2 + 1
        # restructured degree and order
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.l = np.zeros((n_harm,), dtype=np.int64)
        temp.m = np.zeros((n_harm,), dtype=np.int64)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        # copy date variables if applicable
        if date:
            temp.time = np.copy(self.time)
            temp.month = np.copy(self.month)
        # restructured spherical harmonic arrays
        if (self.clm.ndim == 2):
            temp.clm = np.zeros((n_harm))
            temp.slm = np.zeros((n_harm))
        else:
            n = self.clm.shape[-1]
            temp.clm = np.zeros((n_harm, n))
            temp.slm = np.zeros((n_harm, n))
        # create counter variable lm
        lm = 0
        for m in range(0,self.mmax + 1):# MMAX+1 to include MMAX
            for l in range(m,self.lmax + 1):# LMAX+1 to include LMAX
                temp.l[lm] = np.int64(l)
                temp.m[lm] = np.int64(m)
                if self.clm.ndim == 2:
                    temp.clm[lm] = self.clm[l, m]
                    temp.slm[lm] = self.slm[l, m]
                else:
                    temp.clm[lm, :] = self.clm[l, m, :]
                    temp.slm[lm, :] = self.slm[l, m, :]
                # add 1 to lm counter variable
                lm += 1
        # update flattened attribute
        temp.flattened = True
        # return the flattened arrays
        return temp

    def expand(self, date=True):
        """
        Expand flattened ``harmonics`` into matrices

        Parameters
        ----------
        date: bool, default True
            ``harmonics`` objects contain date information
        """
        # number of harmonics
        n_harm = len(self.l)
        # restructured degree and order
        #n_harm = (self.lmax**2 + 3*self.lmax - (self.lmax - self.mmax)**2 -
        #    (self.lmax - self.mmax))//2 + 1
        # restructured degree and order
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        # copy date variables if applicable
        if date:
            temp.time = np.copy(self.time)
            temp.month = np.copy(self.month)
        # restructured spherical harmonic matrices
        if self.clm.ndim == 1:
            temp.clm = np.zeros((self.lmax + 1, self.mmax + 1))
            temp.slm = np.zeros((self.lmax + 1, self.mmax + 1))
        else:
            n = self.clm.shape[-1]
            temp.clm = np.zeros((self.lmax + 1,self.mmax + 1, n))
            temp.slm = np.zeros((self.lmax + 1,self.mmax + 1, n))
        # create counter variable lm
        for lm in range(n_harm):
            l = self.l[lm]
            m = self.m[lm]
            if self.clm.ndim == 1:
                temp.clm[l, m] = self.clm[lm]
                temp.slm[l, m] = self.slm[lm]
            else:
                temp.clm[l, m, :] = self.clm[lm, :]
                temp.slm[l, m, :] = self.slm[lm, :]
        # update flattened attribute
        temp.flattened = False
        # assign degree and order fields
        temp.update_dimensions()
        # return the expanded harmonics object
        return temp

    def index(self, indice, date=True):
        """
        Subset a ``harmonics`` object to specific index

        Parameters
        ----------
        indice: int
            index in matrix for subsetting
        date: bool, default True
            ``harmonics`` objects contain date information
        """
        # output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax), mmax=np.copy(self.mmax))
        # subset output harmonics
        temp.clm = self.clm[:, :, indice].copy()
        temp.slm = self.slm[:, :, indice].copy()
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
        # assign degree and order fields
        temp.update_dimensions()
        # return the subsetted object
        return temp

    def subset(self, months):
        """
        Subset a ``harmonics`` object to specific GRACE/GRACE-FO months

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
        months_list = [i for i, m in enumerate(self.month) if m in months]
        # output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax), mmax=np.copy(self.mmax))
        # create output harmonics
        temp.clm = np.zeros((temp.lmax + 1, temp.mmax + 1, n))
        temp.slm = np.zeros((temp.lmax + 1, temp.mmax + 1, n))
        temp.time = np.zeros((n))
        temp.month = np.zeros((n), dtype=np.int64)
        temp.filename = []
        # for each indice
        for t,i in enumerate(months_list):
            temp.clm[:,:, t] = self.clm[:,:, i].copy()
            temp.slm[:,:, t] = self.slm[:,:, i].copy()
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
            # subset filenames if applicable
            if getattr(self, 'filename'):
                if isinstance(self.filename, list):
                    temp.filename.append(str(self.filename[i]))
                elif isinstance(self.filename, str):
                    temp.filename.append(self.filename)
        # assign degree and order fields
        temp.update_dimensions()
        # remove singleton dimensions if importing a single value
        return temp.squeeze()

    def truncate(self, lmax, lmin=0, mmax=None):
        """
        Truncate or expand a ``harmonics`` object to a new degree and order

        Parameters
        ----------
        lmax: int
            maximum degree of spherical harmonics
        lmin: int, default 0
            minimum degree of spherical harmonics
        mmax: int or NoneType, default None
            maximum order of spherical harmonics
        """
        # output harmonics dimensions
        lmax = np.copy(self.lmax) if (lmax is None) else lmax
        mmax = np.copy(lmax) if (mmax is None) else mmax
        # copy prior harmonics object
        temp = self.copy()
        # set new degree and order
        self.lmax = np.copy(lmax)
        self.mmax = np.copy(mmax) if mmax else np.copy(lmax)
        # truncation levels
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        # create output harmonics
        if temp.ndim == 3:
            # number of months
            n = temp.clm.shape[-1]
            self.clm = np.zeros((self.lmax+1, self.mmax+1, n))
            self.slm = np.zeros((self.lmax+1, self.mmax+1, n))
            self.clm[lmin:l1, :m1,:] = temp.clm[lmin:l1, :m1,:].copy()
            self.slm[lmin:l1, :m1,:] = temp.slm[lmin:l1, :m1,:].copy()
        else:
            self.clm = np.zeros((self.lmax + 1, self.mmax + 1))
            self.slm = np.zeros((self.lmax + 1, self.mmax + 1))
            self.clm[lmin:l1, :m1] = temp.clm[lmin:l1, :m1].copy()
            self.slm[lmin:l1, :m1] = temp.slm[lmin:l1, :m1].copy()
        # assign degree and order fields
        self.update_dimensions()
        # return the truncated or expanded harmonics object
        return self

    def mean(self, apply=False, indices=Ellipsis):
        """
        Compute mean gravitational field and remove from data if specified

        Parameters
        ----------
        apply: bool, default False
            remove the mean field from the input ``harmonics`` object
        indices: int, default Ellipsis
            indices of input ``harmonics`` object to compute mean
        """
        temp = harmonics(lmax=np.copy(self.lmax), mmax=np.copy(self.mmax))
        # allocate for mean field
        temp.clm = np.zeros((temp.lmax + 1, temp.mmax + 1))
        temp.slm = np.zeros((temp.lmax + 1, temp.mmax + 1))
        # Computes the mean for each spherical harmonic degree and order
        for m in range(0, temp.mmax + 1):# MMAX+1 to include l
            for l in range(m, temp.lmax + 1):# LMAX+1 to include LMAX
                # calculate mean static field
                temp.clm[l, m] = np.mean(self.clm[l, m, indices])
                temp.slm[l, m] = np.mean(self.slm[l, m, indices])
                # calculating the time-variable gravity field by removing
                # the static component of the gravitational field
                if apply:
                    self.clm[l, m, :] -= temp.clm[l, m]
                    self.slm[l, m, :] -= temp.slm[l, m]
        # calculate mean of temporal variables
        for key in ['time','month']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.mean(val[indices]))
            except:
                continue
        # assign degree and order fields
        temp.update_dimensions()
        # return the mean field
        return temp

    def scale(self, var):
        """
        Multiply a ``harmonics`` object by a constant

        Parameters
        ----------
        var: float or np.ndarray
            scalar value to which the ``harmonics`` object will be multiplied
        """
        # assign degree and order fields
        self.update_dimensions()
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        # multiply by a single constant or a time-variable scalar
        if np.ndim(var) == 0:
            temp.clm = var*self.clm
            temp.slm = var*self.slm
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            temp.clm = np.zeros((temp.lmax + 1, temp.mmax + 1, len(var)))
            temp.slm = np.zeros((temp.lmax + 1, temp.mmax + 1, len(var)))
            for i,v in enumerate(var):
                temp.clm[:, :, i] = v*self.clm
                temp.slm[:, :, i] = v*self.slm
        elif (np.ndim(var) == 1) and (self.ndim == 3):
            for i,v in enumerate(var):
                temp.clm[:, :, i] = v*self.clm[:, :, i]
                temp.slm[:, :, i] = v*self.slm[:, :, i]
        # assign degree and order fields
        temp.update_dimensions()
        return temp

    def power(self, power):
        """
        Raise a ``harmonics`` object to a power

        Parameters
        ----------
        var: float
            power to which the ``harmonics`` object will be raised
        """
        # assign degree and order fields
        self.update_dimensions()
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        for key in ['clm','slm']:
            val = getattr(self, key)
            setattr(temp, key, np.power(val,power))
        # assign degree and order fields
        temp.update_dimensions()
        return temp

    def drift(self, t, epoch=2003.3):
        """
        Integrate a ``harmonics`` rate field over time to calculate drift

        Parameters
        ----------
        t: np.ndarray
            times for calculating drift
        epoch: float
            reference epoch for times
        """
        # assign degree and order fields
        self.update_dimensions()
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        # allocate for drift field
        temp.clm = np.zeros((temp.lmax+1,temp.mmax+1,len(t)))
        temp.slm = np.zeros((temp.lmax+1,temp.mmax+1,len(t)))
        # copy time variables and calculate GRACE/GRACE-FO months
        temp.time = np.copy(t)
        temp.month = calendar_to_grace(temp.time)
        # adjust months to fix special cases if necessary
        temp.month = adjust_months(temp.month)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        # calculate drift
        for i,ti in enumerate(t):
            temp.clm[:,:,i] = self.clm*(ti - epoch)
            temp.slm[:,:,i] = self.slm*(ti - epoch)
        # assign degree and order fields
        temp.update_dimensions()
        return temp

    def convolve(self, var):
        """
        Convolve ``harmonics`` with a degree-dependent array

        Parameters
        ----------
        var: np.ndarray
            degree dependent array for convolution
        """
        # assign degree and order fields
        self.update_dimensions()
        # check if a single field or a temporal field
        if self.ndim == 2:
            for l in range(0, self.lmax + 1):# LMAX+1 to include LMAX
                self.clm[l, :] *= var[l]
                self.slm[l, :] *= var[l]
        else:
            for i,t in enumerate(self.time):
                for l in range(0, self.lmax + 1):# LMAX+1 to include LMAX
                    self.clm[l, :, i] *= var[l]
                    self.slm[l, :, i] *= var[l]
        # return the convolved field
        return self

    def destripe(self, **kwargs):
        """
        Filters spherical harmonic coefficients for correlated "striping"
        errors following [Swenson2006]_

        Parameters
        ----------
        kwargs: dict
            keyword arguments for ``destripe_harmonics``

        References
        ----------
        .. [Swenson2006] S. Swenson and J. Wahr,
            "Post-processing removal of correlated errors in GRACE data",
            *Geophysical Research Letters*, 33(L08402), (2006).
            `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_
        """
        # assign degree and order fields
        self.update_dimensions()
        temp = harmonics(lmax=np.copy(self.lmax), mmax=np.copy(self.mmax))
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        # get filenames if applicable
        if getattr(self, 'filename'):
            temp.filename = copy.copy(self.filename)
        # check if a single field or a temporal field
        if self.ndim == 2:
            Ylms = destripe_harmonics(self.clm, self.slm,
                LMIN=1, LMAX=self.lmax, MMAX=self.mmax, **kwargs)
            temp.clm = Ylms['clm'].copy()
            temp.slm = Ylms['slm'].copy()
        else:
            n = self.shape[-1]
            temp.clm = np.zeros((self.lmax+1, self.mmax+1, n))
            temp.slm = np.zeros((self.lmax+1, self.mmax+1, n))
            for i in range(n):
                Ylms = destripe_harmonics(self.clm[:, :, i], self.slm[:, :, i],
                    LMIN=1, LMAX=self.lmax, MMAX=self.mmax, **kwargs)
                temp.clm[:, :, i] = Ylms['clm'].copy()
                temp.slm[:, :, i] = Ylms['slm'].copy()
        # assign degree and order fields
        temp.update_dimensions()
        # return the destriped field
        return temp

    @reify
    def amplitude(self):
        """
        Degree amplitude of the spherical harmonics
        """
        # temporary matrix for squared harmonics
        temp = self.power(2)
        # check if a single field or a temporal field
        if self.ndim == 2:
            # allocate for degree amplitudes
            amp = np.zeros((self.lmax + 1))
            for l in range(self.lmax + 1):
                # truncate at mmax
                m = np.arange(0, temp.mmax + 1)
                # degree amplitude of spherical harmonic degree
                amp[l] = np.sqrt(np.sum(temp.clm[l, m] + temp.slm[l, m]))
        else:
            # allocate for degree amplitudes
            n = self.shape[-1]
            amp = np.zeros((self.lmax + 1, n))
            for l in range(self.lmax + 1):
                # truncate at mmax
                m = np.arange(0, temp.mmax + 1)
                # degree amplitude of spherical harmonic degree
                var = temp.clm[l, m, :] + temp.slm[l, m, :]
                amp[l, :] = np.sqrt(np.sum(var, axis=0))
        # return the degree amplitudes
        return amp

    @property
    def dtype(self):
        """Main data type of ``harmonics`` object"""
        return self.clm.dtype

    @property
    def shape(self):
        """Dimensions of ``harmonics`` object
        """
        return np.shape(self.clm)

    @property
    def ndim(self):
        """Number of dimensions in ``harmonics`` object
        """
        return np.ndim(self.clm)

    @reify
    def ilm(self):
        """
        Complex form of the spherical harmonics
        """
        return self.clm - self.slm*1j

    def __len__(self):
        """Number of months
        """
        return len(self.month)

    def __iter__(self):
        """Iterate over GRACE/GRACE-FO months
        """
        self.__index__ = 0
        return self


    def __next__(self):
        """Get the next month of data
        """
        temp = harmonics(lmax=np.copy(self.lmax), mmax=np.copy(self.mmax))
        try:
            temp.time = self.time[self.__index__].copy()
            temp.month = self.month[self.__index__].copy()
            temp.clm = self.clm[:,:,self.__index__].copy()
            temp.slm = self.slm[:,:,self.__index__].copy()
        except IndexError as exc:
            raise StopIteration from exc
        # subset filename if applicable
        if getattr(self, 'filename'):
            if isinstance(self.filename, (list, tuple, np.ndarray)):
                temp.filename = str(self.filename[self.__index__])
            elif isinstance(self.filename, str):
                temp.filename = copy.copy(self.filename)
        # add to index
        self.__index__ += 1
        return temp

    def gap_fill(self, apply=False, interpolate=1):
        """
        Fill the missing months with a linear interpolation, the interpolation is made on month number, it's imprecise
        Options:
            apply: apply to the object if True, else return a new instance
            interpolate: 0 = fill gap with 0, 1 = linear interpolation
        """
        temp = self.copy()
        missing_month = self.month[-1] - self.month[0] - len(self.month) + 1

        temp.clm = np.zeros((self.lmax + 1, self.mmax + 1, len(self.time) + missing_month))
        temp.slm = np.zeros((self.lmax + 1, self.mmax + 1, len(self.time) + missing_month))
        temp.time = np.zeros(len(self.time) + missing_month)
        temp.month = np.arange(self.month[0], self.month[-1] + 1)

        # initialize index and count variables
        index = 0
        cmp = 0
        for i in range(int(self.month[0]), int(self.month[-1]) + 1):
            if i in self.month: # if month in original object, copy time and data
                cmp_miss_mon = 0 # variable for following missing months
                temp.time[index] = self.time[index - cmp]
                temp.clm[:, :, index] = self.clm[:, :, index - cmp]
                temp.slm[:, :, index] = self.slm[:, :, index - cmp]
            else: # fill values with a linear interpolation
                cmp += 1
                cmp_miss_mon += 1
                # y(t) = (y2 - y1)/(x2 - x1)*t + y1
                temp.time[index] = (self.time[index - cmp + 1] - self.time[index - cmp]) / (
                            self.month[index - cmp + 1] - self.month[index - cmp]) * cmp_miss_mon + self.time[index - cmp]

                if interpolate == 1:
                    temp.clm[:, :, index] = (self.clm[:, :, index - cmp + 1] - self.clm[:, :, index - cmp]) / \
                                            (self.month[index - cmp + 1] - self.month[index - cmp]) * cmp_miss_mon \
                                            + self.clm[:, :, index - cmp]
                    temp.slm[:, :, index] = (self.slm[:, :, index - cmp + 1] - self.slm[:, :, index - cmp]) / \
                                            (self.month[index - cmp + 1] - self.month[index - cmp]) * cmp_miss_mon \
                                            + self.slm[:, :, index - cmp]

                elif interpolate == 0:
                    temp.clm[:, :, index] = 0
                    temp.clm[:, :, index] = 0

            index += 1

        # -- assign ndim and shape attributes
        temp.update_dimensions()

        if apply:
            self.clm = temp.clm
            self.slm = temp.slm
            self.time = temp.time
            self.month = temp.month

            self.update_dimensions()

        return temp

    def plot_correlation(self, l, m, save_path=False):
        """
        Plot correlation between spherical harmonic coefficients of the object
        Inputs:
            l first degree of spherical harmonics
            m second degree of spherical harmonics

        Options:
            save_path : if not False, give a path to save the figure
        """
        mat_c = np.zeros((self.lmax, self.lmax))
        if m:
            mat_s = np.zeros((self.lmax, self.lmax))
        for i in range(self.lmax):
            for j in range(i+1):
                mat_c[i, i - j] = abs(np.mean((self.clm[l, m]-np.mean(self.clm[l, m]))*(self.clm[i, j]-np.mean(self.clm[i, j])))/\
                                   np.sqrt(np.mean((self.clm[l, m]-np.mean(self.clm[l, m]))**2))/\
                                   np.sqrt(np.mean((self.clm[i, j]-np.mean(self.clm[i, j]))**2)))

                if j:
                    mat_c[i - j, i] = abs(np.mean((self.clm[l, m]-np.mean(self.clm[l, m]))*(self.slm[i, j]-np.mean(self.slm[i, j])))/\
                                       np.sqrt(np.mean((self.clm[l, m]-np.mean(self.clm[l, m]))**2))/\
                                       np.sqrt(np.mean((self.slm[i, j]-np.mean(self.slm[i, j]))**2)))

                if m:
                    mat_s[i, i - j] = abs(np.mean(
                        (self.slm[l, m] - np.mean(self.slm[l, m])) * (self.clm[i, j] - np.mean(self.clm[i, j]))) / \
                                       np.sqrt(np.mean((self.slm[l, m] - np.mean(self.slm[l, m]))**2)) / \
                                       np.sqrt(np.mean((self.clm[i, j] - np.mean(self.clm[i, j]))**2)))

                    if j:
                        mat_s[i - j, i] = abs(np.mean(
                            (self.slm[l, m] - np.mean(self.slm[l, m])) * (self.slm[i, j] - np.mean(self.slm[i, j]))) / \
                                           np.sqrt(np.mean((self.slm[l, m] - np.mean(self.slm[l, m]))**2)) / \
                                           np.sqrt(np.mean((self.slm[i, j] - np.mean(self.slm[i, j]))**2)))

        plt.figure()
        plt.matshow(mat_c)
        plt.colorbar()
        plt.title('Correlation of each spherical harmonics with $C_{' + str(l) + ',' + str(m)+ '}$')

        if save_path:
            if os.path.isdir(save_path):
                plt.savefig(os.path.join(save_path, 'C' + str(l) + str(m) + '_correlation.png'))
            else:
                plt.savefig(save_path[:-3] + 'c' + save_path[-3:])

        if m:
            plt.figure()
            plt.matshow(mat_s)
            plt.colorbar()
            plt.title('Correlation of each spherical harmonics with $S_{' + str(l) + ',' + str(m) + '}$')

            if save_path:
                if os.path.isdir(save_path):
                    plt.savefig(os.path.join(save_path, 'S' + str(l) + str(m) + '_correlation.png'))
                else:
                    plt.savefig(save_path[:-3] + 's' + save_path[-3:])
        plt.show()


    def plot_coefficient(self, l, m, dates=[], ylms=[], label=[''], color=[], save_path=False):
        """
        Plot Cl,m and Sl,m harmonic coefficients
        Inputs:
            l first degree of spherical harmonics
            m second degree of spherical harmonics
        Options:
            dates: list with limits of the xaxis in year
            ylms: list of Harmonics objects to plot with the instance
            label: list of label for each Harmonics objects with element 0 representing the current Harmonics object
            save_path : if not False, give a path to save the figure
        """
        #-- figure for Cl,m
        plt.figure()
        ax = plt.gca()
        plt.title("Normalized spherical harmonics coefficient $C_{" + str(l) + "," + str(m) + "}$")
        if len(ylms):
            if len(color):
                plt.plot(self.time, self.clm[l, m, :], label=label[0], color=color[0])
            else:
                plt.plot(self.time, self.clm[l, m, :], label=label[0])
        else:
            plt.plot(self.time, self.clm[l, m, :], label="$C_{" + str(l) + "," + str(m) + "}$")

        try:
            for i in range(len(ylms)):
                if len(color):
                    plt.plot(ylms[i].time, ylms[i].clm[l, m, :], label=label[i + 1], color=color[i + 1])
                else:
                    plt.plot(ylms[i].time, ylms[i].clm[l, m, :], label=label[i + 1])
        except IndexError:
            raise IndexError("The list of labels is incomplete for correct plotting")

        plt.xlabel("Time (year)")
        plt.legend()
        ax.yaxis.offsetText.set_horizontalalignment('right')
        if dates:
            plt.xlim(dates)
        plt.grid()

        if save_path:
            if os.path.isdir(save_path):
                plt.savefig(os.path.join(save_path, 'C' + str(l) + str(m) + '_coefficient.png'))
            else:
                plt.savefig(save_path[:-4] + 'c' + save_path[-4:])

        if m:
            #-- figure for Sl,m
            plt.figure()
            ax = plt.gca()
            plt.title("Normalized spherical harmonic coefficient $S_{" + str(l) + "," + str(m) + "}$")
            if len(ylms):
                if len(color):
                    plt.plot(self.time, self.slm[l, m, :], label=label[0], color=color[0])
                else:
                    plt.plot(self.time, self.slm[l, m, :], label=label[0])
            else:
                plt.plot(self.time, self.slm[l, m, :], label="$S_{" + str(l) + "," + str(m) + "}$")

            try:
                for i in range(len(ylms)):
                    if len(color):
                        plt.plot(ylms[i].time, ylms[i].slm[l, m, :], label=label[i + 1], color=color[i + 1])
                    else:
                        plt.plot(ylms[i].time, ylms[i].slm[l, m, :], label=label[i + 1])
            except IndexError:
                raise IndexError("The list of labels is incomplete for correct plotting")

            plt.xlabel("Time (year)")
            plt.legend()
            ax.yaxis.offsetText.set_horizontalalignment('right')
            if dates:
                plt.xlim(dates)
            plt.grid()

            if save_path:
                if os.path.isdir(save_path):
                    plt.savefig(os.path.join(save_path, 'S' + str(l) + str(m) + '_coefficient.png'))
                else:
                    plt.savefig(save_path[:-4] + 's' + save_path[-4:])

        plt.show()

    def plot_fft(self, l, m, save_path=False, fmax=6):
        """
        Plot Cl,m and Sl,m harmonic coefficients fast fourrier transform
        Inputs:
            l first degree of spherical harmonics
            m second degree of spherical harmonics

        Options:
            save_path : if not False, give a path to save the figure
            fmax : maximal frequency (default to 6 for period > 2 months)
        """
        #-- compute fft and create x monthly frequency
        N = len(self.time)
        cf = sc.fft.fft(self.clm[l, m, :])[0:N // 2]
        sf = sc.fft.fft(self.slm[l, m, :])[0:N // 2]
        xf = np.linspace(0.0, 12/2, N // 2)

        # -- figure for Cl,m and Sl,m
        plt.figure()
        plt.title("Fourier transform of the normalized spherical harmonic coefficients $C_{" + str(l) + "," + str(
            m) + "}$ et $S_{" + str(
            l) + "," + str(m) + "}$")
        plt.plot(xf[xf <= fmax], 2.0 / N * np.abs(cf[xf <= fmax]), label="$C_{" + str(l) + "," + str(m) + "}$")
        if m:
            plt.plot(xf[xf <= fmax], 2.0 / N * np.abs(sf[xf <= fmax]), label="$S_{" + str(l) + "," + str(m) + "}$")


        plt.xlabel("Frequency ($year^{-1}$)")
        plt.ylabel("Power")
        plt.grid()
        plt.legend()

        if save_path:
            if os.path.isdir(save_path):
                plt.savefig(os.path.join(save_path, 'CS' + str(l) + str(m) + '_fft.png'))
            else:
                plt.savefig(save_path)

        plt.show()

    def plot_wavelets(self, l, m, s0=0, j1=None, pad=1, lag1=0, plot_coi=True, mother='MORLET', param=-1, func_plot=np.abs, save_path=False):
        """
        Plot Cl,m and Sl,m wavelet analysis based on (Torrence and Compo, 1998)

        Inputs:
            l first degree of spherical harmonics
            m second degree of spherical harmonics

        Options:
            s0 : minimal period of the wavelets, should be higher than 2*dt
            pad : boolean for the zero padding of the series
            lag1 : caracteristic of the noise: 0 for a white noise (default), 0.72 for a red noise
            plot_coi : boolean to display the cone of interest in the figure
            mother : name of the wavelet, can be MORLET, DOG or PAUL
            param : param of the wavelet, -1 is the default value for each wavelet
            func_plot : funtion for reducing the wave, can be np.abs, np.angle, np.real or np.imag
            save_path : if not False, give a path to save the figure
        """
        # len of the data
        ndata = self.time.shape[0]
        # compute the mean time delta of the object
        dt = np.mean((self.time[1:] - self.time[:-1]))

        # resolution of the wavelet
        dj = 0.005

        if not s0:
            s0 = 4 * dt  # min scale of the wavelets

        # max resolution of the wavelet, fixed for GRACE
        if j1 is None:
            j1 = np.log2(11/s0)/dj

        siglvl = 0.95

        # compute wavelets analysis of Cl,m and Sl,m
        wavec = wv.wavelet(self.clm[l,m], dt, pad, dj, s0, j1, mother, param)[0]
        waves, period, scale, coi = wv.wavelet(self.slm[l,m], dt, pad, dj, s0, j1, mother, param)

        # compute significativity of the wavelets
        signifc = wv.wave_signif(self.clm[l,m], dt, scale, lag1=lag1, siglvl=siglvl, mother=mother, param=param)
        signifs = wv.wave_signif(self.slm[l,m], dt, scale, lag1=lag1, siglvl=siglvl, mother=mother, param=param)

        # compute wavelet significance test at a level of confidence siglvl%
        sig95c = np.abs(wavec**2) / [s * np.ones(ndata) for s in signifc]
        sig95s = np.abs(waves**2) / [s * np.ones(ndata) for s in signifs]

        # Wavelet spectrum for fft plot
        global_wsc = (np.sum(np.abs(wavec ** 2).conj().transpose(), axis=0) / ndata)
        global_wss = (np.sum(np.abs(waves ** 2).conj().transpose(), axis=0) / ndata)

        # compute fft of the signal
        fft_sigc = np.fft.fft(self.clm[l,m])
        sxxc = np.abs((fft_sigc * np.conj(fft_sigc)) / ndata)[int(np.ceil(ndata / 2)):]
        fft_sigs = np.fft.fft(self.slm[l, m])
        sxxs = np.abs((fft_sigs * np.conj(fft_sigs)) / ndata)[int(np.ceil(ndata / 2)):]

        # compute frequency
        f = -np.fft.fftfreq(ndata)[int(np.ceil(ndata / 2)):]

        # prepare yticks
        yticks = []
        for i in [0.5, 1, 2, 4, 6, 10, 15]:
            if np.min(period) <= i <= np.max(period):
                yticks.append(i)

        # create figure Cl,m
        fig = plt.figure(constrained_layout=True, figsize=(12, 6), dpi=200)
        spec = matplotlib.gridspec.GridSpec(ncols=2, nrows=1, wspace=0.02, width_ratios=[3, 1])
        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[1], sharey=ax0)
        axs = [ax0, ax1]
        plt.setp(axs[1].get_yticklabels(), visible=False)

        # plot wavelet
        im = axs[0].contourf(self.time, period, func_plot(wavec), 100)
        axs[0].contour(self.time, period, sig95c, levels=[1], linewidths=2)

        # plot cone of interest of the wavelet
        if plot_coi:
            axs[0].fill(np.concatenate((self.time[:1] - 0.0001, self.time, self.time[-1:] + 0.0001,
                                        self.time[-1:] + 0.0001, self.time[:1] - 0.0001, self.time[:1] - 0.0001)),
                        np.concatenate(([np.min(period)], coi, [np.min(period)], period[-1:], period[-1:],
                                        [np.min(period)])), 'r', alpha=0.2, hatch='/')
            axs[0].plot(self.time, coi, 'r--', lw=1.4)

        fig.colorbar(im, ax=axs[0], location='left')
        axs[0].invert_yaxis()
        axs[0].set_yscale('log', base=2)
        axs[0].set_ylabel('Period (year)')
        axs[0].set_yticks(yticks)
        axs[0].set_ylim(np.max(period), np.min(period))
        axs[0].set_xlabel('Time (year)')
        axs[0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axs[0].set_title('Wavelet Power Spectrum')

        # plot fft analysis at the right of the figure
        axs[1].plot(sxxc, 1 / f * dt, 'gray', label='Fourier spectrum')
        axs[1].plot(global_wsc, period, 'b', label='Wavelet spectrum')
        axs[1].plot(np.array(signifc), period, 'g--', label='95% confidence spectrum')
        axs[1].set_xlabel('Power')
        axs[1].set_title('Global Wavelet Spectrum')

        plt.legend(loc='upper right')

        if save_path:
            if os.path.isdir(save_path):
                plt.savefig(os.path.join(save_path, 'C' + str(l) + str(m) + '_wavelet.png'))
            else:
                plt.savefig(save_path[:-4] + 'c' + save_path[-4:])

        if m:
            # create figure Sl,m
            fig = plt.figure(constrained_layout=True, figsize=(12, 6), dpi=200)
            spec = matplotlib.gridspec.GridSpec(ncols=2, nrows=1, wspace=0.02, width_ratios=[3, 1])
            ax0 = fig.add_subplot(spec[0])
            ax1 = fig.add_subplot(spec[1], sharey=ax0)
            axs = [ax0, ax1]
            plt.setp(axs[1].get_yticklabels(), visible=False)

            # plot wavelet
            im = axs[0].contourf(self.time, period, np.abs(waves), 100)
            axs[0].contour(self.time, period, sig95s, levels=[1], linewidths=2)

            # plot cone of interest of the wavelet
            if plot_coi:
                axs[0].fill(np.concatenate((self.time[:1] - 0.0001, self.time, self.time[-1:] + 0.0001,
                                            self.time[-1:] + 0.0001, self.time[:1] - 0.0001, self.time[:1] - 0.0001)),
                        np.concatenate(([s0], coi, [s0], period[-1:], period[-1:], [s0])), 'r', alpha=0.2, hatch='/')
                axs[0].plot(self.time, coi, 'r--', lw=1.4)

            fig.colorbar(im, ax=axs[0], location='left')
            axs[0].invert_yaxis()
            axs[0].set_yscale('log', base=2)
            axs[0].set_ylabel('Period (year)')
            axs[0].set_ylim(np.max(period), np.min(period))
            axs[0].set_yticks(yticks)
            axs[0].set_xlabel('Time (year)')
            axs[0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            axs[0].set_title('Wavelet Power Spectrum')

            # plot fft analysis at the right of the figure
            axs[1].plot(sxxs, 1 / f * dt, 'gray', label='Fourier spectrum')
            axs[1].plot(global_wss, period, 'b', label='Wavelet spectrum')
            axs[1].plot(np.array(signifs) * np.var(self.clm[l, m]), period, 'g--', label='95% confidence spectrum')
            axs[1].set_xlabel('Power')
            axs[1].set_title('Global Wavelet Spectrum')

            plt.legend(loc='upper right')

            if save_path:
                if os.path.isdir(save_path):
                    plt.savefig(os.path.join(save_path, 'S' + str(l) + str(m) + '_wavelet.png'))
                else:
                    plt.savefig(save_path[:-4] + 's' + save_path[-4:])

        plt.show()
