#!/usr/bin/env python
u"""
harmonics.py
Written by Tyler Sutterley (12/2020)

Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (https://www.h5py.org/)

PROGRAM DEPENDENCIES:
    ncdf_stokes.py: writes output spherical harmonic data to netcdf
    hdf5_stokes.py: writes output spherical harmonic data to HDF5
    ncdf_read_stokes.py: reads spherical harmonic data from netcdf
    hdf5_read_stokes.py: reads spherical harmonic data from HDF5
    read_ICGEM_harmonics.py: reads gravity model coefficients from GFZ ICGEM
    destripe_harmonics.py: filters spherical harmonics for correlated errors

UPDATE HISTORY:
    Updated 12/2020: added verbose option for gfc files
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
import os
import re
import gzip
import zipfile
import numpy as np
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.hdf5_stokes import hdf5_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes
from gravity_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from gravity_toolkit.destripe_harmonics import destripe_harmonics

class harmonics(object):
    """
    Data class for reading, writing and processing spherical harmonic data
    """
    np.seterr(invalid='ignore')
    def __init__(self, lmax=None, mmax=None):
        self.clm=None
        self.slm=None
        self.time=None
        self.month=None
        self.lmax=lmax
        self.mmax=mmax
        self.l=np.arange(self.lmax+1) if self.lmax else None
        self.m=np.arange(self.mmax+1) if self.mmax else None
        self.shape=None
        self.ndim=None
        self.filename=None

    def case_insensitive_filename(self,filename):
        """
        Searches a directory for a filename without case dependence
        """
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

    def from_ascii(self, filename, date=True, compression=None, verbose=False):
        """
        Read a harmonics object from an ascii file
        Inputs: full path of input ascii file
        Options:
            ascii file contains date information
            ascii file is compressed using gzip
            verbose output of file information
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
        else:
            #-- read input ascii file (.txt, .asc) and split lines
            with open(self.filename,'r') as f:
                file_contents = f.read().splitlines()
        #-- compile regular expression operator for extracting numerical values
        #-- from input ascii files of spherical harmonics
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)
        #-- find maximum degree and order of harmonics
        self.lmax = 0
        self.mmax = 0
        #-- for each line in the file
        for line in file_contents:
            if date:
                l1,m1,clm1,slm1,time = rx.findall(line)
            else:
                l1,m1,clm1,slm1 = rx.findall(line)
            #-- convert line degree and order to integers
            l1,m1 = np.array([l1,m1],dtype=np.int)
            self.lmax = np.copy(l1) if (l1 > self.lmax) else self.lmax
            self.mmax = np.copy(m1) if (m1 > self.mmax) else self.mmax
        #-- output spherical harmonics dimensions array
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        #-- output spherical harmonics data
        self.clm = np.zeros((self.lmax+1,self.mmax+1))
        self.slm = np.zeros((self.lmax+1,self.mmax+1))
        #-- if the ascii file contains date variables
        if date:
            self.time = np.float(time)
            self.month = np.int(12.0*(self.time - 2002.0)) + 1
        #-- extract harmonics and convert to matrix
        #-- for each line in the file
        for line in file_contents:
            if date:
                l1,m1,clm1,slm1,time = rx.findall(line)
            else:
                l1,m1,clm1,slm1 = rx.findall(line)
            #-- convert line degree and order to integers
            ll,mm = np.array([l1,m1],dtype=np.int)
            #-- convert fortran exponentials if applicable
            self.clm[ll,mm] = np.float(clm1.replace('D','E'))
            self.slm[ll,mm] = np.float(slm1.replace('D','E'))
        #-- assign shape and ndim attributes
        self.update_dimensions()
        return self

    def from_netCDF4(self, filename, date=True, compression=None, verbose=False):
        """
        Read a harmonics object from a netCDF4 file
        Inputs: full path of input netCDF4 file
        Options:
            netCDF4 file contains date information
            netCDF4 file is compressed using gzip or zip
            verbose output of file information
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- read data from netCDF4 file
        Ylms = ncdf_read_stokes(self.filename, ATTRIBUTES=False,
            DATE=date, COMPRESSION=compression, VERBOSE=verbose)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        if date:
            self.time = Ylms['time'].copy()
            self.month = Ylms['month'].copy()
        #-- assign shape and ndim attributes
        self.update_dimensions()
        return self

    def from_HDF5(self, filename, date=True, compression=None, verbose=False):
        """
        Read a harmonics object from a HDF5 file
        Inputs: full path of input HDF5 file
        Options:
            HDF5 file contains date information
            HDF5 file is compressed using gzip or zip
            verbose output of file information
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- read data from HDF5 file
        Ylms = hdf5_read_stokes(self.filename, ATTRIBUTES=False,
            DATE=date, COMPRESSION=compression, VERBOSE=verbose)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        if date:
            self.time = Ylms['time'].copy()
            self.month = Ylms['month'].copy()
        #-- assign shape and ndim attributes
        self.update_dimensions()
        return self

    def from_gfc(self, filename, verbose=False):
        """
        Read a harmonics object from a gfc gravity model file from the GFZ ICGEM
        Inputs: full path of input gfc file
        Options:
            verbose output of file information
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- read data from gfc file
        Ylms = read_ICGEM_harmonics(self.filename)
        #-- Output file information
        if verbose:
            print(self.filename)
            print(list(Ylms.keys()))
        #-- copy variables for static gravity model
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.lmax = np.int(Ylms['max_degree'])
        self.mmax = np.int(Ylms['max_degree'])
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        #-- geophysical parameters of gravity model
        self.GM = np.float(Ylms['earth_gravity_constant'])
        self.R = np.float(Ylms['radius'])
        self.tide = Ylms['tide_system']
        #-- assign shape and ndim attributes
        self.update_dimensions()
        return self

    def from_index(self, filename, format=None, date=True, sort=True):
        """
        Read a harmonics object from an index of ascii, netCDF4 or HDF5 files
        Inputs: full path of index file to be read into a harmonics object
        Options:
            format of files in index (ascii, netCDF4 or HDF5)
            ascii, netCDF4, or HDF5 contains date information
            sort harmonics objects by date information
        """
        #-- set filename
        self.case_insensitive_filename(filename)
        #-- Read index file of input spherical harmonics
        with open(self.filename,'r') as f:
            file_list = f.read().splitlines()
        #-- create a list of harmonic objects
        h = []
        #-- for each file in the index
        for i,f in enumerate(file_list):
            if (format == 'ascii'):
                #-- ascii (.txt)
                h.append(harmonics().from_ascii(os.path.expanduser(f),date=date))
            elif (format == 'netCDF4'):
                #-- netcdf (.nc)
                h.append(harmonics().from_netCDF4(os.path.expanduser(f),date=date))
            elif (format == 'HDF5'):
                #-- HDF5 (.H5)
                h.append(harmonics().from_HDF5(os.path.expanduser(f),date=date))
        #-- create a single harmonic object from the list
        return self.from_list(h,date=date,sort=sort)

    def from_list(self, object_list, date=True, sort=True):
        """
        Build a sorted harmonics object from a list of other harmonics objects
        Inputs: list of harmonics object to be merged
        Options:
            harmonics objects contain date information
            sort harmonics objects by date information
        """
        #-- number of harmonic objects in list
        n = len(object_list)
        #-- indices to sort data objects if harmonics list contain dates
        if date and sort:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
        else:
            list_sort = np.arange(n)
        #-- truncate to maximum degree and order
        self.lmax = np.min([d.lmax for d in object_list])
        self.mmax = np.min([d.mmax for d in object_list])
        #-- output degree and order
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        #-- create output harmonics
        self.clm = np.zeros((self.lmax+1,self.mmax+1,n))
        self.slm = np.zeros((self.lmax+1,self.mmax+1,n))
        #-- create list of files
        self.filename = []
        #-- output dates
        if date:
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int)
        #-- for each indice
        for t,i in enumerate(list_sort):
            self.clm[:,:,t] = object_list[i].clm[:self.lmax+1,:self.mmax+1]
            self.slm[:,:,t] = object_list[i].slm[:self.lmax+1,:self.mmax+1]
            if date:
                self.time[t] = object_list[i].time[:].copy()
                self.month[t] = object_list[i].month[:].copy()
            #-- append filename to list
            if getattr(object_list[i], 'filename'):
                self.filename.append(object_list[i].filename)
        #-- assign shape and ndim attributes
        self.update_dimensions()
        #-- return the single harmonic object
        return self

    def from_dict(self, d):
        """
        Convert a dict object to a harmonics object
        Inputs: dictionary object to be converted
        """
        #-- assign dictionary variables to self
        for key in ['l','m','clm','slm','time','month']:
            try:
                setattr(self, key, d[key].copy())
            except (AttributeError, KeyError):
                pass
        #-- maximum degree and order
        self.lmax = np.max(d['l'])
        self.mmax = np.max(d['m'])
        #-- assign shape and ndim attributes
        self.update_dimensions()
        return self

    def to_ascii(self, filename, date=True):
        """
        Write a harmonics object to ascii file
        Inputs: full path of output ascii file
        Options: harmonics objects contain date information
        """
        self.filename = os.path.expanduser(filename)
        #-- open the output file
        fid = open(self.filename, 'w')
        if date:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e} {4:10.4f}'
        else:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e}'
        #-- write to file for each spherical harmonic degree and order
        for m in range(0, self.mmax+1):
            for l in range(m, self.lmax+1):
                args = (l, m, self.clm[l,m], self.slm[l,m], self.time)
                print(file_format.format(*args), file=fid)
        #-- close the output file
        fid.close()

    def to_netCDF4(self, filename, date=True, **kwargs):
        """
        Write a harmonics object to netCDF4 file
        Inputs: full path of output netCDF4 file
        Options: harmonics objects contain date information
        **kwargs: keyword arguments for ncdf_stokes
        """
        self.filename = os.path.expanduser(filename)
        if 'TIME_UNITS' not in kwargs.keys():
            kwargs['TIME_UNITS'] = 'years'
        if 'TIME_LONGNAME' not in kwargs.keys():
            kwargs['TIME_LONGNAME'] = 'Date_in_Decimal_Years'
        ncdf_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=self.filename, DATE=date, **kwargs)

    def to_HDF5(self, filename, date=True, **kwargs):
        """
        Write a harmonics object to HDF5 file
        Inputs: full path of output HDF5 file
        Options: harmonics objects contain date information
        **kwargs: keyword arguments for hdf5_stokes
        """
        self.filename = os.path.expanduser(filename)
        if 'TIME_UNITS' not in kwargs.keys():
            kwargs['TIME_UNITS'] = 'years'
        if 'TIME_LONGNAME' not in kwargs.keys():
            kwargs['TIME_LONGNAME'] = 'Date_in_Decimal_Years'
        hdf5_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=self.filename, DATE=date, **kwargs)

    def update_dimensions(self):
        """
        Update the dimensions of the spatial object
        """
        self.ndim = self.clm.ndim
        self.shape = self.clm.shape
        return self

    def add(self, temp):
        """
        Add two harmonics objects
        Inputs: harmonic object to be added
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (self.ndim == 2):
            self.clm[:l1,:m1] += temp.clm[:l1,:m1]
            self.slm[:l1,:m1] += temp.slm[:l1,:m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1,:m1,i] += temp.clm[:l1,:m1]
                self.slm[:l1,:m1,i] += temp.slm[:l1,:m1]
        else:
            self.clm[:l1,:m1,:] += temp.clm[:l1,:m1,:]
            self.slm[:l1,:m1,:] += temp.slm[:l1,:m1,:]
        return self

    def subtract(self, temp):
        """
        Subtract one harmonics object from another
        Inputs: harmonic object to be subtracted
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (self.ndim == 2):
            self.clm[:l1,:m1] -= temp.clm[:l1,:m1]
            self.slm[:l1,:m1] -= temp.slm[:l1,:m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1,:m1,i] -= temp.clm[:l1,:m1]
                self.slm[:l1,:m1,i] -= temp.slm[:l1,:m1]
        else:
            self.clm[:l1,:m1,:] -= temp.clm[:l1,:m1,:]
            self.slm[:l1,:m1,:] -= temp.slm[:l1,:m1,:]
        return self

    def multiply(self, temp):
        """
        Multiply two harmonics objects
        Inputs: harmonic object to be multiplied
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (self.ndim == 2):
            self.clm[:l1,:m1] *= temp.clm[:l1,:m1]
            self.slm[:l1,:m1] *= temp.slm[:l1,:m1]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[:l1,:m1,i] *= temp.clm[:l1,:m1]
                self.slm[:l1,:m1,i] *= temp.slm[:l1,:m1]
        else:
            self.clm[:l1,:m1,:] *= temp.clm[:l1,:m1,:]
            self.slm[:l1,:m1,:] *= temp.slm[:l1,:m1,:]
        return self

    def divide(self, temp):
        """
        Divide one harmonics object from another
        Inputs: harmonic object to be divided
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp.update_dimensions()
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        #-- indices for cosine spherical harmonics (including zonals)
        lc,mc = np.tril_indices(l1, m=m1)
        #-- indices for sine spherical harmonics (excluding zonals)
        m0 = np.nonzero(mc != 0)
        ls,ms = (lc[m0],mc[m0])
        if (self.ndim == 2):
            self.clm[lc,mc] /= temp.clm[lc,mc]
            self.slm[ls,ms] /= temp.slm[ls,ms]
        elif (self.ndim == 3) and (temp.ndim == 2):
            for i,t in enumerate(self.time):
                self.clm[lc,mc,i] /= temp.clm[lc,mc]
                self.slm[ls,ms,i] /= temp.slm[ls,ms]
        else:
            self.clm[lc,mc,:] /= temp.clm[lc,mc,:]
            self.slm[ls,ms,:] /= temp.slm[ls,ms,:]
        return self

    def copy(self):
        """
        Copy a harmonics object to a new harmonics object
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        #-- try to assign variables to self
        for key in ['clm','slm','time','month','shape','ndim','filename']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.copy(val))
            except AttributeError:
                pass
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        return temp

    def zeros_like(self):
        """
        Create a harmonics object using the dimensions of another
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        #-- assign variables to self
        for key in ['clm','slm','time','month']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.zeros_like(val))
            except AttributeError:
                pass
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a harmonics object if non-existent
        """
        #-- change time dimensions to be iterable
        self.time = np.atleast_1d(self.time)
        self.month = np.atleast_1d(self.month)
        #-- output harmonics with a third dimension
        if (self.ndim == 2):
            self.clm = self.clm[:,:,None]
            self.slm = self.slm[:,:,None]
        #-- reassign ndim and shape attributes
        self.update_dimensions()
        return self

    def squeeze(self):
        """
        Remove singleton dimensions from a harmonics object
        """
        #-- squeeze singleton dimensions
        self.time = np.squeeze(self.time)
        self.month = np.squeeze(self.month)
        self.clm = np.squeeze(self.clm)
        self.slm = np.squeeze(self.slm)
        #-- reassign ndim and shape attributes
        self.update_dimensions()
        return self

    def flatten(self, date=True):
        """
        Flatten harmonics matrices into arrays
        Options: harmonics objects contain date information
        """
        n_harm = (self.lmax**2 + 3*self.lmax - (self.lmax-self.mmax)**2 -
            (self.lmax-self.mmax))//2 + 1
        #-- restructured degree and order
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.l = np.zeros((n_harm,), dtype=np.int32)
        temp.m = np.zeros((n_harm,), dtype=np.int32)
        #-- copy date variables if applicable
        if date:
            temp.time = np.copy(self.time)
            temp.month = np.copy(self.month)
        #-- restructured spherical harmonic arrays
        if (self.clm.ndim == 2):
            temp.clm = np.zeros((n_harm))
            temp.slm = np.zeros((n_harm))
        else:
            n = self.clm.shape[-1]
            temp.clm = np.zeros((n_harm,n))
            temp.slm = np.zeros((n_harm,n))
        #-- create counter variable lm
        lm = 0
        for m in range(0,self.mmax+1):#-- MMAX+1 to include MMAX
            for l in range(m,self.lmax+1):#-- LMAX+1 to include LMAX
                temp.l[lm] = np.int(l)
                temp.m[lm] = np.int(m)
                if (self.clm.ndim == 2):
                    temp.clm[lm] = self.clm[l,m]
                    temp.slm[lm] = self.slm[l,m]
                else:
                    temp.clm[lm,:] = self.clm[l,m,:]
                    temp.slm[lm,:] = self.slm[l,m,:]
                #-- add 1 to lm counter variable
                lm += 1
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- return the flattened arrays
        return temp

    def expand(self, date=True):
        """
        Expand flattened harmonics into matrices
        Options: harmonics objects contain date information
        """
        n_harm = (self.lmax**2 + 3*self.lmax - (self.lmax-self.mmax)**2 -
            (self.lmax-self.mmax))//2 + 1
        #-- restructured degree and order
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        #-- copy date variables if applicable
        if date:
            temp.time = np.copy(self.time)
            temp.month = np.copy(self.month)
        #-- restructured spherical harmonic matrices
        if (self.clm.ndim == 1):
            temp.clm = np.zeros((self.lmax+1,self.mmax+1))
            temp.slm = np.zeros((self.lmax+1,self.mmax+1))
        else:
            n = self.clm.shape[-1]
            temp.clm = np.zeros((self.lmax+1,self.mmax+1,n))
            temp.slm = np.zeros((self.lmax+1,self.mmax+1,n))
        #-- create counter variable lm
        for lm in range(n_harm):
            l = self.l[lm]
            m = self.m[lm]
            if (self.clm.ndim == 1):
                temp.clm[l,m] = self.clm[lm]
                temp.slm[l,m] = self.slm[lm]
            else:
                temp.clm[l,m,:] = self.clm[lm,:]
                temp.slm[l,m,:] = self.slm[lm,:]
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- return the expanded harmonics object
        return temp

    def index(self, indice, date=True):
        """
        Subset a harmonics object to specific index
        Inputs: indice in matrix to subset
        Options: harmonics objects contain date information
        """
        #-- output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        #-- subset output harmonics
        temp.clm = self.clm[:,:,indice].copy()
        temp.slm = self.slm[:,:,indice].copy()
        #-- subset output dates
        if date:
            temp.time = self.time[indice].copy()
            temp.month = self.month[indice].copy()
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- subset filenames
        if getattr(self, 'filename'):
            temp.filename = self.filename[indice]
        return temp

    def subset(self, months):
        """
        Subset a harmonics object to specific GRACE/GRACE-FO months
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
        #-- output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        #-- create output harmonics
        temp.clm = np.zeros((temp.lmax+1,temp.mmax+1,n))
        temp.slm = np.zeros((temp.lmax+1,temp.mmax+1,n))
        temp.time = np.zeros((n))
        temp.month = np.zeros((n),dtype=np.int)
        temp.filename = []
        #-- for each indice
        for t,i in enumerate(months_list):
            temp.clm[:,:,t] = self.clm[:,:,i].copy()
            temp.slm[:,:,t] = self.slm[:,:,i].copy()
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
            if getattr(self, 'filename'):
                temp.filename.append(self.filename[i])
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- remove singleton dimensions if importing a single value
        return temp.squeeze()

    def truncate(self, lmax, lmin=0, mmax=None):
        """
        Truncate or expand a harmonics object to a new degree and order
        Inputs: lmax maximum degree of spherical harmonics
        Options: lmin minimum degree of spherical harmonics
            mmax maximum order of spherical harmonics
        """
        #-- output harmonics object
        mmax = np.copy(lmax) if (mmax is None) else mmax
        #-- copy prior harmonics object
        temp = self.copy()
        #-- set new degree and order
        self.lmax = np.copy(lmax)
        self.mmax = np.copy(mmax) if mmax else np.copy(lmax)
        #-- truncation levels
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        #-- create output harmonics
        if (temp.ndim == 3):
            #-- number of months
            n = temp.clm.shape[-1]
            self.clm = np.zeros((self.lmax+1,self.mmax+1,n))
            self.slm = np.zeros((self.lmax+1,self.mmax+1,n))
            self.clm[lmin:l1,:m1,:] = temp.clm[lmin:l1,:m1,:].copy()
            self.slm[lmin:l1,:m1,:] = temp.slm[lmin:l1,:m1,:].copy()
        else:
            self.clm = np.zeros((self.lmax+1,self.mmax+1))
            self.slm = np.zeros((self.lmax+1,self.mmax+1))
            self.clm[lmin:l1,:m1] = temp.clm[lmin:l1,:m1].copy()
            self.slm[lmin:l1,:m1] = temp.slm[lmin:l1,:m1].copy()
        #-- reassign ndim and shape attributes
        self.update_dimensions()
        #-- return the truncated or expanded harmonics object
        return self

    def mean(self, apply=False):
        """
        Compute mean gravitational field and remove from data if specified
        Options: apply to remove the mean field from the input harmonics
        """
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        #-- allocate for mean field
        temp.clm = np.zeros((temp.lmax+1,temp.mmax+1))
        temp.slm = np.zeros((temp.lmax+1,temp.mmax+1))
        #-- Computes the mean for each spherical harmonic degree and order
        for m in range(0,temp.mmax+1):#-- MMAX+1 to include l
            for l in range(m,temp.lmax+1):#-- LMAX+1 to include LMAX
                #-- calculate mean static field
                temp.clm[l,m] = np.mean(self.clm[l,m,:])
                temp.slm[l,m] = np.mean(self.slm[l,m,:])
                #-- calculating the time-variable gravity field by removing
                #-- the static component of the gravitational field
                if apply:
                    self.clm[l,m,:] -= temp.clm[l,m]
                    self.slm[l,m,:] -= temp.slm[l,m]
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- return the mean field
        return temp

    def scale(self, var):
        """
        Multiply a harmonics object by a constant
        Inputs: scalar value to which the harmonics object will be multiplied
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        #-- multiply by a single constant or a time-variable scalar
        if (np.ndim(var) == 0):
            temp.clm = var*self.clm
            temp.slm = var*self.slm
        elif (np.ndim(var) == 1) and (self.ndim == 2):
            temp.clm = np.zeros((temp.lmax+1,temp.mmax+1,len(var)))
            temp.slm = np.zeros((temp.lmax+1,temp.mmax+1,len(var)))
            for i,v in enumerate(var):
                temp.clm[:,:,i] = v*self.clm
                temp.slm[:,:,i] = v*self.slm
        elif (np.ndim(var) == 1) and (self.ndim == 3):
            for i,v in enumerate(var):
                temp.clm[:,:,i] = v*self.clm[:,:,i]
                temp.slm[:,:,i] = v*self.slm[:,:,i]
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        return temp

    def power(self, power):
        """
        Raise a harmonics object to a power
        Inputs: power to which the harmonics object will be raised
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        for key in ['clm','slm']:
            val = getattr(self, key)
            setattr(temp, key, np.power(val,power))
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        return temp

    def convolve(self, var):
        """
        Convolve spherical harmonics with a degree-dependent array
        Inputs: degree dependent array for convolution
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        #-- check if a single field or a temporal field
        if (self.ndim == 2):
            for l in range(0,self.lmax+1):#-- LMAX+1 to include LMAX
                self.clm[l,:] *= var[l]
                self.slm[l,:] *= var[l]
        else:
            for i,t in enumerate(self.time):
                for l in range(0,self.lmax+1):#-- LMAX+1 to include LMAX
                    self.clm[l,:,i] *= var[l]
                    self.slm[l,:,i] *= var[l]
        #-- return the convolved field
        return self

    def destripe(self, **kwargs):
        """
        Filters spherical harmonic coefficients for correlated "striping" errors
        Options: keyword arguments for destripe_harmonics
        """
        #-- reassign shape and ndim attributes
        self.update_dimensions()
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        #-- check if a single field or a temporal field
        if (self.ndim == 2):
            Ylms = destripe_harmonics(self.clm, self.slm,
                LMIN=1, LMAX=self.lmax, MMAX=self.mmax, **kwargs)
            temp.clm = Ylms['clm'].copy()
            temp.slm = Ylms['slm'].copy()
        else:
            n = self.shape[-1]
            temp.clm = np.zeros((self.lmax+1,self.mmax+1,n))
            temp.slm = np.zeros((self.lmax+1,self.mmax+1,n))
            for i in range(n):
                Ylms = destripe_harmonics(self.clm[:,:,i], self.slm[:,:,i],
                    LMIN=1, LMAX=self.lmax, MMAX=self.mmax, **kwargs)
                temp.clm[:,:,i] = Ylms['clm'].copy()
                temp.slm[:,:,i] = Ylms['slm'].copy()
        #-- assign ndim and shape attributes
        temp.update_dimensions()
        #-- return the destriped field
        return temp
