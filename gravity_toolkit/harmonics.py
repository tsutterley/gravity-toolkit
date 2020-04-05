#!/usr/bin/env python
u"""
harmonics.py
Written by Tyler Sutterley (04/2020)

Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (http://www.h5py.org/)

PROGRAM DEPENDENCIES:
    ncdf_stokes.py: writes output spherical harmonic data to netcdf
    hdf5_stokes.py: writes output spherical harmonic data to HDF5
    ncdf_read_stokes.py: reads spherical harmonic data from netcdf
    hdf5_read_stokes.py: reads spherical harmonic data from HDF5
    read_ICGEM_harmonics.py: reads gravity model coefficients from GFZ ICGEM
    destripe_harmonics.py: filters spherical harmonics for correlated errors

UPDATE HISTORY:
    Updated 04/2020: added from_gfc to read gravity model coefficients from GFZ
        add to_ascii and iterate over temporal fields in convolve and destripe
        make date optional for harmonic read functions.  add more math functions
        add option to sort if reading from an index or merging a list
    Written 03/2020
"""
import os
import re
import numpy as np
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.hdf5_stokes import hdf5_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes
from gravity_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from gravity_toolkit.destripe_harmonics import destripe_harmonics

class harmonics(object):
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
        self.filename=None

    def from_ascii(self, filename, date=True):
        """
        Read a harmonics object from an ascii file
        Inputs: full path of input ascii file
        Options: ascii file contains date information
        """
        self.filename = filename
        #-- read input ascii file (.txt) and split lines
        with open(os.path.expanduser(filename),'r') as f:
            file_contents = f.read().splitlines()
        #-- compile regular expression operator for extracting numerical values
        #-- from input ascii files of spherical harmonics
        regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[EeD][+-]?\d+)?'
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
            l1,m1 = np.array([l1,m1],dtype=np.int)
            #-- convert fortran exponentials if applicable
            self.clm[ll,mm] = np.float(clm1.replace('D','E'))
            self.slm[ll,mm] = np.float(slm1.replace('D','E'))
        return self

    def from_netCDF4(self, filename, date=True):
        """
        Read a harmonics object from a netCDF4 file
        Inputs: full path of input netCDF4 file
        Options: netCDF4 file contains date information
        """
        self.filename = filename
        Ylms = ncdf_read_stokes(os.path.expanduser(filename),
            ATTRIBUTES=False, DATE=date)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        if date:
            self.time = Ylms['time'].copy()
            self.month = Ylms['month'].copy()
        return self

    def from_HDF5(self, filename, date=True):
        """
        Read a harmonics object from a HDF5 file
        Inputs: full path of input HDF5 file
        Options: HDF5 file contains date information
        """
        self.filename = filename
        Ylms = hdf5_read_stokes(os.path.expanduser(filename),
            ATTRIBUTES=False, DATE=date)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        if date:
            self.time = Ylms['time'].copy()
            self.month = Ylms['month'].copy()
        return self

    def from_gfc(self, filename):
        """
        Read a harmonics object from a gfc gravity model file from the GFZ ICGEM
        Inputs: full path of input gfc file
        """
        self.filename = filename
        Ylms = read_ICGEM_harmonics(os.path.expanduser(filename))
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.lmax = np.int(Ylms['max_degree'])
        self.mmax = np.int(Ylms['max_degree'])
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        self.GM = np.float(Ylms['earth_gravity_constant'])
        self.R = np.float(Ylms['radius'])
        self.tide = Ylms['tide_system']
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
        self.filename = filename
        #-- Read index file of input spherical harmonics
        with open(os.path.expanduser(filename),'r') as f:
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
            harmonics objects contains date information
            sort harmonics objects by date information
        """
        #-- number of harmonic objects in list
        n = len(object_list)
        #-- indices to sort data objects if harmonics list contain dates
        if date and sort:
            list_sort = np.argsort([d.time for d in object_list],axis=None)
            self.time = np.zeros((n))
            self.month = np.zeros((n),dtype=np.int)
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
        #-- for each indice
        for t,i in enumerate(list_sort):
            self.clm[:,:,t] = object_list[i].clm[:self.lmax+1,:self.mmax+1]
            self.slm[:,:,t] = object_list[i].slm[:self.lmax+1,:self.mmax+1]
            if date:
                self.time[t] = object_list[i].time[:].copy()
                self.month[t] = object_list[i].month[:].copy()
        #-- return the single harmonic object
        return self

    def from_dict(self, d):
        """
        Convert a dict object to a harmonics object
        Inputs: dictionary object to be converted
        """
        #-- find valid keys within dictionary
        k = [k for k in ['l','m','clm','slm','time','month'] if k in d.keys()]
        #-- assign variables to self
        for key in k:
            setattr(self, key, d[key].copy())
        #-- maximum degree and order
        self.lmax = np.max(d['l'])
        self.mmax = np.max(d['m'])
        return self

    def to_ascii(self, filename, date=True):
        """
        Write a harmonics object to ascii file
        Inputs: full path of output ascii file
        Options: harmonics objects contains date information
        """
        self.filename = filename
        #-- open the output file
        fid = open(os.path.expanduser(filename), 'w')
        if date:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e} {4:10.4f}'
        else:
            file_format = '{0:5d} {1:5d} {2:+21.12e} {3:+21.12e}'
        #-- write to file for each spherical harmonic degree and order
        for m in range(0, self.mmax+1):
            for l in range(m, self.lmax+1):
                args = (l, m, Ylms.clm[l,m], Ylms.slm[l,m], Ylms.time)
                print(file_format.format(*args), file=fid)
        #-- close the output file
        fid.close()

    def to_netCDF4(self, filename, date=True):
        """
        Write a harmonics object to netCDF4 file
        Inputs: full path of output netCDF4 file
        Options: harmonics objects contains date information
        """
        self.filename = filename
        ncdf_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=os.path.expanduser(filename), TIME_UNITS='years',
            TIME_LONGNAME='Date_in_Decimal_Years', VERBOSE=False, DATE=date)

    def to_HDF5(self, filename, date=True):
        """
        Write a harmonics object to HDF5 file
        Inputs: full path of output HDF5 file
        Options: harmonics objects contains date information
        """
        self.filename = filename
        hdf5_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=os.path.expanduser(filename), TIME_UNITS='years',
            TIME_LONGNAME='Date_in_Decimal_Years', VERBOSE=False, DATE=date)

    def add(self, temp):
        """
        Add two harmonics objects
        Inputs: harmonic object to be added
        """
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (np.ndim(self.clm) == 2):
            self.clm[:l1,:m1] += temp.clm[:l1,:m1]
            self.slm[:l1,:m1] += temp.slm[:l1,:m1]
        elif (np.ndim(self.clm) == 3) and (np.ndim(temp.clm) == 2):
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
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (np.ndim(self.clm) == 2):
            self.clm[:l1,:m1] -= temp.clm[:l1,:m1]
            self.slm[:l1,:m1] -= temp.slm[:l1,:m1]
        elif (np.ndim(self.clm) == 3) and (np.ndim(temp.clm) == 2):
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
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if (np.ndim(self.clm) == 2):
            self.clm[:l1,:m1] *= temp.clm[:l1,:m1]
            self.slm[:l1,:m1] *= temp.slm[:l1,:m1]
        elif (np.ndim(self.clm) == 3) and (np.ndim(temp.clm) == 2):
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
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        #-- indices for cosine spherical harmonics (including zonals)
        lc,mc = np.tril_indices(l1, m=m1)
        #-- indices for sine spherical harmonics (excluding zonals)
        m0 = np.nonzero(mc != 0)
        ls,ms = (lc[m0],mc[m0])
        if (np.ndim(self.clm) == 2):
            self.clm[lc,mc] /= temp.clm[lc,mc]
            self.slm[ls,ms] /= temp.slm[ls,ms]
        elif (np.ndim(self.clm) == 3) and (np.ndim(temp.clm) == 2):
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
        #-- assign variables to self
        for key in ['clm','slm','time','month']:
            val = getattr(self, key)
            setattr(temp, key, np.copy(val))
        return temp

    def expand_dims(self):
        """
        Add a singleton dimension to a harmonics object if non-existent
        """
        #-- change time dimensions to be iterable
        if (np.ndim(self.time) == 0):
            self.time = np.array([self.time])
            self.month = np.array([self.month])
        #-- output harmonics with a third dimension
        if (np.ndim(self.clm) == 2):
            self.clm = self.clm[:,:,None]
            self.slm = self.slm[:,:,None]
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
        return self

    def index(self, indice):
        """
        Subset a harmonics object to specific index
        Inputs: indice in matrix to subset
        """
        #-- output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        #-- subset output harmonics
        temp.clm = self.clm[:,:,indice].copy()
        temp.slm = self.slm[:,:,indice].copy()
        temp.time = self.time[indice].copy()
        temp.month = self.month[indice].copy()
        return temp

    def subset(self, months):
        """
        Subset a harmonics object to specific GRACE/GRACE-FO months
        Inputs: GRACE/GRACE-FO months
        """
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
        #-- for each indice
        for t,i in enumerate(months_list):
            temp.clm[:,:,t] = self.clm[:,:,i].copy()
            temp.slm[:,:,t] = self.slm[:,:,i].copy()
            temp.time[t] = self.time[i].copy()
            temp.month[t] = self.month[i].copy()
        return temp

    def truncate(self, lmax, lmin=0, mmax=None):
        """
        Truncate or expand a harmonics object to a new degree and order
        Inputs: lmax maximum degree of spherical harmonics
        Option: lmin minimum degree of spherical harmonics
            mmax maximum order of spherical harmonics
        """
        #-- output harmonics object
        mmax = np.copy(lmax) if (mmax == None) else mmax
        #-- copy prior harmonics object
        temp = self.copy()
        #-- set new degree and order
        self.lmax = np.copy(lmax)
        self.mmax = np.copy(mmax) if mmax else np.copy(lmax)
        #-- truncation levels
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        #-- create output harmonics
        if (np.ndim(temp.clm) == 3):
            #-- number of months
            n = len(temp.month)
            self.clm = np.zeros((self.lmax+1,self.mmax+1,n))
            self.slm = np.zeros((self.lmax+1,self.mmax+1,n))
            self.clm[lmin:l1,:m1,:] = temp.clm[lmin:l1,:m1,:].copy()
            self.slm[lmin:l1,:m1,:] = temp.slm[lmin:l1,:m1,:].copy()
        else:
            self.clm = np.zeros((self.lmax+1,self.mmax+1))
            self.slm = np.zeros((self.lmax+1,self.mmax+1))
            self.clm[lmin:l1,:m1] = temp.clm[lmin:l1,:m1].copy()
            self.slm[lmin:l1,:m1] = temp.slm[lmin:l1,:m1].copy()
        #-- return the truncated or expanded harmonics object
        return self

    def mean(self, apply=False):
        """
        Compute mean gravitational field and remove from monthly if specified
        Option: apply to remove the mean field from the input harmonics
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
        #-- return the mean field
        return temp

    def scale(self, var):
        """
        Multiply a harmonics object by a constant
        Inputs: scalar value to which the harmonics object will be multiplied
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        for key in ['clm','slm']:
            val = getattr(self, key)
            setattr(temp, key, val*var)
        return temp

    def power(self, pow):
        """
        Raise a harmonics object to a power
        Inputs: power to which the harmonics object will be raised
        """
        temp = harmonics(lmax=self.lmax, mmax=self.mmax)
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        for key in ['clm','slm']:
            val = getattr(self, key)
            setattr(temp, key, np.power(val,pow))
        return temp

    def convolve(self, var):
        """
        Convolve spherical harmonics with a degree-dependent array
        Inputs: degree dependent array for convolution
        """
        #-- check if a single field or a temporal field
        if (np.ndim(self.clm) == 2):
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

    def destripe(self):
        """
        Filters spherical harmonic coefficients for correlated "striping" errors
        """
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        #-- check if a single field or a temporal field
        if (np.ndim(self.clm) == 2):
            Ylms = destripe_harmonics(self.clm, self.slm,
                LMIN=1, LMAX=self.lmax, MMAX=self.mmax)
            temp.clm = Ylms['clm'].copy()
            temp.slm = Ylms['slm'].copy()
        else:
            nt = len(self.time)
            temp.clm = np.zeros((self.lmax+1,self.mmax+1,nt))
            temp.slm = np.zeros((self.lmax+1,self.mmax+1,nt))
            for i in range(nt):
                Ylms = destripe_harmonics(self.clm[:,:,i], self.slm[:,:,i],
                    LMIN=1, LMAX=self.lmax, MMAX=self.mmax)
                temp.clm[:,:,i] = Ylms['clm'].copy()
                temp.slm[:,:,i] = Ylms['slm'].copy()
        #-- return the destriped field
        return temp
