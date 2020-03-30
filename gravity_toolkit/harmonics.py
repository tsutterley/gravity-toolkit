#!/usr/bin/env python
u"""
harmonics.py
Written by Tyler Sutterley (03/2020)

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
    destripe_harmonics.py: filters spherical harmonics for correlated  errors

UPDATE HISTORY:
    Written 03/2020
"""
import os
import re
import numpy as np
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.hdf5_stokes import hdf5_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes
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

    def from_ascii(self, filename):
        """
        read a harmonics object from an ascii file
        """
        self.filename = filename
        Ylms = np.loadtxt(os.path.expanduser(input_file), dtype={'names':
            ('l','m','clm','slm','time'), 'formats':('i','i','f8','f8','f8')})
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        self.clm = np.zeros((self.lmax+1,self.mmax+1))
        self.slm = np.zeros((self.lmax+1,self.mmax+1))
        self.time = np.copy(Ylms['time'][0])
        self.month = np.int(12.0*(self.time - 2002.0)) + 1
        for ll,mm,clm,slm in zip(Ylms['l'],Ylms['m'],Ylms['clm'],Ylms['slm']):
            self.clm[ll,mm] = clm.copy()
            self.slm[ll,mm] = slm.copy()
        return self

    def from_netCDF4(self, filename):
        """
        read a harmonics object from a netCDF4 file
        """
        self.filename = filename
        Ylms = ncdf_read_stokes(os.path.expanduser(filename),
            ATTRIBUTES=False, DATE=True)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.time = Ylms['time'].copy()
        self.month = Ylms['month'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        return self

    def from_HDF5(self, filename):
        """
        read a harmonics object from a HDF5 file
        """
        self.filename = filename
        Ylms = hdf5_read_stokes(os.path.expanduser(filename),
            ATTRIBUTES=False, DATE=True)
        self.clm = Ylms['clm'].copy()
        self.slm = Ylms['slm'].copy()
        self.time = Ylms['time'].copy()
        self.month = Ylms['month'].copy()
        self.l = Ylms['l'].copy()
        self.m = Ylms['m'].copy()
        self.lmax = np.max(Ylms['l'])
        self.mmax = np.max(Ylms['m'])
        return self

    def from_index(self, filename, format=None):
        """
        read a harmonics object from an index of ascii, netCDF4 or HDF5 files
        """
        self.filename = filename
        #-- Read index file of input spherical harmonics
        with open(os.path.expanduser(filename),'r') as f:
            file_list = f.read().splitlines()
        #-- create a list of harmonic objects
        hlist = []
        #-- for each file in the index
        for i,f in enumerate(file_list):
            if (format == 'ascii'):
                #-- ascii (.txt)
                hlist.append(harmonics().from_ascii(os.path.expanduser(f)))
            elif (format == 'netCDF4'):
                #-- netcdf (.nc)
                hlist.append(harmonics().from_netCDF4(os.path.expanduser(f)))
            elif (format == 'HDF5'):
                #-- HDF5 (.H5)
                hlist.append(harmonics().from_HDF5(os.path.expanduser(f)))
        #-- create a single harmonic object from the list
        return self.from_list(hlist)

    def from_list(self, object_list):
        """
        build a sorted harmonics object from a list of other harmonics objects
        """
        #-- number of harmonic objects in list
        n = len(object_list)
        #-- indices to sort data objects
        list_sort = np.argsort([d.time for d in object_list],axis=None)
        #-- truncate to maximum degree and order
        self.lmax = np.min([d.lmax for d in object_list])
        self.mmax = np.min([d.mmax for d in object_list])
        #-- output degree and order
        self.l = np.arange(self.lmax+1)
        self.m = np.arange(self.mmax+1)
        #-- create output harmonics
        self.clm = np.zeros((self.lmax+1,self.mmax+1,n))
        self.slm = np.zeros((self.lmax+1,self.mmax+1,n))
        self.time = np.zeros((n))
        self.month = np.zeros((n),dtype=np.int)
        #-- for each indice
        for t,i in enumerate(list_sort):
            self.clm[:,:,t] = object_list[i].clm[:self.lmax+1,:self.mmax+1]
            self.slm[:,:,t] = object_list[i].slm[:self.lmax+1,:self.mmax+1]
            self.time[t] = object_list[i].time[:].copy()
            self.month[t] = object_list[i].month[:].copy()
        return self

    def from_dict(self, d):
        """
        convert a dict object to a harmonics object
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

    def add(self, temp):
        """
        add two harmonics objects
        """
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if np.ndim(self.clm == 2):
            self.clm[:l1,:m1] += temp.clm[:l1,:m1]
            self.slm[:l1,:m1] += temp.slm[:l1,:m1]
        else:
            self.clm[:l1,:m1,:] += temp.clm[:l1,:m1,:]
            self.slm[:l1,:m1,:] += temp.slm[:l1,:m1,:]
        return self

    def subtract(self, temp):
        """
        subtract one harmonics object from another
        """
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        if np.ndim(self.clm == 2):
            self.clm[:l1,:m1] -= temp.clm[:l1,:m1]
            self.slm[:l1,:m1] -= temp.slm[:l1,:m1]
        else:
            self.clm[:l1,:m1,:] -= temp.clm[:l1,:m1,:]
            self.slm[:l1,:m1,:] -= temp.slm[:l1,:m1,:]
        return self

    def index(self, indice):
        """
        subset a harmonics object to specific index
        """
        #-- output harmonics object
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        #-- subset output harmonics
        temp.clm = self.clm[:,:,indice]
        temp.slm = self.slm[:,:,indice]
        temp.time = self.time[indice].copy()
        temp.month = self.month[indice].copy()
        return temp

    def subset(self, months):
        """
        subset a harmonics object to specific GRACE/GRACE-FO months
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

    def truncate(self, lmax, mmax=None):
        """
        truncate a harmonics object to a new degree and order
        """
        #-- number of months
        n = len(self.month)
        #-- output harmonics object
        mmax = np.copy(lmax) if (mmax == None) else mmax
        temp = harmonics(lmax=lmax,mmax=mmax)
        #-- date variables
        temp.time = np.copy((self.time))
        temp.month = np.copy((self.month))
        #-- create output harmonics
        temp.clm = np.zeros((temp.lmax+1,temp.mmax+1,n))
        temp.slm = np.zeros((temp.lmax+1,temp.mmax+1,n))
        #-- truncation levels
        l1 = self.lmax+1 if (temp.lmax > self.lmax) else temp.lmax+1
        m1 = self.mmax+1 if (temp.mmax > self.mmax) else temp.mmax+1
        temp.clm[:l1,:m1,:] = self.clm[:l1,:m1,:].copy()
        temp.slm[:l1,:m1,:] = self.slm[:l1,:m1,:].copy()
        return temp

    def mean(self, apply=False):
        """
        Compute mean gravitational field and remove from monthly if specified
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

    def convolve(self, var):
        """
        Convolve spherical harmonics with a degree-dependent array
        """
        for l in range(0,self.lmax+1):#-- LMAX+1 to include LMAX
            self.clm[l,:] *= var[l]
            self.slm[l,:] *= var[l]
        #-- return the convolved field
        return self

    def destripe(self):
        """
        Filters spherical harmonic coefficients for correlated "striping" errors
        """
        temp = harmonics(lmax=np.copy(self.lmax),mmax=np.copy(self.mmax))
        Ylms = destripe_harmonics(self.clm, self.slm, LMIN=1,
            LMAX=self.lmax, MMAX=self.mmax)
        temp.clm = Ylms['clm'].copy()
        temp.slm = Ylms['slm'].copy()
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        return temp

    def units(self, hl, kl, ll):
        """
        Calculates degree dependent factors for converting units
        """
        # Earth Parameters
        # Average Density of the Earth [g/cm^3]
        rho_e = 5.517
        # Average Radius of the Earth [cm]
        rad_e = 6.371e8
        # WGS84 Gravitational Constant of the Earth [cm^3/s^2]
        GM_e = 3986004.418e14
        # Gravitational Constant of the Earth's atmosphere
        GM_atm = 3.5e14
        # Gravitational Constant of the Earth (w/o atm)
        GM = GM_e - GM_atm
        # standard gravitational acceleration (World Meteorological Organization)
        g_wmo = 980.665
        # degree dependent coefficients
        # norm, geodesy normalized spherical harmonics
        self.norm = np.ones((self.lmax+1))
        # cmwe, centimeters water equivalent
        self.cmwe = rho_e*rad_e*(2.0*self.l+1.0)/(1.0 +kl[self.l])/3.0
        # mmGH, millimeters geoid height
        self.mmGH = np.ones((self.lmax+1))*(10.0*rad_e)
        # mmCU, millimeters elastic crustal deformation
        self.mmCU = 10.0*rad_e*hl[self.l]/(1.0 +kl[self.l])
        # microGal, microGal gravity perturbations
        self.microGal = 1.e6*GM*(self.l+1.0)/(rad_e**2.0)
        # mbar, millibar equivalent surface pressure
        self.mbar = g_wmo*rho_e*rad_e*(2.0*self.l+1.0)/(1.0 +kl[self.l])/3e3
        # Pa, pascals equivalent surface pressure
        self.Pa = g_wmo*rho_e*rad_e*(2.0*self.l+1.0)/(1.0 +kl[self.l])/30.0
        # return the degree dependent unit conversions
        return self

    def to_netCDF4(self, filename):
        """
        write a harmonics object to netCDF4 file
        """
        self.filename = filename
        ncdf_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=os.path.expanduser(filename), TIME_UNITS='years',
            TIME_LONGNAME='Date_in_Decimal_Years', VERBOSE=False, DATE=True)

    def to_HDF5(self, filename):
        """
        write a harmonics object to HDF5 file
        """
        self.filename = filename
        hdf5_stokes(self.clm, self.slm, self.l, self.m, self.time, self.month,
            FILENAME=os.path.expanduser(filename), TIME_UNITS='years',
            TIME_LONGNAME='Date_in_Decimal_Years', VERBOSE=False, DATE=True)
