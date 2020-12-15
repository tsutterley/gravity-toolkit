#!/usr/bin/env python
u"""
mascon_reconstruct.py
Written by Tyler Sutterley (12/2020)

Calculates the equivalent spherical harmonics from a mascon time series

INPUTS:
    parameter files containing specific variables for each analysis

COMMAND LINE OPTIONS:
    --help: list the command line options
    -P X, --np X: run in parallel with X number of processes
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -r X, --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -M X, --mode X: permissions mode of the files created

PYTHON DEPENDENCIES:
    Python: a general-purpose, high-level programming language
        http://www.python.org/
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

PROGRAM DEPENDENCIES:
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using harmonics class for spherical harmonic operations
    Updated 04/2020: updates to reading load love numbers
        reading land-sea mask as a parameter
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: use a different land-sea mask for calculating ocean_Ylms
    Updated 03/2018: include a flag if using atmospheric ECMWF "jump" files
    Updated 02/2017: added parameter to redistribute mascon mass over the ocean
    Updated 06/2016: using __future__ print function
    Updated 02/2016: use getopt parameters to set number of PROCESSES to run
    Updated 06/2015: output harmonics for each mascon
    Updated 05/2015: added parameter MMAX for LMAX != MMAX.
        minor update to have parameters converted in function.
        added ascii and HDF5 output option
    Updated 01/2015: added error handling for multiprocessing threads
    Updated 10/2014: Distribute computing with multiprocessing module
    Updated 09/2014: Converted to function with main args
    Updated 05/2014
"""
from __future__ import print_function

import sys
import os
import argparse
import numpy as np
import multiprocessing
import traceback

from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.units import units
from gravity_toolkit.utilities import get_data_path

#-- PURPOSE: keep track of multiprocessing threads
def info(title):
    print(os.path.basename(sys.argv[0]))
    print(title)
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: read load love numbers for the range of spherical harmonic degrees
def load_love_numbers(LMAX, LOVE_NUMBERS=0, REFERENCE='CF'):
    """
    Reads PREM load Love numbers for the range of spherical harmonic degrees
    and applies isomorphic parameters

    Arguments
    ---------
    LMAX: maximum spherical harmonic degree

    Keyword arguments
    -----------------
    LOVE_NUMBERS: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth

    Returns
    -------
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement
    """
    #-- load love numbers file
    if (LOVE_NUMBERS == 0):
        #-- PREM outputs from Han and Wahr (1995)
        #-- https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
        love_numbers_file = get_data_path(['data','love_numbers'])
        header = 2
        columns = ['l','hl','kl','ll']
    elif (LOVE_NUMBERS == 1):
        #-- PREM outputs from Gegout (2005)
        #-- http://gemini.gsfc.nasa.gov/aplo/
        love_numbers_file = get_data_path(['data','Load_Love2_CE.dat'])
        header = 3
        columns = ['l','hl','ll','kl']
    elif (LOVE_NUMBERS == 2):
        #-- PREM outputs from Wang et al. (2012)
        #-- https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREM-LLNs-truncated.dat'])
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']    
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, LMAX=LMAX, HEADER=header,
        COLUMNS=columns, REFERENCE=REFERENCE, FORMAT='tuple')
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- PURPOSE: Reconstruct spherical harmonic fields from the mascon
#-- time series calculated in calc_mascon
def mascon_reconstruct(parameters,LOVE_NUMBERS=0,REFERENCE=None,MODE=0o775):
    #-- convert parameters into variables
    #-- Data processing center
    PROC = parameters['PROC']
    #-- Data Release
    DREL = parameters['DREL']
    #-- GRACE dataset
    DSET = parameters['DSET']
    #-- Date Range
    START_MON = np.int(parameters['START'])
    END_MON = np.int(parameters['END'])
    #-- spherical harmonic parameters
    LMAX = np.int(parameters['LMAX'])
    #-- maximum spherical harmonic order
    if (parameters['MMAX'].title() == 'None'):
        MMAX = np.copy(LMAX)
    else:
        MMAX = np.int(parameters['MMAX'])
    #-- gaussian smoothing radius
    RAD = np.float(parameters['RAD'])
    #-- filtered coefficients for stripe effects
    DESTRIPE = parameters['DESTRIPE'] in ('Y','y')
    #-- ECMWF jump corrections
    ATM = parameters['ATM'] in ('Y','y')
    #-- Glacial Isostatic Adjustment file to read
    GIA = parameters['GIA'] if (parameters['GIA'].title() != 'None') else None
    GIA_FILE = os.path.expanduser(parameters['GIA_FILE'])
    #-- input/output data format (ascii, netCDF4, HDF5)
    DATAFORM = parameters['DATAFORM']
    #-- index of mascons spherical harmonics
    #-- path.expanduser = tilde expansion of path
    MASCON_INDEX = os.path.expanduser(parameters['MASCON_INDEX'])
    #-- mascon distribution over the ocean
    MASCON_OCEAN = parameters['MASCON_OCEAN'] in ('Y','y')

    #-- for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else '_{0}'.format(DSET)
    #-- atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    #-- Gaussian smoothing string for radius RAD
    gw_str = '_r{0:0.0f}km'.format(RAD) if (RAD != 0) else ''
    #-- input GIA spherical harmonic datafiles
    GIA_Ylms_rate = read_GIA_model(GIA_FILE,GIA=GIA,LMAX=LMAX,MMAX=MMAX)
    gia_str = '_{0}'.format(GIA_Ylms_rate['title']) if GIA else ''
    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- filter grace coefficients flag
    ds_str = '_FL' if DESTRIPE else ''
    #-- output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- GRACE output filename prefix
    #-- mascon directory for GRACE product, processing and date range
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])

    #-- create initial reconstruct index for calc_mascon.py
    fid = open(os.path.expanduser(parameters['RECONSTRUCT_INDEX']),'w')
    HOME = os.path.expanduser('~')
    #-- output file format
    file_format = '{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}_{9:03d}-{10:03d}.{11}'

    #-- read load love numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)
    #-- Earth Parameters
    factors = units(lmax=LMAX).harmonic(hl,kl,ll)
    #-- Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    #-- Average Radius of the Earth [cm]
    rad_e = factors.rad_e
    #-- Read Ocean function and convert to Ylms for redistribution
    if MASCON_OCEAN:
        #-- read Land-Sea Mask and convert to spherical harmonics
        LANDMASK = os.path.expanduser(parameters['LANDMASK'])
        ocean_Ylms = ocean_stokes(LANDMASK,LMAX,MMAX=MMAX,LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        #-- not distributing uniformly over ocean
        ocean_str = ''

    #-- input mascon spherical harmonic datafiles
    with open(MASCON_INDEX,'r') as f:
        mascon_files = f.read().splitlines()
    for k,fi in enumerate(mascon_files):
        #-- read mascon spherical harmonics
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            Ylms=harmonics().from_ascii(os.path.expanduser(fi),date=False)
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            Ylms=harmonics().from_netCDF4(os.path.expanduser(fi),date=False)
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            Ylms=harmonics().from_HDF5(os.path.expanduser(fi),date=False)
        #-- Calculating the total mass of each mascon (1 cmwe uniform)
        total_area = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
        #-- distribute MASCON mass uniformly over the ocean
        if MASCON_OCEAN:
            #-- calculate ratio between total mascon mass and
            #-- a uniformly distributed cm of water over the ocean
            ratio = Ylms.clm[0,0]/ocean_Ylms['clm'][0,0]
            #-- for each spherical harmonic
            for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    #-- remove ratio*ocean Ylms from mascon Ylms
                    #-- note: x -= y is equivalent to x = x - y
                    Ylms.clm[l,m] -= ratio*ocean_Ylms['clm'][l,m]
                    Ylms.slm[l,m] -= ratio*ocean_Ylms['slm'][l,m]
        #-- truncate mascon spherical harmonics to d/o LMAX/MMAX
        Ylms = Ylms.truncate(lmax=LMAX, mmax=MMAX)
        #-- mascon base is the file without directory or suffix
        mascon_base = os.path.basename(fi)
        mascon_base = os.path.splitext(mascon_base)[0]
        #-- if lower case, will capitalize
        mascon_base = mascon_base.upper()
        #-- if mascon name contains degree and order info, remove
        mascon_name = mascon_base.replace('_L{0:d}'.format(LMAX),'')

        #-- input filename format (for both LMAX==MMAX and LMAX != MMAX cases):
        #-- mascon name, GRACE dataset, GIA model, LMAX, (MMAX,)
        #-- Gaussian smoothing, filter flag, remove reconstructed fields flag
        #-- output GRACE error file
        args = (mascon_name,dset_str,gia_str.upper(),atm_str,ocean_str,
            LMAX,order_str,gw_str,ds_str)
        file_input = '{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}.txt'.format(*args)
        mascon_data_input = np.loadtxt(os.path.join(DIRECTORY,file_input))
        nmon = np.shape(mascon_data_input)[0]

        #-- convert mascon time-series from Gt to cmwe
        mascon_sigma = 1e15*mascon_data_input[:,2]/total_area
        #-- mascon time-series Ylms
        mascon_Ylms = Ylms.scale(mascon_sigma)
        mascon_Ylms.time = mascon_data_input[:,1].copy()
        mascon_Ylms.month = mascon_data_input[:,0].astype(np.int)

        #-- output to file: no ascii option
        args = (mascon_name,dset_str,gia_str.upper(),atm_str,ocean_str,
            LMAX,order_str,gw_str,ds_str,START_MON,END_MON,suffix[DATAFORM])
        FILE = file_format.format(*args)
        #-- output harmonics to file
        if (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            mascon_Ylms.to_netCDF4(os.path.join(DIRECTORY,FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            mascon_Ylms.to_HDF5(os.path.join(DIRECTORY,FILE))
        #-- print file name to index
        print(os.path.join(DIRECTORY,FILE).replace(HOME,'~'), file=fid)
        #-- change the permissions mode
        os.chmod(os.path.join(DIRECTORY,FILE),MODE)
    #-- close the reconstruct index
    fid.close()
    #-- change the permissions mode of the index file
    os.chmod(os.path.expanduser(parameters['RECONSTRUCT_INDEX']),MODE)

#-- PURPOSE: define the analysis for multiprocessing
def define_analysis(parameter_file,LOVE_NUMBERS=0,REFERENCE=None,MODE=0o775):
    #-- keep track of multiprocessing threads
    info(os.path.basename(parameter_file))

    #-- variable with parameter definitions
    parameters = {}
    #-- Opening parameter file and assigning file ID number (fid)
    fid = open(os.path.expanduser(parameter_file), 'r')
    #-- for each line in the file will extract the parameter (name and value)
    for fileline in fid:
        #-- Splitting the input line between parameter name and value
        part = fileline.split()
        #-- filling the parameter definition variable
        parameters[part[0]] = part[1]
    #-- close the parameter file
    fid.close()

    #-- try to run the analysis with listed parameters
    try:
        #-- run the reconstruction function with chosen parameters
        mascon_reconstruct(parameters, LOVE_NUMBERS=LOVE_NUMBERS,
            REFERENCE=REFERENCE, MODE=MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser()
    #-- command line parameters
    parser.add_argument('parameters',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Parameter files containing specific variables for each analysis')
    #-- number of processes to run in parallel
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to run in parallel')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference','-r',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- use parameter files from system arguments listed after the program.
    if (args.np == 0):
        #-- run directly as series if PROCESSES = 0
        #-- for each entered parameter file
        for f in args.parameters:
            define_analysis(f,LOVE_NUMBERS=args.love,
                REFERENCE=args.reference,MODE=args.mode)
    else:
        #-- run in parallel with multiprocessing Pool
        pool = multiprocessing.Pool(processes=args.np)
        #-- for each entered parameter file
        for f in args.parameters:
            kwds=dict(LOVE_NUMBERS=args.love,REFERENCE=args.reference,
                MODE=args.mode)
            pool.apply_async(define_analysis,args=(f,),kwds=kwds)
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()


#-- run main program
if __name__ == '__main__':
    main()
