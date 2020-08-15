#!/usr/bin/env python
u"""
grace_spatial_maps.py
Written by Tyler Sutterley (08/2020)

Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports
    monthly spatial fields

Will correct with the specified GIA model group, destripe/smooth/process,
    and export the data in specified units

Spatial output units: cm w.e., mm geoid height, mm elastic uplift,
    microgal gravity perturbation or surface pressure (Pa)

CALLING SEQUENCE:
    python grace_spatial_maps.py input_parameters_file

    Can also input several parameter files in series:
    python grace_spatial_maps.py parameter_file1 parameter_file2

    Can be run in parallel with the python multiprocessing package:
    python grace_spatial_maps.py --np=2 parameter_file1 parameter_file2
    python grace_spatial_maps.py -P 2 parameter_file1 parameter_file2

    Can output a log file listing the input parameters and output files:
    python grace_spatial_maps.py --log parameter_file
    python grace_spatial_maps.py -l parameter_file

SYSTEM ARGUMENTS README:
    program is run as:
    python run_grace_process_input.py inp1 inp2 inp3
        where inp1, inp2 and inp3 are different inputs

        firstinput=sys.argv[1] (in this case inp1)
        secondinput=sys.argv[2] (in this case inp2)
        thirdinput=sys.argv[3] (in this case inp3)

    As python is base 0, sys.argv[0] is equal to run_grace_process_input.py
        (which is useful in some applications, but not for this program)

    For this program, the system arguments are parameter files
    The program reads the parameter file, which is separated by column as:
        Column 1: parameter name (such as LMAX)
        Column 2: parameter (e.g. 60)
        Column 3: comments (which are discarded)
    The parameters are stored in a python dictionary (variables indexed by keys)
        the keys are the parameter name (for LMAX: parameters['LMAX'] == 60)

INPUTS:
    parameter file (entered after program)
        parameter files contain specific variables for each analysis

COMMAND LINE OPTIONS:
    --help: list the command line options
    -P X, --np=X: Run in parallel with X number of processes
    -l, --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
    -M X, --mode=X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Replaces Degree 1 with with Swenson values (if specified)
        Replaces C20 and C30 with SLR values (if specified)
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    plm_holmes.py: Computes fully normalized associated Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: converts a land-sea mask to a series of spherical harmonics
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    geocenter.py: converts between spherical harmonics and geocenter variations
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    ncdf_read_stokes.py: reads spherical harmonic netcdf files
    ncdf_stokes.py: writes output spherical harmonic data to netcdf
    hdf5_read_stokes.py: reads spherical harmonic HDF5 files
    hdf5_stokes.py: writes output spherical harmonic data to HDF5
    spatial.py: spatial data class for reading, writing and processing data
    ncdf_read.py: reads input spatial data from netCDF4 files
    hdf5_read.py: reads input spatial data from HDF5 files
    ncdf_write.py: writes output spatial data to netCDF4
    hdf5_write.py: writes output spatial data to HDF5
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for output operations
    Updated 05/2020: for public release
"""
from __future__ import print_function

import sys
import os
import re
import time
import numpy as np
import getopt
import multiprocessing
import traceback

from gravity_toolkit.grace_input_months import grace_input_months
from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.spatial import spatial
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
def load_love_numbers(LMAX, REFERENCE='CF'):
    #-- load love numbers file
    love_numbers_file = get_data_path(['data','love_numbers'])
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    if (LMAX > 696):
        #-- Creates arrays of kl, hl, and ll Love Numbers
        hl = np.zeros((LMAX+1))
        kl = np.zeros((LMAX+1))
        ll = np.zeros((LMAX+1))
        hl[:697],kl[:697],ll[:697] = read_love_numbers(love_numbers_file,
            FORMAT='tuple', REFERENCE=REFERENCE)
        #-- for degrees greater than 696
        for l in range(697,LMAX+1):
            hl[l] = 2.0*hl[l-1] - hl[l-2]#-- linearly extrapolating hl
            kl[l] = 2.0*kl[l-1] - kl[l-2]#-- linearly extrapolating kl
            ll[l] = 2.0*ll[l-1] - ll[l-2]#-- linearly extrapolating ll
    else:
        #-- read arrays of kl, hl, and ll Love Numbers
        hl,kl,ll=read_love_numbers(love_numbers_file, REFERENCE=REFERENCE)
    #-- return a tuple of load love numbers
    return (hl,kl,ll)

#-- PURPOSE: import GRACE files for a given months range
#-- Converts the GRACE/GRACE-FO harmonics applying the specified procedures
def grace_spatial_maps(base_dir, parameters, VERBOSE, MODE):
    #-- Data processing center
    PROC = parameters['PROC']
    #-- Data Release
    DREL = parameters['DREL']
    #-- GRACE dataset
    DSET = parameters['DSET']
    #-- Date Range and missing months
    start_mon = np.int(parameters['START'])
    end_mon = np.int(parameters['END'])
    missing = np.array(parameters['MISSING'].split(','),dtype=np.int)
    #-- minimum degree
    LMIN = np.int(parameters['LMIN'])
    #-- maximum degree and order
    LMAX = np.int(parameters['LMAX'])
    if (parameters['MMAX'].title() == 'None'):
        MMAX = np.copy(LMAX)
    else:
        MMAX = np.int(parameters['MMAX'])
    #-- SLR C2,0 and C3,0
    SLR_C20 = parameters['SLR_C20']
    SLR_C30 = parameters['SLR_C30']
    #-- Degree 1 correction
    DEG1 = parameters['DEG1']
    MODEL_DEG1 = parameters['MODEL_DEG1'] in ('Y','y')
    #-- ECMWF jump corrections
    ATM = parameters['ATM'] in ('Y','y')
    #-- Pole Tide correction from Wahr et al. (2015)
    POLE_TIDE = parameters['POLE_TIDE'] in ('Y','y')
    #-- Glacial Isostatic Adjustment file to read
    GIA = parameters['GIA'] if (parameters['GIA'].title() != 'None') else None
    GIA_FILE = os.path.expanduser(parameters['GIA_FILE'])
    #-- remove a set of spherical harmonics from the GRACE data
    REMOVE_INDEX = parameters['REMOVE_INDEX']
    REDISTRIBUTE_REMOVED = parameters['REDISTRIBUTE_REMOVED'] in ('Y','y')
    #-- smoothing radius
    RAD = np.int(parameters['RAD'])
    #-- destriped coefficients
    DESTRIPE = parameters['DESTRIPE'] in ('Y','y')
    #-- output spatial units
    UNITS = np.int(parameters['UNITS'])
    #-- output degree spacing
    #-- can enter dlon and dlat as [dlon,dlat] or a single value
    DDEG = np.squeeze(np.array(parameters['DDEG'].split(','),dtype='f'))
    #-- output degree interval (0:360, 90:-90) or (degree spacing/2)
    INTERVAL = np.int(parameters['INTERVAL'])
    #-- output data format (1: ascii, 2: netcdf, 3: HDF5)
    DATAFORM = np.int(parameters['DATAFORM'])
    #-- output directory and base filename
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    FILENAME = parameters['FILENAME']

    #-- recursively create output directory if not currently existing
    if (not os.access(DIRECTORY, os.F_OK)):
        os.makedirs(DIRECTORY, mode=MODE, exist_ok=True)

    #-- list object of output files for file logs (full path)
    output_files = []

    #-- file information
    suffix = ['txt', 'nc', 'H5'][DATAFORM-1]
    format_str = ['ascii','netCDF4','HDF5'][DATAFORM-1]

    #-- read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = load_love_numbers(LMAX, REFERENCE='CF')

    #-- Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gauss_weights(RAD,LMAX)
        gw_str = '_r{0:0.0f}km'.format(RAD)
    else:
        #-- else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    #-- flag for spherical harmonic order
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''

    #-- reading GRACE months for input range with grace_input_months.py
    #-- replacing SLR and Degree 1 if specified
    #-- correcting for Pole-Tide and Atmospheric Jumps if specified
    Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        start_mon, end_mon, missing, SLR_C20, DEG1, MMAX=MMAX, SLR_C30=SLR_C30,
        MODEL_DEG1=MODEL_DEG1, ATM=ATM, POLE_TIDE=POLE_TIDE)
    #-- convert to harmonics object and remove mean if specified
    GRACE_Ylms = harmonics().from_dict(Ylms)
    GRACE_Ylms.directory = Ylms['directory']
    GRACE_Ylms.title = Ylms['title']
    #-- mean parameters
    MEAN = parameters['MEAN'] in ('Y','y')
    #-- use a mean file for the static field to remove
    if (parameters['MEAN_FILE'].title() == 'None'):
        mean_Ylms = GRACE_Ylms.mean(apply=MEAN)
    else:
        #-- data form for input mean file (1: ascii, 2: netcdf, 3: HDF5)
        MEANFORM = np.int(parameters['MEANFORM'])
        if (MEANFORM == 1):
            mean_Ylms=harmonics().from_ascii(parameters['MEAN_FILE'],date=False)
        if (MEANFORM == 2):
            mean_Ylms=harmonics().from_netCDF4(parameters['MEAN_FILE'],date=False)
        if (MEANFORM == 3):
            mean_Ylms=harmonics().from_HDF5(parameters['MEAN_FILE'],date=False)
        #-- remove the input mean
        if MEAN:
            GRACE_Ylms.subtract(mean_Ylms)
    #-- date information of GRACE/GRACE-FO coefficients
    nfiles = len(GRACE_Ylms.time)

    #-- filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        #-- destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
    else:
        #-- using standard GRACE/GRACE-FO harmonics
        ds_str = ''

    #-- input GIA spherical harmonic datafiles
    GIA_Ylms_rate = read_GIA_model(GIA_FILE,GIA=GIA,LMAX=LMAX,MMAX=MMAX)
    #-- calculate the monthly mass change from GIA
    GIA_Ylms = GRACE_Ylms.zeros_like()
    GIA_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    GIA_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    #-- monthly GIA calculated by gia_rate*time elapsed
    #-- finding change in GIA each month
    for t in range(nfiles):
        GIA_Ylms.clm[:,:,t] = GIA_Ylms_rate['clm']*(GIA_Ylms.time[t]-2003.3)
        GIA_Ylms.slm[:,:,t] = GIA_Ylms_rate['slm']*(GIA_Ylms.time[t]-2003.3)

    #-- Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_REMOVED:
        #-- read Land-Sea Mask and convert to spherical harmonics
        LANDMASK = os.path.expanduser(parameters['LANDMASK'])
        ocean_Ylms = ocean_stokes(LANDMASK,LMAX,MMAX=MMAX,LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        ocean_str = ''

    #-- input spherical harmonic datafiles to be removed from the GRACE data
    #-- Remove sets of Ylms from the GRACE data before returning
    remove_Ylms = GRACE_Ylms.zeros_like()
    remove_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    remove_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    if (parameters['REMOVE_INDEX'].title() != 'None'):
        #-- for each file index separated by commas
        for REMOVE_INDEX in parameters['REMOVE_INDEX'].split(','):
            Ylms = harmonics().from_index(REMOVE_INDEX, format_str)
            #-- reduce to GRACE/GRACE-FO months and truncate to degree and order
            Ylms = Ylms.subset(GRACE_Ylms.month).truncate(lmax=LMAX,mmax=MMAX)
            #-- distribute removed Ylms uniformly over the ocean
            if REDISTRIBUTE_REMOVED:
                #-- calculate ratio between total removed mass and
                #-- a uniformly distributed cm of water over the ocean
                ratio = Ylms.clm[0,0,:]/ocean_Ylms['clm'][0,0]
                #-- for each spherical harmonic
                for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
                    for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                        #-- remove the ratio*ocean Ylms from Ylms
                        #-- note: x -= y is equivalent to x = x - y
                        Ylms.clm[l,m,:] -= ratio*ocean_Ylms['clm'][l,m]
                        Ylms.slm[l,m,:] -= ratio*ocean_Ylms['slm'][l,m]
            #-- filter removed coefficients
            if DESTRIPE:
                Ylms = Ylms.destripe()
            #-- add data for month t and INDEX_FILE to the total
            #-- remove_clm and remove_slm matrices
            #-- redistributing the mass over the ocean if specified
            remove_Ylms.add(Ylms)

    #-- Output spatial data object
    grid = spatial()
    #-- Output Degree Spacing
    if (np.ndim(DDEG) == 0):
        #-- dlon == dlat
        dlon = DDEG
        dlat = DDEG
    else:
        #-- dlon != dlat
        dlon = DDEG[0]
        dlat = DDEG[1]
    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (-180:180,90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
        grid.lon = -180 + dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        #-- (Degree spacing)/2
        grid.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)

    #-- Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM,dPLM = plm_holmes(LMAX,np.cos(theta))

    #-- Earth Parameters
    factors = units(lmax=LMAX).harmonic(hl,kl,ll)
    #-- output spatial units
    unit_list = ['cmwe', 'mmGH', 'mmCU', u'\u03BCGal', 'mbar']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']
    #-- dfactor is the degree dependent coefficients
    #-- for specific spherical harmonic output units
    if (UNITS == 1):
        #-- 1: cmwe, centimeters water equivalent
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).cmwe
    elif (UNITS == 2):
        #-- 2: mmGH, millimeters geoid height
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmGH
    elif (UNITS == 3):
        #-- 3: mmCU, millimeters elastic crustal deformation
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mmCU
    elif (UNITS == 4):
        #-- 4: micGal, microGal gravity perturbations
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).microGal
    elif (UNITS == 5):
        #-- 5: mbar, millibar equivalent surface pressure
        dfactor = units(lmax=LMAX).harmonic(hl,kl,ll).mbar

    #-- output file format
    file_format = '{0}{1}_L{2:d}{3}{4}{5}_{6:03d}.{7}'
    #-- converting harmonics to truncated, smoothed coefficients in units
    #-- combining harmonics to calculate output spatial fields
    for i,grace_month in enumerate(GRACE_Ylms.month):
        #-- GRACE/GRACE-FO harmonics for time t
        Ylms = GRACE_Ylms.index(i)
        #-- Remove GIA rate for time
        Ylms.subtract(GIA_Ylms.index(i))
        #-- Remove monthly files to be removed
        Ylms.subtract(remove_Ylms.index(i))
        #-- smooth harmonics and convert to output units
        Ylms.convolve(dfactor*wt)
        #-- convert spherical harmonics to output spatial grid
        grid.data = harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, MMAX=MMAX, PLM=PLM).T
        #-- copy time variables for month
        grid.time = np.copy(Ylms.time)
        grid.month = np.copy(Ylms.month)

        #-- output monthly files to ascii, netCDF4 or HDF5
        args=(FILENAME,unit_list[UNITS-1],LMAX,order_str,gw_str,
            ds_str,grace_month,suffix)
        FILE=os.path.join(DIRECTORY,file_format.format(*args))
        if (DATAFORM == 1):
            #-- ascii (.txt)
            grid.to_ascii(FILE, date=True, verbose=VERBOSE)
        elif (DATAFORM == 2):
            #-- netCDF4
            grid.to_netCDF4(FILE, date=True, verbose=VERBOSE,
                units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
                title='GRACE/GRACE-FO Spatial Data')
        elif (DATAFORM == 3):
            #-- HDF5
            grid.to_HDF5(FILE, date=True, verbose=VERBOSE,
                units=unit_list[UNITS-1], longname=unit_name[UNITS-1],
                title='GRACE/GRACE-FO Spatial Data')
        #-- set the permissions mode of the output files
        os.chmod(FILE, MODE)
        #-- add file to list
        output_files.append(FILE)

    #-- return the list of output files
    return output_files

#-- PURPOSE: print a file log for the GRACE analysis
#-- lists: the parameter file, the parameters and the output files
def output_log_file(parameters,output_files):
    #-- format: GRACE_processing_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- check if run from entering parameters or from parameter files
    if parameters['PARAMETER_FILE'] is not None:
        #-- print parameter file on top
        print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
            os.path.abspath(parameters['PARAMETER_FILE'])), file=fid)
    else:
        print('PARAMETERS:', file=fid)
    #-- print parameter values sorted alphabetically
    for p in sorted(list(set(parameters.keys())-set(['PARAMETER_FILE']))):
        print('{0}: {1}'.format(p, parameters[p]), file=fid)
    #-- print output files
    print('\n\nOUTPUT FILES:', file=fid)
    for f in output_files:
        print('{0}'.format(f), file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE analysis
#-- lists: the parameter file, the parameters and the error
def output_error_log_file(parameters):
    #-- format: GRACE_processing_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- check if run from entering parameters or from parameter files
    if parameters['PARAMETER_FILE'] is not None:
        #-- print parameter file on top
        print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
            os.path.abspath(parameters['PARAMETER_FILE'])), file=fid)
    else:
        print('PARAMETERS:', file=fid)
    #-- print parameter values sorted alphabetically
    for p in sorted(list(set(parameters.keys())-set(['PARAMETER_FILE']))):
        print('{0}: {1}'.format(p, parameters[p]), file=fid)
    #-- print traceback error
    print('\n\nTRACEBACK ERROR:', file=fid)
    traceback.print_exc(file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: open a unique log file adding a numerical instance if existing
def create_unique_logfile(filename):
    #-- split filename into fileBasename and fileExtension
    fileBasename, fileExtension = os.path.splitext(filename)
    #-- create counter to add to the end of the filename if existing
    counter = 1
    while counter:
        try:
            #-- open file descriptor only if the file doesn't exist
            fd = os.open(filename, os.O_CREAT | os.O_EXCL | os.O_RDWR)
        except OSError:
            pass
        else:
            return os.fdopen(fd, 'w+')
        #-- new filename adds counter the between fileBasename and fileExtension
        filename = '{0}_{1:d}{2}'.format(fileBasename, counter, fileExtension)
        counter += 1

#-- PURPOSE: define the analysis for multiprocessing
def define_analysis(parameter_file,base_dir,LOG,VERBOSE,MODE):
    #-- keep track of multiprocessing threads
    info(os.path.basename(parameter_file))

    #-- variable with parameter definitions
    parameters = {}
    parameters['PARAMETER_FILE'] = parameter_file
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
        #-- run GRACE/GRACE-FO spatial algorithm with parameters
        output_files = grace_spatial_maps(base_dir,parameters,VERBOSE,MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
        if LOG:#-- write failed job completion log file
            output_error_log_file(parameters)
    else:
        if LOG:#-- write successful job completion log file
            output_log_file(parameters,output_files)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -P X, --np=X\t\tRun in parallel with X number of processes')
    print(' -V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermissions mode of the files created')
    print(' -l, --log\t\tOutput log file for each job')
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
    print('    Successful log file format: {0}'.format(LOGFILE))
    LOGFILE = 'GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
    print('    Failed log file format: {0}\n'.format(LOGFILE))

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program and run the analyses
    #-- with the specific parameters.
    long_options = ['help','directory=','np=','log','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:P:lVM:',long_options)

    #-- command line parameters
    base_dir = os.getcwd()
    PROCESSES = 0
    LOG = False
    #-- verbose output of processing run
    VERBOSE = False
    #-- permissions mode of the files created (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("-P","--np"):
            PROCESSES = np.int(arg)
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- raise exception if no parameter files entered
    if not arglist:
        raise IOError('No Parameter File Specified')

    #-- use parameter files from system arguments listed after the program.
    if (PROCESSES == 0):
        #-- run directly as series if PROCESSES = 0
        for f in arglist:
            define_analysis(os.path.expanduser(f),base_dir,LOG,VERBOSE,MODE)
    else:
        #-- run in parallel with multiprocessing Pool
        pool = multiprocessing.Pool(processes=PROCESSES)
        #-- for each parameter file
        for f in arglist:
            pool.apply_async(define_analysis, args=(os.path.expanduser(f),
                base_dir,LOG,VERBOSE,MODE))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()

#-- run main program
if __name__ == '__main__':
    main()
