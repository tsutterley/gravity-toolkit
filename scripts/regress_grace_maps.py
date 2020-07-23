#!/usr/bin/env python
u"""
regress_grace_maps.py
Written by Tyler Sutterley (06/2020)

Reads in GRACE/GRACE-FO spatial files from grace_spatial_maps.py and
    fits a regression model at each grid point

CALLING SEQUENCE:
    python regress_grace_maps.py --start=4 --end=175 --order=2 parameter_file

COMMAND LINE OPTIONS:
    -P X, --np=X: Run in parallel with X number of processes
    -S X, --start=X: starting GRACE month for time series regression
    -E X, --end=X: ending GRACE month for time series regression
    --order=X: regression fit polynomial order
    --cycles=X: regression fit cyclical terms
    -M X, --mode=X: permissions mode of the output files
    -V, --verbose: verbose output of processing run
    -l, --log: output log file for each job

SYSTEM ARGUMENTS README:
    program is run as:
    python regress_grace_maps.py parameter_file

    As python is base 0, sys.argv[0] is equal to regress_grace_maps.py
        (which is useful in some applications, but not for this program)

    For this program, the system arguments are parameter files
    The program reads the parameter file, which is separated by column as:
        Column 1: parameter name (such as LMAX)
        Column 2: parameter (e.g. 60)
        Column 3: comments (which are discarded)
    The parameters are stored in a python dictionary (variables indexed by keys)
        the keys are the parameter name (for LMAX: parameters['LMAX'] == 60)

INPUTS:
    parameter file containing specific variables for the analysis

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py:  Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    tsregress.py: calculates trend coefficients using least-squares
    tsamplitude.py: calculates the amplitude and phase of a harmonic function
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5

UPDATE HISTORY:
    Updated 06/2020: using spatial data class for input and output operations
    Updated 01/2020: output seasonal amplitude and phase
    Updated 10/2019: changing Y/N flags to True/False
    Updated 01/2019: include 161-day S2 tidal aliasing terms in regression
    Updated 12/2018: added parallel processing with multiprocessing
    Updated 11/2018: using future division for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
        can run for different sets of months using the --start and --end options
        can run for different fit types using the --fit option
        use output_data wrapper function for writing data to file
    Updated 04/2018: changed the names of the output log files
    Updated 03/2018: added option --mode to set the output file permissions
        added option for output in pascals (UNITS=5)
    Updated 06/2016: using __future__ print function, MMAX for LMAX != MMAX
    Updated 06/2015: added output_files for log files
    Written 09/2013
"""
from __future__ import print_function, division

import sys
import os
import time
import getopt
import inspect
import traceback
import numpy as np
import multiprocessing

from gravity_toolkit.tsregress import tsregress
from gravity_toolkit.tsamplitude import tsamplitude
from gravity_toolkit.spatial import spatial

#-- PURPOSE: keep track of multiprocessing threads
def info(title):
    print(os.path.basename(sys.argv[0]))
    print(title)
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- program module to run with specified parameters
def regress_grace_maps(parameters,ORDER,CYCLES,VERBOSE,MODE):
    #-- convert parameters into variables
    #-- Data processing center
    PROC = parameters['PROC']
    #-- Data Release
    DREL = parameters['DREL']
    #-- GRACE/GRACE-FO dataset
    DSET = parameters['DSET']
    #-- GRACE/GRACE-FO months
    START_MON = np.int(parameters['START'])
    END_MON = np.int(parameters['END'])
    MISSING = np.array(parameters['MISSING'].split(','),dtype=np.int)
    months = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    nmon = len(months)
    #-- maximum degree and order
    LMAX = np.int(parameters['LMAX'])
    if (parameters['MMAX'].title() == 'None'):
        MMAX = np.copy(LMAX)
    else:
        MMAX = np.int(parameters['MMAX'])
    #-- Glacial Isostatic Adjustment file to read
    GIA = parameters['GIA'] if (parameters['GIA'].title() != 'None') else None
    GIA_FILE = os.path.expanduser(parameters['GIA_FILE'])
    #-- remove a set of spherical harmonics from the GRACE data
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
    #-- output filename suffix
    suffix = ['txt', 'nc', 'H5'][DATAFORM-1]

    #-- flag for spherical harmonic order
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- Calculating the Gaussian smoothing for radius RAD
    gw_str = '_r{0:0.0f}km'.format(RAD) if (RAD != 0) else ''
    #-- destriped GRACE/GRACE-FO coefficients
    ds_str = '_FL' if DESTRIPE else ''
    #-- distributing removed mass uniformly over ocean
    ocean_str = '_OCN' if REDISTRIBUTE_REMOVED else ''
    #-- input and output spatial units
    unit_list = ['cmwe', 'mmGH', 'mmCU', u'\u03BCGal', 'mbar']
    unit_name = ['Equivalent Water Thickness', 'Geoid Height',
        'Elastic Crustal Uplift', 'Gravitational Undulation',
        'Equivalent Surface Pressure']

    #-- input file format
    input_format = '{0}{1}_L{2:d}{3}{4}{5}_{6:03d}.{7}'
    #-- output file format
    output_format = '{0}{1}_L{2:d}{3}{4}{5}_{6}{7}_{8:03d}-{9:03d}.{10}'

    #-- Output Degree Spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    #-- Output Degree Interval
    if (INTERVAL == 1):
        #-- (-180:180,90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
    elif (INTERVAL == 2):
        #-- (Degree spacing)/2
        nlon = np.int(360.0/dlon)
        nlat = np.int(180.0/dlat)

    #-- Setting output parameters for each fit type
    coef_str = ['x{0:d}'.format(o) for o in range(ORDER+1)]
    unit_suffix = [' yr^{0:d}'.format(-o) if o else '' for o in range(ORDER+1)]
    if (ORDER == 0):#-- Mean
        unit_longname = ['Mean']
    elif (ORDER == 1):#-- Trend
        unit_longname = ['Constant','Trend']
    elif (ORDER == 2):#-- Quadratic
        unit_longname = ['Constant','Linear','Quadratic']
    #-- filename strings for cyclical terms
    cyclic_str = {}
    cyclic_str['SEMI'] = ['SS','SC']
    cyclic_str['ANN'] = ['AS','AC']
    cyclic_str['S2'] = ['S2S','S2C']
    #-- unit longnames for cyclical terms
    cyclic_longname = {}
    cyclic_longname['SEMI'] = ['Semi-Annual Sine', 'Semi-Annual Cosine']
    cyclic_longname['ANN'] = ['Annual Sine', 'Annual Cosine']
    cyclic_longname['S2'] = ['S2 Tidal Alias Sine', 'S2 Tidal Alias Cosine']
    amp_str = []
    for i,c in enumerate(CYCLES):
        if (c == 0.5):
            flag = 'SEMI'
        elif (c == 1.0):
            flag = 'ANN'
        elif (c == (161.0/365.25)):
            flag = 'S2'
        coef_str.extend(cyclic_str[flag])
        unit_longname.extend(cyclic_longname[flag])
        unit_suffix.extend(['',''])
        amp_str.append(flag)

    #-- input data spatial object
    spatial_list = []
    for t,grace_month in enumerate(months):
        #-- input GRACE/GRACE-FO spatial file
        fi = input_format.format(FILENAME,unit_list[UNITS-1],LMAX,
            order_str,gw_str,ds_str,grace_month,suffix)
        #-- read GRACE/GRACE-FO spatial file
        if (DATAFORM == 1):
            dinput = spatial(spacing=[dlon,dlat],nlon=nlon,
                nlat=nlat).from_ascii(os.path.join(DIRECTORY,fi))
        elif (DATAFORM == 2):
            #-- netcdf (.nc)
            dinput = spatial().from_netCDF4(os.path.join(DIRECTORY,fi))
        elif (DATAFORM == 3):
            #-- HDF5 (.H5)
            dinput = spatial().from_HDF5(os.path.join(DIRECTORY,fi))
        #-- append to spatial list
        dinput.month[:] = grace_month
        nlat,nlon = dinput.shape
        spatial_list.append(dinput)

    #-- concatenate list to single spatial object
    grid = spatial().from_list(spatial_list)
    spatial_list = None

    #-- Fitting seasonal components
    ncomp = len(coef_str)
    ncycles = 2*len(CYCLES)

    #-- Allocating memory for output variables
    out = dinput.zeros_like()
    out.data = np.zeros((nlat,nlon,ncomp))
    out.error = np.zeros((nlat,nlon,ncomp))
    out.mask = np.ones((nlat,nlon,ncomp),dtype=np.bool)
    #-- Fit Significance
    FS = {}
    #-- SSE: Sum of Squares Error
    #-- AIC: Akaike information criterion
    #-- BIC: Bayesian information criterion
    #-- R2Adj: Adjusted Coefficient of Determination
    for key in ['SSE','AIC','BIC','R2Adj']:
        FS[key] = dinput.zeros_like()

    #-- valid values for ocean function
    #-- calculate the regression coefficients and fit significance
    for i in range(nlat):
        for j in range(nlon):
            #-- Calculating the regression coefficients
            tsbeta = tsregress(grid.time, grid.data[i,j,:],
                ORDER=ORDER, CYCLES=CYCLES, CONF=0.95)
            #-- save regression components
            for k in range(0, ncomp):
                out.data[i,j,k] = tsbeta['beta'][k]
                out.error[i,j,k] = tsbeta['error'][k]
                out.mask[i,j,k] = False
            #-- Fit significance terms
            #-- Degrees of Freedom
            nu = tsbeta['DOF']
            #-- Converting Mean Square Error to Sum of Squares Error
            FS['SSE'].data[i,j] = tsbeta['MSE']*nu
            FS['AIC'].data[i,j] = tsbeta['AIC']
            FS['BIC'].data[i,j] = tsbeta['BIC']
            FS['R2Adj'].data[i,j] = tsbeta['R2Adj']

    #-- list of output files
    output_files = []
    #-- Output spatial files
    for i in range(0,ncomp):
        #-- output spatial file name
        f1 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
            gw_str,ds_str,coef_str[i],'',START_MON,END_MON,suffix)
        f2 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
            gw_str,ds_str,coef_str[i],'_ERROR',START_MON,END_MON,suffix)
        #-- full attributes
        UNITS_TITLE = '{0}{1}'.format(unit_list[UNITS-1],unit_suffix[i])
        LONGNAME = unit_name[UNITS-1]
        FILE_TITLE = 'GRACE/GRACE-FO_Spatial_Data_{0}'.format(unit_longname[i])
        #-- output regression fit to file
        output = out.index(i, date=False)
        output_data(output, FILENAME=os.path.join(DIRECTORY,f1),
            DATAFORM=DATAFORM, UNITS=UNITS_TITLE, LONGNAME=LONGNAME,
            TITLE=FILE_TITLE, VERBOSE=VERBOSE, MODE=MODE)
        output_data(output, FILENAME=os.path.join(DIRECTORY,f2),
            DATAFORM=DATAFORM, UNITS=UNITS_TITLE, LONGNAME=LONGNAME,
            TITLE=FILE_TITLE, KEY='error', VERBOSE=VERBOSE, MODE=MODE)
        #-- add output files to list object
        output_files.append(os.path.join(DIRECTORY,f1))
        output_files.append(os.path.join(DIRECTORY,f2))

    #-- if fitting coefficients with cyclical components
    if (ncycles > 0):
        #-- output spatial titles for amplitudes
        amp_title = {'ANN':'Annual Amplitude','SEMI':'Semi-Annual Amplitude',
            'S2':'S2 Tidal Alias Amplitude'}
        ph_title = {'ANN':'Annual Phase','SEMI':'Semi-Annual Phase',
            'S2':'S2 Tidal Alias Phase'}

        #-- output amplitude and phase of cyclical components
        for i,flag in enumerate(amp_str):
            #-- Indice pointing to the cyclical components
            j = 1 + ORDER + 2*i
            #-- Allocating memory for output amplitude and phase
            amp = dinput.zeros_like()
            ph = dinput.zeros_like()
            #-- calculating amplitude and phase of spatial field
            amp.data,ph.data = tsamplitude(out.data[:,:,j],out.data[:,:,j+1])
            #-- convert phase from -180:180 to 0:360
            ii,jj = np.nonzero(ph.data < 0)
            ph.data[ii,jj] += 360.0
            #-- Amplitude Error
            comp1 = out.error[:,:,j]*out.data[:,:,j]/amp.data
            comp2 = out.error[:,:,j+1]*out.data[:,:,j+1]/amp.data
            amp.error = np.sqrt(comp1**2 + comp2**2)
            #-- Phase Error (degrees)
            comp1 = out.error[:,:,j]*out.data[:,:,j+1]/(amp.data**2)
            comp2 = out.error[:,:,j+1]*out.data[:,:,j]/(amp.data**2)
            ph.error = (180.0/np.pi)*np.sqrt(comp1**2 + comp2**2)

            #-- output file names for amplitude, phase and errors
            f3 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
                gw_str,ds_str,flag,'',START_MON,END_MON,suffix)
            f4 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
                gw_str,ds_str,flag,'_PHASE',START_MON,END_MON,suffix)
            #-- output spatial error file name
            f5 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
                gw_str,ds_str,flag,'_ERROR',START_MON,END_MON,suffix)
            f6 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
                gw_str,ds_str,flag,'_PHASE_ERROR',START_MON,END_MON,suffix)
            #-- full attributes
            AMP_UNITS = unit_list[UNITS-1]
            PH_UNITS = 'degrees'
            LONGNAME = unit_name[UNITS-1]
            AMP_TITLE = 'GRACE/GRACE-FO_Spatial_Data_{0}'.format(amp_title[flag])
            PH_TITLE = 'GRACE/GRACE-FO_Spatial_Data_{0}'.format(ph_title[flag])
            #-- Output seasonal amplitude and phase to files
            output_data(amp, FILENAME=os.path.join(DIRECTORY,f3),
                DATAFORM=DATAFORM, UNITS=AMP_UNITS, LONGNAME=LONGNAME,
                TITLE=AMP_TITLE, VERBOSE=VERBOSE, MODE=MODE)
            output_data(ph, FILENAME=os.path.join(DIRECTORY,f4),
                DATAFORM=DATAFORM, UNITS=PH_UNITS, LONGNAME='Phase',
                TITLE=PH_TITLE, VERBOSE=VERBOSE, MODE=MODE)
            #-- Output seasonal amplitude and phase error to files
            output_data(amp, FILENAME=os.path.join(DIRECTORY,f5),
                DATAFORM=DATAFORM, UNITS=AMP_UNITS, LONGNAME=LONGNAME,
                TITLE=AMP_TITLE, KEY='error', VERBOSE=VERBOSE, MODE=MODE)
            output_data(ph, FILENAME=os.path.join(DIRECTORY,f6),
                DATAFORM=DATAFORM, UNITS=PH_UNITS, LONGNAME='Phase',
                TITLE=PH_TITLE, KEY='error', VERBOSE=VERBOSE, MODE=MODE)
            #-- add output files to list object
            output_files.append(os.path.join(DIRECTORY,f3))
            output_files.append(os.path.join(DIRECTORY,f4))
            output_files.append(os.path.join(DIRECTORY,f5))
            output_files.append(os.path.join(DIRECTORY,f6))

    #-- Output fit significance
    signif_longname = {'SSE':'Sum of Squares Error',
        'AIC':'Akaike information criterion',
        'BIC':'Bayesian information criterion',
        'R2Adj':'Adjusted Coefficient of Determination'}
    #-- for each fit significance term
    for key,fs in FS.items():
        #-- output file names for fit significance
        signif_str = '{0}_'.format(key)
        f7 = output_format.format(FILENAME,unit_list[UNITS-1],LMAX,order_str,
            gw_str,ds_str,signif_str,coef_str[ORDER],START_MON,END_MON,suffix)
        #-- full attributes
        LONGNAME = signif_longname[key]
        #-- output fit significance to file
        output_data(fs, FILENAME=os.path.join(DIRECTORY,f7),
            DATAFORM=DATAFORM, UNITS=key, LONGNAME=LONGNAME,
            TITLE=nu, VERBOSE=VERBOSE, MODE=MODE)
        #-- add output files to list object
        output_files.append(os.path.join(DIRECTORY,f7))

    #-- return the list of output files
    return output_files

#-- PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, KEY='data', DATAFORM=None,
    UNITS=None, LONGNAME=None, TITLE=None, VERBOSE=False, MODE=0o775):
    output = data.copy()
    setattr(output,'data',getattr(data,KEY))
    if (DATAFORM == 1):
        #-- ascii (.txt)
        output.to_ascii(FILENAME,date=False,verbose=VERBOSE)
    elif (DATAFORM == 2):
        #-- netcdf (.nc)
        output.to_netCDF4(FILENAME,date=False,verbose=VERBOSE,
            units=UNITS,longname=LONGNAME,title=TITLE)
    elif (DATAFORM == 3):
        #-- HDF5 (.H5)
        output.to_HDF5(FILENAME,date=False,verbose=VERBOSE,
            units=UNITS,longname=LONGNAME,title=TITLE)
    #-- change the permissions mode of the output file
    os.chmod(FILENAME, MODE)

#-- PURPOSE: print a file log for the GRACE/GRACE-FO regression
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

#-- PURPOSE: print a error file log for the GRACE/GRACE-FO regression
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

#-- PURPOSE: define the analysis
def define_analysis(parameter_file,START,END,ORDER,CYCLES,VERBOSE,MODE,LOG):
    #-- keep track of progress
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

    #-- can change the date range from the defaults (change parameters for logs)
    parameters['START'] = np.copy(START) if START else parameters['START']
    parameters['END'] = np.copy(END) if END else parameters['END']
    #-- set regression fit type as parameter for logs
    parameters['ORDER'] = '{0:d}'.format(ORDER)
    parameters['CYCLES'] = ','.join(['{0:4f}'.format(c) for c in CYCLES])
    #-- try to run the program with listed parameters
    try:
        #-- run GRACE/GRACE-FO spatial regression program with parameters
        output_files = regress_grace_maps(parameters,ORDER,CYCLES,VERBOSE,MODE)
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
    print(' -P X, --np=X\tRun in parallel with X number of processes')
    print(' -S X, --start=X\tStarting GRACE month for time series regression')
    print(' -E X, --end=X\t\tEnding GRACE month for time series regression')
    print(' --order=X\t\tRegression fit polynomial order:')
    print('\t0: Mean\n\t1: Trend\n\t2: Quadratic')
    print(' --cycles=X\t\tRegression fit cyclical terms:')
    print(' -M X, --mode=X\t\tPermission mode of directories and files')
    print(' -V, --verbose\t\tVerbose output of processing run')
    print(' -l, --log\t\tOutput log file for each job')
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE='GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
    FAILED_LOGFILE='GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
    print('    Successful log file format: {0}'.format(LOGFILE))
    print('    Failed log file format: {0}\n'.format(FAILED_LOGFILE))

#-- This is the main part of the program that calls the individual modules
#-- If no parameter file is listed as an argument: will exit with an error
def main():
    #-- Read the system arguments listed after the program and run the analyses
    #-- with the specific parameters.
    lopt = ['help','np=','start=','end=','order=','cycles=','verbose','mode=','log']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hP:S:E:F:VM:l',lopt)

    #-- command line parameters
    PROCESSES = 0
    START = None
    END = None
    ORDER = 2
    CYCLES = np.array([0.5,1.0,161.0/365.25])
    LOG = False
    VERBOSE = False
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-P","--np"):
            PROCESSES = np.int(arg)
        elif opt in ("-S","--start"):
            START = np.int(arg)
        elif opt in ("-E","--end"):
            END = np.int(arg)
        elif opt in ("--order"):
            ORDER = np.int(arg)
        elif opt in ("--cycles"):
            CYCLES = np.array([eval(i) for i in arg.split(',')],dtype=np.float)
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg,8)
        elif opt in ("-l","--log"):
            LOG = True

    #-- raise exception if no parameter files entered
    if not arglist:
        raise IOError('No Parameter File Specified')

    #-- use parameter files from system arguments listed after the program.
    if (PROCESSES == 0):
        #-- run directly as series if PROCESSES = 0
        for f in arglist:
            define_analysis(os.path.expanduser(f),START,END,ORDER,CYCLES,
                VERBOSE,MODE,LOG)
    else:
        #-- run in parallel with multiprocessing Pool
        pool = multiprocessing.Pool(processes=PROCESSES)
        #-- for each parameter file
        for f in arglist:
            pool.apply_async(define_analysis,args=(os.path.expanduser(f),
                START,END,ORDER,CYCLES,VERBOSE,MODE,LOG))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()

#-- run main program
if __name__ == '__main__':
    main()
