#!/usr/bin/env python
u"""
grace_mean_harmonics.py
Written by Tyler Sutterley (04/2020)

Calculates the temporal mean of the GRACE/GRACE-FO spherical harmonics
    for a given date range from a set of parameters

CALLING SEQUENCE:
    Uses input parameter files:
    python grace_mean_harmonics.py parameter_file

    Can also input several parameter files in series:
    python grace_mean_harmonics.py parameter_file1 parameter_file2

    Can be run in parallel with the python multiprocessing package:
    python grace_mean_harmonics.py --np=2 parameter_file1 parameter_file2
    python grace_mean_harmonics.py -P 2 parameter_file1 parameter_file2

    Can output a log file listing the input parameters and output files:
    python grace_mean_harmonics.py --log parameter_file
    python grace_mean_harmonics.py -l parameter_file

SYSTEM ARGUMENTS README:
    program is run as:
    python grace_mean_harmonics.py inp1 inp2 inp3
        where inp1, inp2 and inp3 are different inputs

        firstinput=sys.argv[1] (in this case inp1)
        secondinput=sys.argv[2] (in this case inp2)
        thirdinput=sys.argv[3] (in this case inp3)

    As python is base 0, sys.argv[0] is equal to grace_mean_harmonics.py
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
    or parameters entered from the command line

COMMAND LINE OPTIONS:
    -D X, --directory=X: GRACE/GRACE-FO working data directory
    -P X, --np=X: Run in parallel with X number of processes
    -l, --log: Output a log file listing output files

PYTHON DEPENDENCIES:
    Python: a general-purpose, high-level programming language
        http://www.python.org/
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: can replace C30 with coefficients from SLR
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 03/2018: edited output netCDF4 and HDF5 attributes
    Updated 04/2017: changed no input file exception to IOError
    Updated 05/2016: using __future__ print function
    Updated 02/2016: use getopt parameters to set options for the
        number of processes to run in parallel and to output a log file
    Updated 12/2015: output unique log file
    Updated 08/2015: changed sys.exit to a raise exception instance
        Added pole tide and GAE/GAF/GAG correction parameters
    Updated 05/2015: added MMAX parameter to study new 2015 60X30 fields
        DATAFORM now split for INPUT_DATAFORM and OUTPUT_DATAFORM
    Updated 01/2015: added error handling for multiprocessing threads
    Updated 11/2014: added HDF5 dataform option
    Updated 10/2014: updated for distributed computing of tasks
        with the multiprocessing module
    Written 05/2014
"""
from __future__ import print_function

import sys
import os
import time
import getopt
import numpy as np
import multiprocessing
import traceback
from gravity_toolkit.grace_input_months import grace_input_months
from gravity_toolkit.harmonics import harmonics

#-- PURPOSE: keep track of multiprocessing threads
def info(title):
    print(os.path.basename(sys.argv[0]))
    print(title)
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: import GRACE/GRACE-FO files for a given months range
#-- calculate the mean of the spherical harmonics and output to file
def grace_mean_harmonics(base_dir, parameters):
    #-- Data processing center
    PROC = parameters['PROC']
    #-- Data Release
    DREL = parameters['DREL']
    #-- GRACE dataset
    DSET = parameters['DSET']
    #-- pole tide corrections from Wahr et al (2015)
    POLE_TIDE = parameters['POLE_TIDE'] in ('Y','y')
    #-- ATM corrections
    ATM = parameters['ATM'] in ('Y','y')

    #-- maximum degree and order
    LMAX = np.int(parameters['LMAX'])
    #-- maximum spherical harmonic order
    if (parameters['MMAX'].title() == 'None'):
        MMAX = np.copy(LMAX)
    else:
        MMAX = np.int(parameters['MMAX'])

    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''

    #-- Date Range and missing months
    START_MON = np.int(parameters['START'])
    END_MON = np.int(parameters['END'])
    MISSING = np.array(parameters['MISSING'].split(','),dtype=np.int)
    #-- SLR C2,0 and C3,0
    SLR_C20 = parameters['SLR_C20']
    SLR_C30 = parameters['SLR_C30']
    #-- Degree 1 correction
    DEG1 = parameters['DEG1']

    #-- data formats for output: 1) ascii, 2) netcdf, 3) HDF5
    DATAFORM = np.int(parameters['DATAFORM'])
    suffix = ['txt', 'nc', 'H5'][DATAFORM-1]

    #-- reading GRACE months for input range with grace_input_months.py
    #-- replacing C20 and C30 with SLR values and updating Degree 1 if specified
    #-- correcting for Pole Tide Drift and Atmospheric Jumps if specified
    input_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START_MON, END_MON, MISSING, SLR_C20, DEG1, MMAX=MMAX, MODEL_DEG1=False,
        SLR_C30=SLR_C30, POLE_TIDE=POLE_TIDE, ATM=ATM)
    grace_Ylms = harmonics().from_dict(input_Ylms)
    #-- output directory
    if (parameters['DIRECTORY'].title() == 'None'):
        #-- if not entering via parameter file use directory for product
        DIRECTORY = os.path.expanduser(input_Ylms['directory'])
    else:
        DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- descriptor string for processing parameters
    grace_str = input_Ylms['title']
    #-- calculate mean Ylms
    mean_Ylms = grace_Ylms.mean()
    mean_Ylms.time = np.mean(grace_Ylms.time)
    mean_Ylms.month = np.mean(grace_Ylms.month)

    #-- default output filename if not entering via parameter file
    if (parameters['FILENAME'].title() == 'None'):
        file_format = '{0}_{1}_{2}_MEAN_CLM{3}_L{4:d}{5}_{6:03d}-{7:03d}.{8}'
        FILENAME = file_format.format(PROC, DREL, DSET, grace_str, LMAX,
            order_str, START_MON, END_MON, suffix)
    else:
        FILENAME = parameters['FILENAME']

    #-- output spherical harmonics for the static field
    if (DATAFORM == 1):
        #-- output mean field to ascii
        mean_Ylms.to_ascii(os.path.join(DIRECTORY,FILENAME))
    elif (DATAFORM == 2):
        #-- output mean field to netCDF4
        mean_Ylms.to_netCDF4(os.path.join(DIRECTORY,FILENAME))
    elif (DATAFORM == 3):
        #-- output mean field to HDF5
        mean_Ylms.to_HDF5(os.path.join(DIRECTORY,FILENAME))

    #-- return the output file
    return os.path.join(DIRECTORY,FILENAME)

#-- PURPOSE: print a file log for the GRACE/GRACE-FO mean program
#-- lists: the parameter file, the parameters and the output file
def output_log_file(parameters,output_file):
    #-- format: GRACE_mean_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- print parameter file on top
    print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
        os.path.abspath(parameters['PARAMETER_FILE'])),file=fid)
    #-- print parameter values sorted alphabetically
    for p in sorted(list(set(parameters.keys())-set(['PARAMETER_FILE']))):
        print('{0}: {1}'.format(p, parameters[p]), file=fid)
    #-- print output files
    print('\n\nOUTPUT FILE:',file=fid)
    print('{0}'.format(output_file),file=fid)
    #-- close the log file
    fid.close()

#-- PURPOSE: print a error file log for the GRACE/GRACE-FO mean program
#-- lists: the parameter file, the parameters and the error
def output_error_log_file(parameters):
    #-- format: GRACE_mean_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_failed_run_{0}_PID-{1:d}.log'.format(*args)
    DIRECTORY = os.path.expanduser(parameters['DIRECTORY'])
    #-- create a unique log and open the log file
    fid = create_unique_logfile(os.path.join(DIRECTORY,LOGFILE))
    #-- print parameter file on top
    print('PARAMETER FILE:\n{0}\n\nPARAMETERS:'.format(
        os.path.abspath(parameters['PARAMETER_FILE'])),file=fid)
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
def define_analysis(parameter_file,base_dir,LOG):
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
        #-- run mean algorithm with parameters
        output_file = grace_mean_harmonics(base_dir, parameters)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
        if LOG:#-- write failed job completion log file
            output_error_log_file(parameters)
    else:
        #-- write successful job completion log file
        if LOG:
            output_log_file(parameters,output_file)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking GRACE/GRACE-FO data directory')
    print(' -P X, --np=X\tRun in parallel with X number of processes')
    print(' -l, --log\tOutput log file for each job')
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_mean_run_{0}_PID-{1:d}.log'.format(*args)
    print('    Successful log file format: {0}'.format(LOGFILE))
    FAILED_LOGFILE = 'GRACE_mean_failed_run_{0}_PID-{1:d}.log'.format(*args)
    print('    Failed log file format: {0}\n'.format(FAILED_LOGFILE))

#-- This is the main part of the program that calls the individual modules
#-- If no parameter file is listed as an argument: will exit with an error
def main():
    #-- Read the system arguments listed after the program and run the analyses
    #-- with the specific parameters
    long_options = ['help','directory=','np=','log']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:P:l',long_options)

    #-- command line parameters
    base_dir = os.getcwd()
    PROCESSES = 0
    LOG = False
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

    #-- raise exception if no parameter files entered
    if not arglist:
        raise IOError('No Parameter File Specified')

    if (PROCESSES == 0):
        #-- run directly as series if PROCESSES = 0
        for parameter_file in arglist:
            define_analysis(parameter_file,base_dir,LOG)
    else:
        #-- run in parallel with multiprocessing Pool
        pool = multiprocessing.Pool(processes=PROCESSES)
        #-- for each parameter file
        for fi in arglist:
            pool.apply_async(define_analysis, args=(fi,base_dir,LOG))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()

#-- run main program
if __name__ == '__main__':
    main()
