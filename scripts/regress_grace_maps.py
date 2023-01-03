#!/usr/bin/env python
u"""
regress_grace_maps.py
Written by Tyler Sutterley (01/2023)

Reads in GRACE/GRACE-FO spatial files from grace_spatial_maps.py and
    fits a regression model at each grid point

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -P X, --file-prefix X: prefix string for input and output files
    -S X, --start X: starting GRACE month for time series regression
    -E X, --end X: ending GRACE month for time series regression
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -U X, --units X: output units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: mbar equivalent surface pressure
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    --redistribute-removed: redistribute removed mass fields over the ocean
    --order X: regression fit polynomial order
    --cycles X: regression fit cyclical terms
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    time_series.regress.py: calculates trend coefficients using least-squares
    time_series.amplitude.py: calculates the amplitude and phase of a harmonic
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
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
        added option for output in equivalent surface pressure (UNITS=5)
    Updated 06/2016: using __future__ print function, MMAX for LMAX != MMAX
    Updated 06/2015: added output_files for log files
    Written 09/2013
"""
from __future__ import print_function, division

import sys
import os
import time
import logging
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# program module to run with specified parameters
def regress_grace_maps(LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
    UNITS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    DATAFORM=None,
    REDISTRIBUTE_REMOVED=False,
    ORDER=None,
    CYCLES=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    VERBOSE=0,
    MODE=0o775):

    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]

    # flag for spherical harmonic order
    order_str = f'M{MMAX:d}' if MMAX and (MMAX != LMAX) else ''
    # Calculating the Gaussian smoothing for radius RAD
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # destriped GRACE/GRACE-FO coefficients
    ds_str = '_FL' if DESTRIPE else ''
    # distributing removed mass uniformly over ocean
    ocean_str = '_OCN' if REDISTRIBUTE_REMOVED else ''
    # input and output spatial units
    unit_list = ['cmwe', 'mmGH', 'mmCU', u'\u03BCGal', 'mbar']
    unit_name = ['Equivalent_Water_Thickness', 'Geoid_Height',
        'Elastic_Crustal_Uplift', 'Gravitational_Undulation',
        'Equivalent_Surface_Pressure']

    # input file format
    input_format = '{0}{1}_L{2:d}{3}{4}{5}_{6:03d}.{7}'
    # output file format
    output_format = '{0}{1}_L{2:d}{3}{4}{5}_{6}{7}_{8:03d}-{9:03d}.{10}'

    # GRACE months to read
    months = sorted(set(np.arange(START,END+1)) - set(MISSING))

    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        nlon = np.int64(360.0/dlon)
        nlat = np.int64(180.0/dlat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        lon = np.arange(minlon+dlon/2.0,maxlon+dlon/2.0,dlon)
        lat = np.arange(maxlat-dlat/2.0,minlat-dlat/2.0,-dlat)
        nlon = len(lon)
        nlat = len(lat)

    # Setting output parameters for each fit type
    coef_str = ['x{0:d}'.format(o) for o in range(ORDER+1)]
    unit_suffix = [' yr^{0:d}'.format(-o) if o else '' for o in range(ORDER+1)]
    if (ORDER == 0):# Mean
        unit_longname = ['Mean']
    elif (ORDER == 1):# Trend
        unit_longname = ['Constant','Trend']
    elif (ORDER == 2):# Quadratic
        unit_longname = ['Constant','Linear','Quadratic']
    # filename strings for cyclical terms
    cyclic_str = {}
    cyclic_str['SEMI'] = ['SS','SC']
    cyclic_str['ANN'] = ['AS','AC']
    cyclic_str['S2'] = ['S2S','S2C']
    # unit longnames for cyclical terms
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

    # input data spatial object
    spatial_list = []
    for t,grace_month in enumerate(months):
        # input GRACE/GRACE-FO spatial file
        fargs = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
            gw_str, ds_str, grace_month, suffix)
        input_file = os.path.join(OUTPUT_DIRECTORY,input_format.format(*fargs))
        # read GRACE/GRACE-FO spatial file
        if (DATAFORM == 'ascii'):
            dinput = gravtk.spatial(spacing=[dlon,dlat], nlon=nlon,
                nlat=nlat).from_ascii(input_file)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            dinput = gravtk.spatial().from_netCDF4(input_file)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            dinput = gravtk.spatial().from_HDF5(input_file)
        # append to spatial list
        dinput.month = grace_month
        nlat,nlon = dinput.shape
        spatial_list.append(dinput)

    # concatenate list to single spatial object
    grid = gravtk.spatial().from_list(spatial_list)
    spatial_list = None

    # Fitting seasonal components
    ncomp = len(coef_str)
    ncycles = 2*len(CYCLES)

    # Allocating memory for output variables
    out = dinput.zeros_like()
    out.data = np.zeros((nlat,nlon,ncomp))
    out.error = np.zeros((nlat,nlon,ncomp))
    out.mask = np.ones((nlat,nlon,ncomp),dtype=bool)
    # Fit Significance
    FS = {}
    # SSE: Sum of Squares Error
    # AIC: Akaike information criterion
    # BIC: Bayesian information criterion
    # R2Adj: Adjusted Coefficient of Determination
    for key in ['SSE','AIC','BIC','R2Adj']:
        FS[key] = dinput.zeros_like()

    # calculate the regression coefficients and fit significance
    for i in range(nlat):
        for j in range(nlon):
            # Calculating the regression coefficients
            tsbeta = gravtk.time_series.regress(grid.time, grid.data[i,j,:],
                ORDER=ORDER, CYCLES=CYCLES, CONF=0.95)
            # save regression components
            for k in range(0, ncomp):
                out.data[i,j,k] = tsbeta['beta'][k]
                out.error[i,j,k] = tsbeta['error'][k]
                out.mask[i,j,k] = False
            # Fit significance terms
            # Degrees of Freedom
            nu = tsbeta['DOF']
            # Converting Mean Square Error to Sum of Squares Error
            FS['SSE'].data[i,j] = tsbeta['MSE']*nu
            FS['AIC'].data[i,j] = tsbeta['AIC']
            FS['BIC'].data[i,j] = tsbeta['BIC']
            FS['R2Adj'].data[i,j] = tsbeta['R2Adj']

    # list of output files
    output_files = []
    # Output spatial files
    for i in range(0,ncomp):
        # output spatial file name
        f1 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
            gw_str, ds_str, coef_str[i], '', START, END, suffix)
        f2 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
            gw_str, ds_str, coef_str[i], '_ERROR', START, END, suffix)
        file1 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f1))
        file2 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f2))
        # full attributes
        UNITS_TITLE = f'{unit_list[UNITS-1]}{unit_suffix[i]}'
        LONGNAME = unit_name[UNITS-1]
        FILE_TITLE = f'GRACE/GRACE-FO_Spatial_Data_{unit_longname[i]}'
        # output regression fit to file
        output = out.index(i, date=False)
        output_data(output, FILENAME=file1, DATAFORM=DATAFORM,
            UNITS=UNITS_TITLE, LONGNAME=LONGNAME, TITLE=FILE_TITLE,
            VERBOSE=VERBOSE, MODE=MODE)
        output_data(output, FILENAME=file2, DATAFORM=DATAFORM,
            UNITS=UNITS_TITLE, LONGNAME=LONGNAME, TITLE=FILE_TITLE,
            KEY='error', VERBOSE=VERBOSE, MODE=MODE)
        # add output files to list object
        output_files.append(file1)
        output_files.append(file2)

    # if fitting coefficients with cyclical components
    if (ncycles > 0):
        # output spatial titles for amplitudes
        amp_title = {'ANN':'Annual Amplitude','SEMI':'Semi-Annual Amplitude',
            'S2':'S2 Tidal Alias Amplitude'}
        ph_title = {'ANN':'Annual Phase','SEMI':'Semi-Annual Phase',
            'S2':'S2 Tidal Alias Phase'}

        # output amplitude and phase of cyclical components
        for i,flag in enumerate(amp_str):
            # Indice pointing to the cyclical components
            j = 1 + ORDER + 2*i
            # Allocating memory for output amplitude and phase
            amp = dinput.zeros_like()
            ph = dinput.zeros_like()
            # calculating amplitude and phase of spatial field
            amp.data,ph.data = gravtk.time_series.amplitude(
                out.data[:,:,j], out.data[:,:,j+1]
            )
            # convert phase from -180:180 to 0:360
            ii,jj = np.nonzero(ph.data < 0)
            ph.data[ii,jj] += 360.0
            # Amplitude Error
            comp1 = out.error[:,:,j]*out.data[:,:,j]/amp.data
            comp2 = out.error[:,:,j+1]*out.data[:,:,j+1]/amp.data
            amp.error = np.sqrt(comp1**2 + comp2**2)
            # Phase Error (degrees)
            comp1 = out.error[:,:,j]*out.data[:,:,j+1]/(amp.data**2)
            comp2 = out.error[:,:,j+1]*out.data[:,:,j]/(amp.data**2)
            ph.error = (180.0/np.pi)*np.sqrt(comp1**2 + comp2**2)

            # output file names for amplitude, phase and errors
            f3 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
                gw_str, ds_str, flag, '', START, END, suffix)
            f4 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
                gw_str, ds_str, flag,'_PHASE', START, END, suffix)
            file3 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f3))
            file4 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f4))
            # output spatial error file name
            f5 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
                gw_str, ds_str, flag, '_ERROR', START, END, suffix)
            f6 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
                gw_str, ds_str, flag, '_PHASE_ERROR', START, END, suffix)
            file5 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f5))
            file6 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f6))
            # full attributes
            AMP_UNITS = unit_list[UNITS-1]
            PH_UNITS = 'degrees'
            LONGNAME = unit_name[UNITS-1]
            AMP_TITLE = f'GRACE/GRACE-FO_Spatial_Data_{amp_title[flag]}'
            PH_TITLE = f'GRACE/GRACE-FO_Spatial_Data_{ph_title[flag]}'
            # Output seasonal amplitude and phase to files
            output_data(amp, FILENAME=file3, DATAFORM=DATAFORM,
                UNITS=AMP_UNITS, LONGNAME=LONGNAME, TITLE=AMP_TITLE,
                VERBOSE=VERBOSE, MODE=MODE)
            output_data(ph, FILENAME=file4, DATAFORM=DATAFORM,
                UNITS=PH_UNITS, LONGNAME='Phase', TITLE=PH_TITLE,
                VERBOSE=VERBOSE, MODE=MODE)
            # Output seasonal amplitude and phase error to files
            output_data(amp, FILENAME=file5, DATAFORM=DATAFORM,
                UNITS=AMP_UNITS, LONGNAME=LONGNAME, TITLE=AMP_TITLE,
                KEY='error', VERBOSE=VERBOSE, MODE=MODE)
            output_data(ph, FILENAME=file6, DATAFORM=DATAFORM,
                UNITS=PH_UNITS, LONGNAME='Phase', TITLE=PH_TITLE,
                KEY='error', VERBOSE=VERBOSE, MODE=MODE)
            # add output files to list object
            output_files.append(file3)
            output_files.append(file4)
            output_files.append(file5)
            output_files.append(file6)

    # Output fit significance
    signif_longname = {'SSE':'Sum of Squares Error',
        'AIC':'Akaike information criterion',
        'BIC':'Bayesian information criterion',
        'R2Adj':'Adjusted Coefficient of Determination'}
    # for each fit significance term
    for key,fs in FS.items():
        # output file names for fit significance
        signif_str = f'{key}_'
        f7 = (FILE_PREFIX, unit_list[UNITS-1], LMAX, order_str,
            gw_str, ds_str, signif_str, coef_str[ORDER], START, END, suffix)
        file7 = os.path.join(OUTPUT_DIRECTORY,output_format.format(*f7))
        # full attributes
        LONGNAME = signif_longname[key]
        # output fit significance to file
        output_data(fs, FILENAME=file7, DATAFORM=DATAFORM,
            UNITS=key, LONGNAME=LONGNAME, TITLE=nu,
            VERBOSE=VERBOSE, MODE=MODE)
        # add output files to list object
        output_files.append(os.path.join(OUTPUT_DIRECTORY,f7))

    # return the list of output files
    return output_files

# PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, KEY='data', DATAFORM=None,
    UNITS=None, LONGNAME=None, TITLE=None, VERBOSE=0, MODE=0o775):
    output = data.copy()
    setattr(output,'data',getattr(data, KEY))
    # attributes for output files
    attributes = {}
    attributes['units'] = UNITS
    attributes['longname'] = LONGNAME
    attributes['title'] = TITLE
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # write to output file
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        output.to_ascii(FILENAME, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        output.to_netCDF4(FILENAME, date=False, verbose=VERBOSE,
            **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        output.to_HDF5(FILENAME, date=False, verbose=VERBOSE,
            **attributes)
    # change the permissions mode of the output file
    os.chmod(FILENAME, MODE)

# PURPOSE: print a file log for the GRACE/GRACE-FO regression
def output_log_file(input_arguments, output_files):
    # format: GRACE_processing_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the GRACE/GRACE-FO regression
def output_error_log_file(input_arguments):
    # format: GRACE_processing_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads in GRACE/GRACE-FO spatial files and calculates the
            trends at each grid point following an input regression model
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for spatial files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for input and output files')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month for time series regression')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month for time series regression')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3,4,5],
        help='Output units')
    # output grid parameters
    parser.add_argument('--spacing',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    parser.add_argument('--redistribute-removed',
        default=False, action='store_true',
        help='Redistribute removed mass fields over the ocean')
    # regression parameters
    # 0: mean
    # 1: trend
    # 2: acceleration
    parser.add_argument('--order',
        type=int, default=2,
        help='Regression fit polynomial order')
    # regression fit cyclical terms
    parser.add_argument('--cycles',
        type=float, default=[0.5,1.0,161.0/365.25], nargs='+',
        help='Regression fit cyclical terms')
    # Output log file for each job in forms
    # GRACE_processing_run_2002-04-01_PID-00000.log
    # GRACE_processing_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run regress_grace_maps algorithm with parameters
        output_files = regress_grace_maps(
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            DATAFORM=args.format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            ORDER=args.order,
            CYCLES=args.cycles,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as e:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
