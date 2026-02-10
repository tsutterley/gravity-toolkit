#!/usr/bin/env python
u"""
geocenter_spatial_maps.py
Written by Tyler Sutterley (05/2023)

Reads in GRACE/GRACE-FO geocenter coefficients and exports trends
    in the monthly spatial fields in millimeters water equivalent

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for spatial files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    -k X, --kl X: Degree 1 Gravitational Load Love number
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: output data format
        ascii
        netCDF4
        HDF5
    --slr-c20 X: Replace C20 coefficients with SLR values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GFZ: use values from GFZ
        GSFC: use values from GSFC (TN-14)
    --slr-21 X: Replace C21 and S21 coefficients with SLR values
        CSR: use values from CSR
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC
    --slr-22 X: Replace C22 and S22 coefficients with SLR values
        CSR: use values from CSR
    --slr-c30 X: Replace C30 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC (TN-14)
    --slr-c40 X: Replace C40 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --slr-c50 X: Replace C50 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
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
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces low-degree harmonics with SLR values (if specified)
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    geocenter.py: converts between spherical harmonics and geocenter variations
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    spatial.py: spatial data class for reading, writing and processing data
    time_series.regress.py: calculates trend coefficients using least-squares
    time_series.amplitude.py: calculates the amplitude and phase of a harmonic
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: split S2 tidal aliasing terms into GRACE and GRACE-FO eras
        use pathlib to define and operate on paths
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: option to use solutions with SLR degree 4 zonal harmonics
    Updated 07/2022: create mask for output gridded variables
    Updated 05/2022: use argparse descriptions within documentation
        use command line option to set degree 1 gravitational love number
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: add GSFC low-degree harmonics
        use gravity_toolkit geocenter class for operations
        can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
        add more choices for setting input format of the removed files
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 12/2020: added more love number options and from gfc for mean files
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for output operations
    Updated 10/2019: changing Y/N flags to True/False
        calculate phase of the regressed harmonics
    Updated 01/2019: calculate and output annual and semi-annual phase
        include 161-day S2 tidal aliasing terms in regression
    Written 10/2018
"""
from __future__ import print_function

import sys
import os
import copy
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: import GRACE/GRACE-FO geocenter files for a given months range
# Converts the GRACE/GRACE-FO harmonics applying the specified procedures
def geocenter_spatial_maps(base_dir, PROC, DREL,
    START=None,
    END=None,
    MISSING=None,
    LOVE_NUMBERS=0,
    LOVE_K1=None,
    REFERENCE=None,
    DESTRIPE=False,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C40=None,
    SLR_C50=None,
    DATAFORM=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # input directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    grace_dir = base_dir.joinpath('geocenter')
    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # list object of output files for file logs (full path)
    output_files = []

    # file information
    file_format = '{0}_{1}_{2}{3}{4}{5}_{6}_{7}{8}.{9}'
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # GRACE months
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    months = sorted(set(np.arange(START,END+1)) - set(MISSING))
    # labels for Release-6
    model_str = 'OMCT' if DREL in ('RL04','RL05') else 'MPIOM'
    # labels for each scenario
    FLAGS = []
    # FLAGS.append('')
    # FLAGS.append('_iter')
    FLAGS.append('_SLF_iter')
    # GIA and processing labels
    gia_str = '_AW13_ice6g_GA'
    ds_str = '_FL' if DESTRIPE else ''
    # output flag for low-degree harmonic replacements
    if SLR_21 in ('CSR','GFZ','GSFC'):
        C21_str = f'_w{SLR_21}_21'
    else:
        C21_str = ''
    if SLR_22 in ('CSR','GSFC'):
        C22_str = f'_w{SLR_22}_22'
    else:
        C22_str = ''
    if SLR_C30 in ('GSFC',):
        # C30 replacement now default for all solutions
        C30_str = ''
    elif SLR_C30 in ('CSR','GFZ','LARES'):
        C30_str = f'_w{SLR_C30}_C30'
    else:
        C30_str = ''
    if SLR_C40 in ('CSR','GSFC','LARES'):
        C40_str = f'_w{SLR_C40}_C40'
    else:
        C40_str = ''
    if SLR_C50 in ('CSR','GSFC','LARES'):
        C50_str = f'_w{SLR_C50}_C50'
    else:
        C50_str = ''
    # combine satellite laser ranging flags
    slr_str = ''.join([C21_str,C22_str,C30_str,C40_str,C50_str])
    # degree one coefficient labels
    coef_labels = ['C10','C11','S11']

    # read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(1, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)
    # set gravitational load love number to a specific value
    if LOVE_K1:
        kl[1] = np.copy(LOVE_K1)

    # mmwe, millimeters water equivalent
    unit_label = 'mmwe'
    unit_name = 'Equivalent_Water_Thickness'
    dfactor = gravtk.units(lmax=1).harmonic(hl,kl,ll).mmwe
    # attributes for output files
    attributes = {}
    attributes['field_mapping'] = dict(lon='lon', lat='lat', data='z', time='time')
    attributes['time_units'] = 'years'
    attributes['time_longname'] = 'Date_in_Decimal_Years'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # read the GAE/GAF/GAG atmospheric correction coefficient files
    atm_corr = gravtk.read_ecmwf_corrections(base_dir, 1, months)

    # Output spatial data object
    grid = gravtk.spatial()
    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        grid.lon = -180 + dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        grid.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0, maxlon+dlon/2.0, dlon)
        grid.lat = np.arange(maxlat-dlat/2.0, minlat-dlat/2.0, -dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)

    # Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(2, np.cos(theta))
    # fit coefficients
    fits = ['x1','x2','SS','SC','AS','AC','S2SGRC','S2CGRC','S2SGFO','S2CGFO']

    # for each test run
    for i,input_flag in enumerate(FLAGS):
        # read geocenter file for processing center and model
        fargs = (PROC,DREL,model_str,input_flag,slr_str,gia_str,ds_str)
        grace_file = '{0}_{1}_{2}{3}{4}{5}{6}.txt'.format(*fargs)
        DEG1 = gravtk.geocenter().from_UCI(grace_dir.joinpath(grace_file))
        # indices for months
        kk,=np.nonzero((DEG1.month >= START) & (DEG1.month <= END))
        MEAN = DEG1.mean(indices=kk)
        # if data is Release-5: remove ECMWF jump corrections
        if (DREL == 'RL05'):
            DEG1.C10[kk] -= atm_corr['clm'][1,0,:]
            DEG1.C11[kk] -= atm_corr['clm'][1,1,:]
            DEG1.S11[kk] -= atm_corr['slm'][1,1,:]
        # create dictionary for extracting regressed coefficients
        Ylm = {}
        for k,f in enumerate(fits):
            Ylm[f] = np.zeros((3))
        # calculate regression over each coefficient
        for j,key in enumerate(coef_labels):
            val = getattr(DEG1, key)
            # calculate regression coefficients
            TERMS = gravtk.time_series.aliasing_terms(DEG1.time[kk])
            x1 = gravtk.time_series.regress(DEG1.time[kk], dfactor[1]*val[kk],
                ORDER=1, CYCLES=[0.5,1.0], TERMS=TERMS,
                CONF=0.95, AICc=True)
            x2 = gravtk.time_series.regress(DEG1.time[kk], dfactor[1]*val[kk],
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS,
                CONF=0.95, AICc=True)
            # save coefficients
            for k,f in enumerate(fits):
                Ylm[f][j] = x1['beta'][k+1]
        # extract Ylms and convert to spatial
        var = {}
        for k,f in enumerate(fits):
            clm = np.zeros((2,2))
            slm = np.zeros((2,2))
            clm[1,0] = Ylm[f][0]
            clm[1,1] = Ylm[f][1]
            slm[1,1] = Ylm[f][2]
            var[f] = gravtk.harmonic_summation(clm, slm, grid.lon, grid.lat,
                LMAX=1, MMAX=1, PLM=PLM[:2,:2,:]).T
        # amplitude and phase of cyclical components
        var['SA'],var['SP'] = gravtk.time_series.amplitude(var['SS'],var['SC'])
        var['AA'],var['AP'] = gravtk.time_series.amplitude(var['AS'],var['AC'])
        var['S2AGRC'],var['S2PGRC'] = gravtk.time_series.amplitude(var['S2SGRC'],var['S2CGRC'])
        var['S2AGFO'],var['S2PGFO'] = gravtk.time_series.amplitude(var['S2SGFO'],var['S2CGFO'])

        # out regression coefficients and amplitudes to file
        unit_suffix = [' yr^-1', ' yr^-2', '', '', '']
        for j,key in enumerate(['x1','x2','SA','AA','S2AGRC','S2AGFO']):
            # copy variables to output grid
            grid.data = np.copy(var[key])
            grid.mask = np.zeros_like(grid.data, dtype=bool)
            grid.time = np.copy(MEAN.time)
            # add specific attributes
            attributes['units'] = f'{unit_label}{unit_suffix[j]}'
            attributes['longname'] = copy.copy(unit_name)
            attributes['title'] = copy.copy(key)
            # save to file
            FILE = file_format.format(PROC,DREL,model_str,input_flag,
                slr_str,gia_str,unit_label,key,ds_str,suffix[DATAFORM])
            OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
            grid.to_file(OUTPUT_FILE, format=DATAFORM, **attributes)
            # add file to output list
            output_files.append(OUTPUT_FILE)
            # change the permissions mode
            OUTPUT_FILE.chmod(mode=MODE)

        # output phase to file
        for key in ['SP','AP','S2PGRC','S2PGFO']:
            # convert phase from -180:180 to 0:360
            ix,iy = np.nonzero(var[key] < 0)
            var[key][ix,iy] += 360.0
            # copy variables to output grid
            grid.data = np.copy(var[key])
            grid.time = np.copy(MEAN.time)
            # add specific attributes
            attributes['units'] = 'degrees'
            attributes['longname'] = 'Phase'
            attributes['title'] = copy.copy(key)
            # save to file
            FILE = file_format.format(PROC,DREL,model_str,input_flag,
                slr_str,gia_str,unit_label,key,ds_str,suffix[DATAFORM])
            OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
            grid.to_file(OUTPUT_FILE, format=DATAFORM, **attributes)
            # add file to output list
            output_files.append(OUTPUT_FILE)
            # change the permissions mode
            OUTPUT_FILE.chmod(mode=MODE)
    # return the list of output files
    return output_files

# PURPOSE: print a file log for the geocenter analysis
def output_log_file(input_arguments, output_files):
    # format: geocenter_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'geocenter_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
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

# PURPOSE: print a error file log for the geocenter analysis
def output_error_log_file(input_arguments):
    # format: geocenter_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'geocenter_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
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
        description="""Reads in GRACE/GRACE-FO geocenter coefficients
            and exports trends in the monthly spatial fields in
            millimeters water equivalent
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for spatial files')
    # Data processing center or satellite mission
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    parser.add_argument('--kl','-k',
        type=float, default=0.021,
        help='Degree 1 gravitational Load Love number')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
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
    # replace low degree harmonics with values from Satellite Laser Ranging
    parser.add_argument('--slr-c20',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C20 coefficients with SLR values')
    parser.add_argument('--slr-21',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C21 and S21 coefficients with SLR values')
    parser.add_argument('--slr-22',
        type=str, default=None, choices=['CSR','GSFC'],
        help='Replace C22 and S22 coefficients with SLR values')
    parser.add_argument('--slr-c30',
        type=str, default=None, choices=['CSR','GFZ','GSFC','LARES'],
        help='Replace C30 coefficients with SLR values')
    parser.add_argument('--slr-c40',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C40 coefficients with SLR values')
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    # Output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Output data format')
    # Output log file for each job in forms
    # geocenter_run_2002-04-01_PID-00000.log
    # geocenter_failed_run_2002-04-01_PID-00000.log
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
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run geocenter_spatial_maps algorithm with parameters
        output_files = geocenter_spatial_maps(
            args.directory,
            args.center,
            args.release,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            LOVE_NUMBERS=args.love,
            LOVE_K1=args.kl,
            REFERENCE=args.reference,
            DESTRIPE=args.destripe,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C40=args.slr_c40,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
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
