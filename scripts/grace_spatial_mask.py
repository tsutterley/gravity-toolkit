#!/usr/bin/env python
u"""
grace_spatial_mask.py
Written by Tyler Sutterley (11/2024)

Creates mask files for GRACE/GRACE-FO trends following Velicogna (2014)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -P X, --file-prefix X: prefix string for input and output files
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
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
    --redistribute-removed: redistribute removed mass fields over the ocean
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

UPDATE HISTORY:
    Written 11/2024
"""
import sys
import os
import copy
import logging
import pathlib
import argparse
import traceback
import numpy as np
import scipy.stats
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate the GRACE/GRACE-FO mask files
def grace_spatial_mask(LMAX, RAD,
        START=None,
        END=None,
        MMAX=None,
        DESTRIPE=False,
        UNITS=None,
        DDEG=None,
        INTERVAL=None,
        BOUNDS=None,
        DATAFORM=None,
        REDISTRIBUTE_REMOVED=False,
        OUTPUT_DIRECTORY=None,
        FILE_PREFIX=None,
        VERBOSE=0,
        MODE=0o775
    ):

    # recursively create output directory if not currently existing
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # list object of output files for file logs (full path)
    output_files = []

    # file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # flag for spherical harmonic order
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # Calculating the Gaussian smoothing for radius RAD
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # destriped GRACE/GRACE-FO coefficients
    ds_str = '_FL' if DESTRIPE else ''
    # distributing removed mass uniformly over ocean
    ocean_str = '_OCN' if REDISTRIBUTE_REMOVED else ''
    # input and output spatial units
    # 1: cmwe, centimeters water equivalent
    # 2: mmGH, millimeters geoid height
    # 3: mmCU, millimeters elastic crustal deformation
    # 4: micGal, microGal gravity perturbations
    # 5: mbar, millibars equivalent surface pressure
    units = gravtk.units.bycode(UNITS)
    units_name, units_longname = gravtk.units.get_attributes(units)

    # input GRACE file format
    output_format = '{0}{1}_L{2:d}{3}{4}{5}_{6}_{7:03d}-{8:03d}.{9}'

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
        nlon = np.int64((maxlon-minlon)/dlon)
        nlat = np.int64((maxlat-minlat)/dlat)

    # allocate lists for input variables
    dinput = {}
    # degrees of freedom for each regression
    DOF = dict(x0=0, x1=0, x2=0)

    # read mass trend file and mass trend error files
    for key in ['x1','x2']:
        F1 = output_format.format(FILE_PREFIX, units, LMAX, order_str,
            gw_str, ds_str, key, START, END, suffix[DATAFORM])
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
        field_mapping = dict(lon='lon', lat='lat', data='data', error='error')
        dinput[key] = gravtk.spatial().from_file(INPUT_FILE,
            format=DATAFORM, date=False, spacing=[dlon,dlat],
            nlon=nlon, nlat=nlat, field_mapping=field_mapping)
    # read AIC files
    for key in ['AIC_x0','AIC_x1','AIC_x2']:
        F1 = output_format.format(FILE_PREFIX, units, LMAX, order_str,
            gw_str, ds_str, key, START, END, suffix[DATAFORM])
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
        field_mapping = dict(lon='lon', lat='lat', data='data')
        dinput[key] = gravtk.spatial().from_file(INPUT_FILE,
            format=DATAFORM, date=False, field_mapping=field_mapping,
            spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
    # read SSE files
    for key in ['SSE_x0','SSE_x1','SSE_x2']:
        F1 = output_format.format(FILE_PREFIX, units, LMAX, order_str,
            gw_str, ds_str, key, START, END, suffix[DATAFORM])
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
        field_mapping = dict(lon='lon', lat='lat', data='data')
        dinput[key] = gravtk.spatial().from_file(INPUT_FILE,
            format=DATAFORM, date=False, field_mapping=field_mapping,
            spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
        # save degrees of freedom for fit order
        order = key.replace('SSE_','')
        DOF[order] = int(dinput[key].attributes['ROOT']['title'])

    # output mask files
    output = {}
    output_units = {}
    # invalid value for masked grids
    fill_value = -9999.0
    # masked trend values
    ii,jj = np.nonzero((np.abs(dinput['x1'].data) <= dinput['x1'].error) |
        (dinput['AIC_x1'].data >= dinput['AIC_x0'].data))
    output['MASKED_x1'] = dinput['x1'].copy()
    output['MASKED_x1'].fill_value = fill_value
    output['MASKED_x1'].mask[ii,jj] = True
    output['MASKED_x1'].update_mask()
    output_units['MASKED_x1'] = '{0} yr^{1:d}'.format(units_name,-1)
    # write masked acceleration values to file
    ii,jj = np.nonzero((np.abs(dinput['x2'].data) <= dinput['x2'].error) |
        (dinput['AIC_x2'].data >= dinput['AIC_x1'].data))
    output['MASKED_x2'] = dinput['x2'].scale(2.0)
    output['MASKED_x2'].fill_value = fill_value
    output['MASKED_x2'].mask[ii,jj] = True
    output['MASKED_x2'].update_mask()
    output_units['MASKED_x2'] = '{0} yr^{1:d}'.format(units_name,-2)

    # attributes for output files
    title = 'GRACE/GRACE-FO Spatial Data'
    # output data to file
    for key,val in output.items():
        F2 = output_format.format(FILE_PREFIX, units, LMAX, order_str,
            gw_str, ds_str, key, START, END, suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(F2)
        output_data(val, FILENAME=OUTPUT_FILE, DATAFORM=DATAFORM,
            UNITS=output_units[key], LONGNAME=units_longname,
            TITLE=title, CONF=0.95, VERBOSE=VERBOSE, MODE=MODE)
        # add to output file list
        output_files.append(OUTPUT_FILE)

    # Calculate F-statistics
    output['F_x1'] = dinput['SSE_x1'].zeros_like()
    output['F_x1'].data = (dinput['SSE_x0'].data - dinput['SSE_x1'].data) / \
        (DOF['x0'] - DOF['x1'])*(dinput['SSE_x1'].data / DOF['x1'])
    output['F_x2'] = dinput['SSE_x2'].zeros_like()
    output['F_x2'].data  = (dinput['SSE_x1'].data - dinput['SSE_x2'].data) / \
        (DOF['x1'] - DOF['x2'])*(dinput['SSE_x2'].data / DOF['x2'])
    # calculate F-test p-values
    output['P_x1'] = dinput['SSE_x1'].zeros_like()
    output['P_x1'].data = 1.0 - scipy.stats.f.cdf(output['F_x1'].data,
        DOF['x0'] - DOF['x1'], DOF['x1'])
    output['P_x2'] = dinput['SSE_x2'].zeros_like()
    output['P_x2'].data = 1.0 - scipy.stats.f.cdf(output['F_x2'].data,
        DOF['x1'] - DOF['x2'], DOF['x2'])

    # for each F-test significance term
    signif_longname = {'F_x1':'F Statistic','F_x2':'F Statistic',
        'P_x1':'F-test P-value','P_x2':'F-test P-value'}
    # output data to file
    for key in ['F_x1', 'F_x2', 'P_x1', 'P_x2']:
        # F-test significance term
        # output file names for fit significance
        F3 = output_format.format(FILE_PREFIX, units, LMAX, order_str,
            gw_str, ds_str, key, START, END, suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(F3)
        output_data(output[key], FILENAME=OUTPUT_FILE,
            DATAFORM=DATAFORM, UNITS='1', LONGNAME=signif_longname[key],
            TITLE=title, VERBOSE=VERBOSE, MODE=MODE)
        # add to output file list
        output_files.append(OUTPUT_FILE)

    # return the list of files
    return output_files

# PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, DATAFORM=None, UNITS=None,
    LONGNAME=None, TITLE=None, CONF=0, VERBOSE=0, MODE=0o775):
    # field mapping for output regression data
    field_mapping = {}
    field_mapping['lat'] = 'lat'
    field_mapping['lon'] = 'lon'
    field_mapping['data'] = 'data'
    # add data error to output field mapping
    if hasattr(data, 'error'):
        field_mapping['error'] = 'error'
    # output attributes
    attributes = dict(ROOT={})
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    attributes['data'] = {}
    attributes['data']['description'] = 'Model fit'
    attributes['data']['long_name'] = LONGNAME
    attributes['data']['units'] = UNITS
    attributes['error'] = {}
    attributes['error']['description'] = 'Uncertainty in model fit'
    attributes['error']['long_name'] = LONGNAME
    attributes['error']['units'] = UNITS
    attributes['error']['confidence'] = 100*CONF
    # output global attributes
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # write to output file
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        data.to_ascii(FILENAME, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        data.to_netCDF4(FILENAME, date=False, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attributes,
            title=TITLE, reference=REFERENCE)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        data.to_HDF5(FILENAME, date=False, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attributes,
            title=TITLE, reference=REFERENCE)
    # change the permissions mode of the output file
    FILENAME.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates mask files for GRACE/GRACE-FO trends
            following Velicogna (2014)
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
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
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
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
        # run grace_spatial_mask algorithm with parameters
        output_files = grace_spatial_mask(
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            DATAFORM=args.format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
