#!/usr/bin/env python
u"""
grace_spatial_mean.py
Written by Tyler Sutterley (05/2024)

Calculates the mean GRACE/GRACE-FO trend for a range of GIA outputs

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -f X, --flag X: GRACE/GRACE-FO specific data flag
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
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
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        AW13-IJ05: Geruo A IJ05-R2 GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA files to read
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
    Updated 05/2024: fixed filenames for output masked files
    Updated 06/2023: fixed output unit attributes for each variable
    Updated 05/2023: use pathlib to define and operate on paths
        read from merged data and error file as output from regression program
    Updated 03/2023: updated inputs to spatial from_file function
        use attributes from units class for writing to netCDF4/HDF5 files
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: added AW13 combination IJ05-R2 and ICE6G models
    Updated 01/2022: added F-test statistics from regression fit terms
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        simplified file imports and exports using wrappers in spatial utility
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 07/2020: rewrite of program to use spatial class for operations
    Updated 10/2019: get parameters with getopt
    Updated 12/2018: python3 compatibility updates
    Updated 06/2018: output masked files using method from Velicogna (2014)
    Updated 09/2017: simplified outputs variables for netCDF4 files
    Updated 08/2017: split SLF maps into calculate_mean_SLF_maps.py
        set grid dimensions, calculate mean of ICE-6G VM5 GIA models
    Written 08/2017
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

# GIA dictionaries
gia_mean_str = {}
# Ivins et al. (2013)
gia_mean_str['IJ05-R2'] = '_IJ05_R2_mean'
# Whitehouse et al. (2012)
gia_mean_str['W12a'] = '_W12a_mean'
# Simpson et al. (2009)
gia_mean_str['SM09'] = '_SM09_HUY2_mean'
# Peltier et al. (2015)
gia_mean_str['ICE6G'] = '_ICE6G_VM5_mean'
# Peltier et al. (2018)
gia_mean_str['ICE6G-D'] = '_ICE6G-D_mean'
# A et al. (2013)
gia_mean_str['AW13-ICE6G'] = '_AW13_ICE6G_mean'
gia_mean_str['ascii'] = '_AW13_IJ05_ICE6G_mean'

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate the GRACE/GRACE-FO mean maps for a set of GIA models
def grace_spatial_mean(PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MMAX=None,
    DESTRIPE=False,
    UNITS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    GIA=None,
    GIA_FILES=None,
    DATAFORM=None,
    REDISTRIBUTE_REMOVED=False,
    OUTPUT_DIRECTORY=None,
    FLAG=None,
    VERBOSE=0,
    MODE=0o775):

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
    ff='{0}_{1}_{2}{3}_{4}_{5}_L{6:d}{7}{8}{9}_{10}_{11:03d}-{12:03d}.{13}'

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
    dinput['x1'] = []
    dinput['x2'] = []
    # fit significance terms
    dinput['AIC_x0'] = []
    dinput['AIC_x1'] = []
    dinput['AIC_x2'] = []
    dinput['SSE_x0'] = []
    dinput['SSE_x1'] = []
    dinput['SSE_x2'] = []
    # degrees of freedom for each regression
    DOF = dict(x0=0, x1=0, x2=0)
    # number of rheologies and ice histories to run
    N = len(GIA_FILES)
    # iterate GIA models
    for h,GIA_FILE in enumerate(GIA_FILES):
        # input GIA spherical harmonic datafiles
        GIA_FILE = pathlib.Path(GIA_FILE).expanduser().absolute() if GIA else None
        GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
        gia_str = f'_{GIA_Ylms_rate.title}'
        # read mass trend file and mass trend error files
        for key in ['x1','x2']:
            F1 = ff.format(PROC,DREL,DSET,gia_str,FLAG,units,LMAX,
                order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
            INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
            field_mapping = dict(lon='lon', lat='lat', data='data', error='error')
            temp = gravtk.spatial().from_file(INPUT_FILE,
                format=DATAFORM, date=False, spacing=[dlon,dlat],
                nlon=nlon, nlat=nlat, field_mapping=field_mapping)
            dinput[key].append(temp)
        # read AIC files
        for key in ['AIC_x0','AIC_x1','AIC_x2']:
            F1 = ff.format(PROC,DREL,DSET,gia_str,FLAG,units,LMAX,
                order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
            INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
            field_mapping = dict(lon='lon', lat='lat', data='data')
            temp = gravtk.spatial().from_file(INPUT_FILE,
                format=DATAFORM, date=False, field_mapping=field_mapping,
                spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
            dinput[key].append(temp)
        # read SSE files
        for key in ['SSE_x0','SSE_x1','SSE_x2']:
            F1 = ff.format(PROC,DREL,DSET,gia_str,FLAG,units,LMAX,
                order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
            INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
            field_mapping = dict(lon='lon', lat='lat', data='data')
            temp = gravtk.spatial().from_file(INPUT_FILE,
                format=DATAFORM, date=False, field_mapping=field_mapping,
                spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
            dinput[key].append(temp)
            # save degrees of freedom for fit order
            order = key.replace('SSE_','')
            DOF[order] = int(temp.attributes['ROOT']['title'])

    # create combined spatial objects
    output = {}
    output_units = {}
    # calculate mean GIA-corrected x1 change (for all Earth rheologies)
    x1 = gravtk.spatial().from_list(dinput['x1'],date=False)
    output['x1'] = x1.mean()
    output_units['x1'] = '{0} yr^{1:d}'.format(units_name,-1)
    # calculate mean acceleration x2 change (for all Earth rheologies)
    x2 = gravtk.spatial().from_list(dinput['x2'],date=False)
    output['x2'] = x2.mean().scale(2.0)
    output_units['x2'] = '{0} yr^{1:d}'.format(units_name,-2)
    # GRACE satellite error component
    e1 = x1.copy(); e1.data = np.copy(x1.error)
    output['x1'].error = e1.sum(power=2.0).scale(1.0/N).power(0.5).data
    e2 = x2.copy(); e2.data = np.copy(x2.error)
    output['x2'].error = e2.sum(power=2.0).scale(4.0/N).power(0.5).data
    # significance means
    AICx0 = gravtk.spatial().from_list(dinput['AIC_x0'],date=False).mean()
    AICx1 = gravtk.spatial().from_list(dinput['AIC_x1'],date=False).mean()
    AICx2 = gravtk.spatial().from_list(dinput['AIC_x2'],date=False).mean()
    # calculate residual sum of squares means
    s0 = gravtk.spatial().from_list(dinput['SSE_x0'],date=False)
    SSEx0 = s0.sum(power=2.0).scale(1.0/N).power(0.5)
    s1 = gravtk.spatial().from_list(dinput['SSE_x1'],date=False)
    SSEx1 = s1.sum(power=2.0).scale(1.0/N).power(0.5)
    s2 = gravtk.spatial().from_list(dinput['SSE_x2'],date=False)
    SSEx2 = s2.sum(power=2.0).scale(1.0/N).power(0.5)

    # invalid value for masked grids
    fill_value = -9999.0
    # masked trend values
    ii,jj = np.nonzero((np.abs(output['x1'].data) <= output['x1'].error) |
        (AICx1.data >= AICx0.data))
    output['MASKED_x1'] = output['x1'].copy()
    output['MASKED_x1'].fill_value = fill_value
    output['MASKED_x1'].mask[ii,jj] = True
    output['MASKED_x1'].update_mask()
    output_units['MASKED_x1'] = '{0} yr^{1:d}'.format(units_name,-1)
    # write masked acceleration values to file
    ii,jj = np.nonzero((np.abs(output['x2'].data) <= output['x2'].error) |
        (AICx2.data >= AICx1.data))
    output['MASKED_x2'] = output['x2'].copy()
    output['MASKED_x2'].fill_value = fill_value
    output['MASKED_x2'].mask[ii,jj] = True
    output['MASKED_x2'].update_mask()
    output_units['MASKED_x2'] = '{0} yr^{1:d}'.format(units_name,-2)

    # attributes for output files
    title = 'GRACE/GRACE-FO Spatial Data'
    # output data to file
    for key,val in output.items():
        F2 = ff.format(PROC,DREL,DSET,gia_mean_str[GIA],FLAG,units,
            LMAX,order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(F2)
        output_data(val, FILENAME=OUTPUT_FILE, DATAFORM=DATAFORM,
            UNITS=output_units[key], LONGNAME=units_longname,
            TITLE=title, CONF=0.95, VERBOSE=VERBOSE, MODE=MODE)
        # add to output file list
        output_files.append(OUTPUT_FILE)

    # Calculate F-statistics
    output['F_x1'] = SSEx1.zeros_like()
    output['F_x1'].data = (SSEx0.data - SSEx1.data) / \
        (DOF['x0'] - DOF['x1'])*(SSEx1.data / DOF['x1'])
    output['F_x2'] = SSEx2.zeros_like()
    output['F_x2'].data  = (SSEx1.data - SSEx2.data) / \
        (DOF['x1'] - DOF['x2'])*(SSEx2.data / DOF['x2'])
    # calculate F-test p-values
    output['P_x1'] = SSEx1.zeros_like()
    output['P_x1'].data = 1.0 - scipy.stats.f.cdf(output['F_x1'].data,
        DOF['x0'] - DOF['x1'], DOF['x1'])
    output['P_x2'] = SSEx2.zeros_like()
    output['P_x2'].data = 1.0 - scipy.stats.f.cdf(output['F_x2'].data,
        DOF['x1'] - DOF['x2'], DOF['x2'])

    # for each F-test significance term
    signif_longname = {'F_x1':'F Statistic','F_x2':'F Statistic',
        'P_x1':'F-test P-value','P_x2':'F-test P-value'}
    # output data to file
    for key in ['F_x1', 'F_x2', 'P_x1', 'P_x2']:
        # F-test significance term
        # output file names for fit significance
        F3 = ff.format(PROC,DREL,DSET,gia_mean_str[GIA],FLAG,units,
            LMAX,order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
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
        description="""Calculates the mean GRACE/GRACE-FO trend for a range
            of GIA outputs
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for spatial files')
    # GRACE/GRACE-FO flags
    # 'wSLR_C20_wDEG1'
    # 'rmTWC'
    # 'rmSLF_rmTWC'
    parser.add_argument('--flag','-f',
        type=str, default='rmSLF_rmTWC',
        help='GRACE/GRACE-FO specific data flags')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
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
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        nargs='+', help='GIA files to read')
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
        # run grace_spatial_mean algorithm with parameters
        output_files = grace_spatial_mean(
            args.center,
            args.release,
            args.product,
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
            GIA=args.gia,
            GIA_FILES=args.gia_file,
            DATAFORM=args.format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            OUTPUT_DIRECTORY=args.output_directory,
            FLAG=args.flag,
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
