#!/usr/bin/env python
u"""
grace_spatial_differences.py
Written by Tyler Sutterley (05/2023)

Calculates the impact of sea level fingerprints on the GRACE/GRACE-FO trend

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -f X, --flag X: GRACE/GRACE-FO specific data flags
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
    Updated 06/2023: use field mapping to use specific output variable names
    Updated 05/2023: use pathlib to define and operate on paths
        read from merged data and error file as output from regression program
    Updated 03/2023: updated inputs to spatial from_file function
        use attributes from units class for writing to netCDF4/HDF5 files
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: added AW13 combination IJ05-R2 and ICE6G models
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        simplified file imports and exports using wrappers in spatial utility
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 07/2020: rewrite of program to use spatial class for operations
    Updated 10/2019: added processing center and can get parameters with getopt
    Updated 02/2019: added option for setting GRACE data release
    Updated 12/2018: python3 compatibility updates
    Updated 10/2018: output smoothed correction file
    Written 09/2018
"""
import sys
import os
import copy
import logging
import pathlib
import argparse
import traceback
import numpy as np
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

# PURPOSE: calculate the GRACE/GRACE-FO difference map for each GIA model
def grace_spatial_differences(PROC, DREL, DSET, LMAX, RAD,
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
    FLAGS=None,
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
    x1 = [None] * 2
    e1 = [None] * 2
    AICx0 = [None] * 2
    AICx1 = [None] * 2
    AICx2 = [None] * 2
    # number of rheologies and ice histories to run
    N = len(GIA_FILES)
    # iterate over flags
    for i,FLAG in enumerate(FLAGS):
        # allocate lists for input variables
        dinput = {}
        dinput['x1'] = []
        dinput['AIC_x0'] = []
        dinput['AIC_x1'] = []
        dinput['AIC_x2'] = []
        # iterate GIA models
        for GIA_FILE in GIA_FILES:
            # input GIA spherical harmonic datafiles
            GIA_FILE = pathlib.Path(GIA_FILE).expanduser().absolute() if GIA else None
            GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
            gia_str = f'_{GIA_Ylms_rate.title}'
            # read mass trend file and mass trend error files
            for key in ['x1']:
                F1 = ff.format(PROC,DREL,DSET,gia_str,FLAG,units,
                    LMAX,order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
                INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
                field_mapping = dict(lon='lon', lat='lat', data='data', error='error')
                temp = gravtk.spatial().from_file(INPUT_FILE,
                    format=DATAFORM, date=False, spacing=[dlon,dlat],
                    nlon=nlon, nlat=nlat, field_mapping=field_mapping)
                dinput[key].append(temp)
            # read AIC files
            for key in ['AIC_x0','AIC_x1','AIC_x2']:
                F1 = ff.format(PROC,DREL,DSET,gia_str,FLAG,units,
                    LMAX,order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
                INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
                field_mapping = dict(lon='lon', lat='lat', data='data')
                temp = gravtk.spatial().from_file(INPUT_FILE,
                    format=DATAFORM, date=False, field_mapping=field_mapping,
                    spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
                dinput[key].append(temp)
        # create combined spatial objects
        combined = gravtk.spatial().from_list(dinput['x1'],date=False)
        # calculate mean GIA-corrected x1 change (for all Earth rheologies)
        x1[i] = combined.mean()
        # GRACE satellite error component
        error = combined.copy(); error.data = np.copy(combined.error)
        e1[i] = error.sum(power=2.0).scale(1.0/N).power(0.5)
        # significance means
        AICx0[i] = gravtk.spatial().from_list(dinput['AIC_x0'],date=False).mean()
        AICx1[i] = gravtk.spatial().from_list(dinput['AIC_x1'],date=False).mean()
        AICx2[i] = gravtk.spatial().from_list(dinput['AIC_x2'],date=False).mean()

    # set invalid values for masked variables
    ii,jj = np.nonzero((np.abs(x1[0].data) <= e1[0].data) |
        (np.abs(x1[1].data) <= e1[1].data) |
        (AICx1[0].data >= AICx0[0].data) |
        (AICx1[1].data >= AICx0[1].data))

    # output difference variables
    output = {}
    output_units = {}
    output_longname = {}
    # calculate difference of corrected and uncorrected
    output['CORR'] = x1[0].zeros_like()
    output['CORR'].data = x1[0].data - x1[1].data
    output_units['CORR'] = '{0} yr^{1:d}'.format(units_name,-1)
    output_longname['CORR'] = units_longname
    # masked correction difference
    output['MASKED_CORR'] = output['CORR'].copy()
    output['MASKED_CORR'].mask[ii,jj] = True
    output['MASKED_CORR'].update_mask()
    output_units['MASKED_CORR'] = '{0} yr^{1:d}'.format(units_name,-1)
    output_longname['MASKED_CORR'] = units_longname

    # calculate RMS difference of corrected and uncorrected
    output['RMS'] = output['CORR'].power(2.0).power(0.5)
    output_units['RMS'] = '{0} yr^{1:d}'.format(units_name,-1)
    output_longname['RMS'] = units_longname
    # masked RMS difference
    output['MASKED_RMS'] = output['RMS'].copy()
    output['MASKED_RMS'].mask[ii,jj] = True
    output['MASKED_RMS'].update_mask()
    output_units['MASKED_RMS'] = '{0} yr^{1:d}'.format(units_name,-1)
    output_longname['MASKED_RMS'] = units_longname

    # calculate percent difference of corrected and uncorrected
    output['DIFF'] = x1[0].zeros_like()
    output['DIFF'].data = 100.0*output['RMS'].data/np.abs(x1[0].data)
    output_units['DIFF'] = '%'
    output_longname['DIFF'] = 'Percent'
    # masked percent difference
    output['MASKED_DIFF'] = output['DIFF'].copy()
    output['MASKED_DIFF'].mask[ii,jj] = True
    output['MASKED_DIFF'].update_mask()
    output_units['MASKED_DIFF'] = '%'
    output_longname['MASKED_DIFF'] = 'Percent'

    # field mapping for output regression data
    field_mapping = {}
    field_mapping['lat'] = 'lat'
    field_mapping['lon'] = 'lon'
    field_mapping['data'] = 'data'
    # attributes for output files
    attributes = {}
    attributes['title'] = 'GRACE/GRACE-FO Spatial Data'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # output data to file
    for key,val in output.items():
        # add specific data attributes
        attributes['units'] = copy.copy(output_units[key])
        attributes['longname'] = copy.copy(output_longname[key])
        # output to file
        F2 = ff.format(PROC,DREL,DSET,gia_mean_str[GIA],FLAGS[0],units,
            LMAX,order_str,gw_str,ds_str,key,START,END,suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(F2)
        val.to_file(OUTPUT_FILE, format=DATAFORM, field_mapping=field_mapping,
            date=False, verbose=VERBOSE, **attributes)
        # change the permissions mode
        OUTPUT_FILE.chmod(mode=MODE)
        # add to output file list
        output_files.append(OUTPUT_FILE)

    # return the list of output files
    return output_files

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates the impact of sea level fingerprints on the
            GRACE/GRACE-FO trend
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
        type=str, default=['rmTWC_rmSLF','rmTWC'],
        nargs=2, help='GRACE/GRACE-FO specific data flags')
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
        # run grace_spatial_differences algorithm with parameters
        output_files = grace_spatial_differences(
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
            FLAGS=args.flag,
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
