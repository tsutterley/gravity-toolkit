#!/usr/bin/env python
u"""
sea_level_mean.py
Written by Tyler Sutterley (05/2023)

Calculates the mean sea level fingerprint map for each GIA model

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
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
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    -V, --verbose: verbose output of processing run
    -M X, --mode X: permissions mode of the output files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
        read from merged data and error file as output from regression program
    Updated 03/2023: updated inputs to spatial from_file function
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        simplified file imports and exports using wrappers in spatial utility
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 07/2020: rewrite of program to use spatial class for operations
    Updated 10/2019: changing Y/N flags to True/False
    Updated 02/2019: added option for setting GRACE data release
    Updated 12/2018: python3 compatibility updates
    Updated 09/2017: simplified outputs variables for netCDF4 files
    Updated 08/2017: set grid dimensions, added ICE-6G VM5 GIA models
        calculate mean for different sea level iterations
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
# A et al. (2013)
gia_mean_str['AW13-ICE6G'] = '_AW13_ICE6G_mean'

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate the mean SLF map for each GIA model
def sea_level_mean(PROC, DREL, DSET,
    START=None,
    END=None,
    GIA=None,
    GIA_FILES=None,
    DATAFORM=None,
    REDISTRIBUTE_MASCONS=False,
    OUTPUT_DIRECTORY=None,
    ITERATION=None,
    EXPANSION=None,
    LANDMASK=None,
    VERBOSE=0,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''

    # input and output file formats
    file_format = 'SLF_ITERATION_{0}{1}{2}{3}_L{4:d}_{5}_{6:03d}-{7:03d}.{8}'
    # list of output files
    output_files = []

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF4 file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # create land function
    land_function = np.zeros((nlat, nlon),dtype=np.float64)
    # combine land and island levels for land function
    indy,indx = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indy,indx] = 1.0

    # allocate lists for input variables
    dinput = {}
    dinput['x1'] = []
    dinput['x2'] = []
    dinput['AIC_x0'] = []
    dinput['AIC_x1'] = []
    dinput['AIC_x2'] = []
    # number of rheologies and ice histories to run
    N = len(GIA_FILES)
    # iterate GIA models
    for h,GIA_FILE in enumerate(GIA_FILES):
        # input GIA spherical harmonic datafiles
        GIA_Ylms_rate = gravtk.gia().from_GIA(GIA_FILE, GIA=GIA)
        gia_str = f'_{GIA_Ylms_rate.title}'
        # read mass trend files
        for key in ['x1','x2']:
            F1 = file_format.format(ITERATION, dset_str, gia_str,
                ocean_str, EXPANSION, key, START, END, suffix[DATAFORM])
            INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
            field_mapping = dict(lon='lon', lat='lat', data='data', error='error')
            temp = gravtk.spatial().from_file(INPUT_FILE,
                format=DATAFORM, date=False, field_mapping=field_mapping,
                spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
            dinput[key].append(temp)
        # read AIC files
        for key in ['AIC_x0','AIC_x1','AIC_x2']:
            F1 = file_format.format(ITERATION, dset_str, gia_str,
                ocean_str, EXPANSION, key, START, END, suffix[DATAFORM])
            INPUT_FILE = OUTPUT_DIRECTORY.joinpath(F1)
            field_mapping = dict(lon='lon', lat='lat', data='data')
            temp = gravtk.spatial().from_file(INPUT_FILE,
                format=DATAFORM, date=False, field_mapping=field_mapping,
                spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
            dinput[key].append(temp)

    # create combined spatial objects
    output = {}
    units = {}
    # calculate mean GIA-corrected x1 change (for all Earth rheologies)
    x1 = gravtk.spatial().from_list(dinput['x1'],date=False)
    output['x1'] = x1.mean()
    units['x1'] = '{0} yr^{1:d}'.format('centimeters',-1)
    # calculate mean acceleration x2 change (for all Earth rheologies)
    x2 = gravtk.spatial().from_list(dinput['x2'],date=False)
    output['x2'] = x2.mean().scale(2.0)
    units['x2'] = '{0} yr^{1:d}'.format('centimeters',-2)
    # GRACE satellite error component
    e1 = x1.copy(); e1.data = np.copy(x1.error)
    output['x1'].error = e1.sum(power=2.0).scale(1.0/N).power(0.5).data
    output['x1'].update_mask()
    e2 = x2.copy(); e2.data = np.copy(x2.error)
    output['x2'].error = e2.sum(power=2.0).scale(4.0/N).power(0.5).data
    output['x2'].update_mask()
    # significance means
    AICx0 = gravtk.spatial().from_list(dinput['AIC_x0'],date=False).mean()
    AICx1 = gravtk.spatial().from_list(dinput['AIC_x1'],date=False).mean()
    AICx2 = gravtk.spatial().from_list(dinput['AIC_x2'],date=False).mean()

    # masked trend values
    ii,jj = np.nonzero((np.abs(output['x1'].data) <= output['x1'].error) |
        (AICx1.data >= AICx0.data))
    output['MASKED_x1'] = output['x1'].copy()
    output['MASKED_x1'].mask[ii,jj] = True
    output['MASKED_x1'].update_mask()
    units['MASKED_x1'] = '{0} yr^{1:d}'.format('centimeters',-1)
    # write masked acceleration values to file
    ii,jj = np.nonzero((np.abs(output['x2'].data) <= output['x2'].error) |
        (AICx2.data >= AICx1.data))
    output['MASKED_x2'] = output['x2'].copy()
    output['MASKED_x2'].mask[ii,jj] = True
    output['MASKED_x2'].update_mask()
    units['MASKED_x2'] = '{0} yr^{1:d}'.format('centimeters',-2)

    # attributes for output files
    longname = 'Equivalent_Water_Thickness'
    title = 'Sea_Level_Fingerprint'
    # output data to file
    for key,val in output.items():
        F2 = file_format.format(ITERATION, dset_str, gia_mean_str[GIA],
            ocean_str, EXPANSION, key, START, END, suffix[DATAFORM])
        OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(F2)
        output_data(val, FILENAME=OUTPUT_FILE, DATAFORM=DATAFORM,
            UNITS=units[key], LONGNAME=longname, TITLE=title,
            CONF=0.95, VERBOSE=VERBOSE, MODE=MODE)
        # add to output file list
        output_files.append(OUTPUT_FILE)

    # return the list of output files
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
    attributes['data']['description'] = 'Model_fit'
    attributes['data']['long_name'] = LONGNAME
    attributes['data']['units'] = UNITS
    attributes['error'] = {}
    attributes['error']['description'] = 'Uncertainty_in_model_fit'
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
        description="""Calculates the mean sea level fingerprint map for
            each glacial isostatic adjustment model
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for spatial files')
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
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
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
    # mascon parameters
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # sea level fingerprint parameters
    parser.add_argument('--iteration','-I',
        type=int, default=1,
        help='Sea level fingerprint iteration')
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
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
        # run sea_level_mean algorithm with parameters
        output_files = sea_level_mean(
            args.center,
            args.release,
            args.product,
            START=args.start,
            END=args.end,
            GIA=args.gia,
            GIA_FILES=args.gia_file,
            DATAFORM=args.format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
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
