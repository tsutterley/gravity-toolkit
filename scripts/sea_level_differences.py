#!/usr/bin/env python
u"""
sea_level_difference.py
Written by Tyler Sutterley (05/2023)

Calculates the sea level fingerprint error map

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
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    -R X, --runs X: Number of Monte Carlo iterations
    -V, --verbose: verbose output of processing run
    -M X, --mode X: permissions mode of the output files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
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
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 07/2020: rewrite of program to use spatial class for operations
    Updated 10/2019: changing Y/N flags to True/False
        set processing center and release with command line options
    Updated 12/2018: python3 compatibility updates
    Updated 06/2018: using getopt to set parameters
    Written 08/2017
"""
import sys
import os
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

# PURPOSE: calculate the SLF error map
def sea_level_difference(PROC, DREL, DSET,
    START=None,
    END=None,
    DATAFORM=None,
    REDISTRIBUTE_MASCONS=False,
    OUTPUT_DIRECTORY=None,
    ITERATION=None,
    EXPANSION=None,
    LANDMASK=None,
    RUNS=None,
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
    slf_pattern = 'SLF_MC_ITERATION_{0}{1}{2}_L{3:d}_{4:05d}.{5}'
    file_format = 'SLF_MC_ITERATION_{0}{1}{2}_L{3:d}_{4:03d}-{5:03d}.{6}'
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

    # calculate deviations of monte carlo fields from zero
    VARIANCE = landsea.zeros_like()
    VARIANCE.data = np.zeros((nlat, nlon))
    VARIANCE.mask = land_function.astype(bool)
    VARIANCE.fill_value = -9999.0
    for n in range(RUNS):
        # read SLF files
        F1 = slf_pattern.format(ITERATION,dset_str,ocean_str,EXPANSION,
            n,suffix[DATAFORM])
        FILE1 = OUTPUT_DIRECTORY.joinpath(F1)
        val = gravtk.spatial().from_file(FILE1,
            format=DATAFORM, date=False, spacing=[dlon,dlat],
            nlon=nlon, nlat=nlat)
        VARIANCE.data += val.power(2.0).data
        VARIANCE.mask |= val.mask
    # update mask
    VARIANCE.update_mask()
    # calculate RMS of variance
    ERROR = VARIANCE.scale(1.0/(RUNS-1.0)).power(0.5)
    ERROR.update_mask()

    # attributes for output files
    attributes = {}
    attributes['units'] = 'centimeters'
    attributes['longname'] = 'Equivalent_Water_Thickness'
    attributes['title'] = 'Sea_Level_Fingerprint'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # save to file
    F2 = file_format.format(ITERATION,dset_str,ocean_str,EXPANSION,
        START,END,suffix[DATAFORM])
    FILE2 = OUTPUT_DIRECTORY.joinpath(F2)
    ERROR.to_file(FILE2, format=DATAFORM,
        date=False, verbose=VERBOSE, **attributes)
    # change the permissions mode
    FILE2.chmod(mode=MODE)
    # add to output file list
    output_files.append(FILE2)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates the sea level fingerprint error map
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
    # number of monte carlo iterations
    parser.add_argument('--runs','-R',
        type=int, default=10000,
        help='Number of Monte Carlo iterations')
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
        # run sea_level_difference algorithm with parameters
        output_files = sea_level_difference(
            args.center,
            args.release,
            args.product,
            START=args.start,
            END=args.end,
            DATAFORM=args.format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            RUNS=args.runs,
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
