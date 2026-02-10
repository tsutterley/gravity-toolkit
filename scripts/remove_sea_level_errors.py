#!/usr/bin/env python
u"""
remove_sea_level_errors.py
Written by Tyler Sutterley (05/2023)
Removes sea level load harmonics and spatial maps after running programs

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for mascon files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -l X, --lmax X: maximum spherical harmonic degree
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    --mascon-type X: input load type (DISC, POINT or CAP)
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --runs X: Number of Monte Carlo iterations
    -V, --verbose: verbose output of processing run

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Written 07/2020
"""
from __future__ import print_function

import sys
import os
import re
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

# PURPOSE: Remove sea level spatial maps to free space
def remove_sea_level_errors(PROC, DREL, DSET,
    DATAFORM=None,
    MASCON_TYPE=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    RUNS=0,
    OUTPUT_DIRECTORY=None):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''

    # output file format for input_distribution and output_slf
    file_format='{0}_MC_ITERATION_{1}{2}{3}_L{4:d}_{5:05d}.{6}'
    # for each monte carlo iteration
    for n in range(RUNS):
        # spherical harmonic and spatial fields from sea level programs
        a1=(MASCON_TYPE,ITERATION,dset_str,ocean_str,EXPANSION,n,suffix)
        a2=('SLF',ITERATION,dset_str,ocean_str,EXPANSION,n,suffix)
        # remove sea level harmonics and spatial files
        FILE1 = OUTPUT_DIRECTORY.joinpath(file_format.format(*a1))
        FILE2 = OUTPUT_DIRECTORY.joinpath(file_format.format(*a2))
        FILE1.unlink()
        FILE2.unlink()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Removes sea level error harmonics and spatial
            maps after running programs
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
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
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # input load type (DISC, POINT or CAP)
    parser.add_argument('--mascon-type','-T',
        type=str.upper, default='CAP', choices=['DISC','POINT','CAP'],
        help='Input load type')
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
    # number of monte carlo iterations
    parser.add_argument('--runs','-R',
        type=int, default=10000,
        help='Number of Monte Carlo iterations')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
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
        # run remove_sea_level_errors algorithm with parameters
        remove_sea_level_errors(args.center, args.release, args.product,
            DATAFORM=args.format,
            MASCON_TYPE=args.mascon_type,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            RUNS=args.runs,
            OUTPUT_DIRECTORY=args.output_directory)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
