#!/usr/bin/env python
u"""
make_sea_level_shells.py
Written by Tyler Sutterley (05/2023)

Creates a shell script for running sea level variation code

COMMAND LINE OPTIONS:
    -O X, --output-directory X: Output directory for sea level files
    -P X, --file-prefix X: Prefix string for output files
    -D, --date: Input and output file have date information
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -b X, --body X: Treatment of the body tide Love number
        0: Wahr (1981) and Wahr (1985) values from PREM
        1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    -f X, --fluid X: Treatment of the fluid Love number
        0: Han and Wahr (1989) fluid love number
        1: Munk and MacDonald (1960) secular love number
        2: Munk and MacDonald (1960) fluid love number
        3: Lambeck (1980) fluid love number
    --polar-feedback: Include polar feedback
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    -e X, --expansion X: Spherical harmonic expansion
    --mask X: Land-sea mask for redistributing land water flux
    -V, --verbose: Verbose output of run
    -M X, --mode X: permissions mode of the output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: update paths to programs
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: add option to include dates in output files
    Updated 09/2019: check which python version is running
    Updated 11/2018: can vary the land-sea mask
    Written 09/2018
"""
from __future__ import print_function

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

# PURPOSE: create a shell script for running the sea level program
def make_sea_level_shells(input_file,
    DATE=False,
    LOVE_NUMBERS=0,
    BODY_TIDE_LOVE=0,
    FLUID_LOVE=0,
    REFERENCE=None,
    DATAFORM=None,
    EXPANSION=None,
    POLAR_FEEDBACK=False,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    VERBOSE=0,
    MODE=0o775):

    # create output directory if currently non-existent
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # suffix for input ascii, netcdf and HDF5 files
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # sea level equation program
    child_program = 'run_sea_level_equation.py'
    # spherical harmonic expansion flag
    expansion_flag = f' --lmax {EXPANSION:d}'
    # pole tide feedback
    polar_flag = ' --polar-feedback' if POLAR_FEEDBACK else ''
    # Love number flags
    love_flag = f' --love {LOVE_NUMBERS:d}'
    body_flag = f' --body {BODY_TIDE_LOVE:d}'
    fluid_flag = f' --fluid {FLUID_LOVE:d}'
    # Set gravitational load love number of degree 1 to 0.027
    reference_flag = f' --reference {REFERENCE}'
    # input and output data have date information
    date_flag = ' --date' if DATE else ''
    # land-sea mask to use (if not default 0.5x0.5)
    if LANDMASK:
        mask_flag = f' --mask {LANDMASK}'
    else:
        mask_flag = ''
    # number of iterations in sea level program
    # limit the program iterations at 6
    iter_flag = ' -I 6'
    # input and output format flag
    format_flag = f' --format {DATAFORM}'
    # verbosity flag
    verbosity_flag = ' --verbose' if VERBOSE else ''

    # input spherical harmonic data and get GRACE/GRACE-FO months
    Ylms = gravtk.harmonics().from_index(input_file, format=DATAFORM,
        date=True, sort=True)

    # create sea level shell script
    f1 = f'{FILE_PREFIX}INDEX_L{EXPANSION:d}.sh'
    output_shell_script = OUTPUT_DIRECTORY.joinpath(f1)
    fid = output_shell_script.open(mode='w', encoding='utf8')
    # print the path to the shell script
    logging.info(str(output_shell_script))
    # output file format
    file_format = '{0}L{1:d}_{2:03d}.{3}'
    # formatting string for each line in the shell script
    shell_format='{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11} --mode {12} -V {13} {14}'
    # for each grace month and input file
    for grace_month in sorted(Ylms.month):
        # input file
        input_Ylms = Ylms.subset(grace_month)
        input_load, = input_Ylms.filename
        # output file
        args = (FILE_PREFIX,EXPANSION,grace_month,suffix[DATAFORM])
        output_slf = OUTPUT_DIRECTORY.joinpath(file_format.format(*args))
        # print shell script commands
        args = (child_program, expansion_flag, polar_flag, love_flag,
            body_flag, fluid_flag, reference_flag, date_flag, mask_flag,
            iter_flag, format_flag, verbosity_flag, oct(MODE),
            Ylms.compressuser(input_load),
            gravtk.spatial().compressuser(output_slf))
        print(shell_format.format(*args), file=fid)
    # close the shell script
    fid.close()
    # change the permissions mode of the shell script
    output_shell_script.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a shell script for running sea level code
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path,
        help='Input index file with spherical harmonic data files')
    # output working data directory
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for sea level files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for output files')
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Model harmonics are a time series')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # different treatments of the body tide Love numbers of degree 2
    # 0: Wahr (1981) and Wahr (1985) values from PREM
    # 1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    parser.add_argument('--body','-b',
        type=int, default=0, choices=[0,1],
        help='Treatment of the body tide Love number')
    # different treatments of the fluid Love number of gravitational potential
    # 0: Han and Wahr (1989) fluid love number
    # 1: Munk and MacDonald (1960) secular love number
    # 2: Munk and MacDonald (1960) fluid love number
    # 3: Lambeck (1980) fluid love number
    parser.add_argument('--fluid','-f',
        type=int, default=0, choices=[0,1,2,3],
        help='Treatment of the fluid Love number')
    # option for polar feedback
    parser.add_argument('--polar-feedback',
        default=False, action='store_true',
        help='Include effects of polar feedback')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # sea level fingerprint parameters
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # land-sea mask for redistributing mascon mass and land water flux
    parser.add_argument('--mask',
        type=pathlib.Path,
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
        # run make_sea_level_shells algorithm with parameters
        make_sea_level_shells(args.infile,
            DATE=args.date,
            LOVE_NUMBERS=args.love,
            BODY_TIDE_LOVE=args.body,
            FLUID_LOVE=args.fluid,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            EXPANSION=args.expansion,
            POLAR_FEEDBACK=args.polar_feedback,
            LANDMASK=args.mask,
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
