#!/usr/bin/env python
u"""
copy_parameter_files.py
Written by Tyler Sutterley (05/2023)
copies a set of parameter files for the latest months with updated missing

CALLING SEQUENCE:
    python copy_parameter_files.py parameter_file.txt

    Can also input several parameter files in series:
    python copy_parameter_files.py parameter_file1 parameter_file2

    To set the starting month, ending month parameters:
    python copy_parameter_files.py --start 4 --end 154 parameter_file.txt

    To use the latest start and end months, but different MMAX and geocenter:
    python copy_parameter_files.py --mmax 30 --geocenter SLR parameter_file.txt

    To print the input and output filenames
    python copy_parameter_files.py --verbose parameter_file.txt

INPUTS:
    parameter files containing specific variables for each analysis

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working GRACE data directory
    -S X, --start X: Starting GRACE month for output parameter file
    -E X, --end X: Ending GRACE month for output parameter file
    -l X, --lmax X: Output Maximum Spherical Harmonic Degree
    -m X, --mmax X: Output Maximum Spherical Harmonic Order
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    --geocenter X: Geocenter Product to use in parameter file
        None: No degree 1
        Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
            https://grace.jpl.nasa.gov/data/get-data/geocenter/
        SLR: satellite laser ranging from CSR
            ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
        UCI: Sutterley and Velicogna, Remote Sensing (2019)
            https://www.mdpi.com/2072-4292/11/18/2108
    --slr-c20 X: C20 Product to use in parameter file (CSR/GSFC/GFZ)
    --slr-c30 X: C30 Product to use in parameter file (CSR/GSFC/GFZ/LARES)
    -C, --comments: Keep commented lines in output parameter file
    -O, --overwrite: remove previous parameter files with the new iterations
    -V, --verbose: Print input and output file names

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    grace_find_months.py: finds available months for a GRACE/GRACE-FO dataset

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: changed remove index parameter to remove file parameter
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 05/2020: updated for public release
    Updated 07/2019: can replace C30 with coefficients from SLR
    Updated 11/2018: added GRACE data release
    Updated 09/2017: added option overwrite to remove the old parameter files
    Updated 04/2017: changed no input file exception to IOError
    Updated 12/2016: added comments command-line option to keep commented lines
        use previous geocenter product if not entered.
    Updated 06/2016: using DSET parameter with grace_months function
        and can change DSET to a different dataset (--product parameter)
    Updated 05/2016: using __future__ print function, added VERBOSE option
    Updated 03/2016: partial rewrite to use getopt to set output parameters
    Updated 02/2016: added lines for MMAX==LMAX
    Written 12/2015
"""
from __future__ import print_function

import sys
import os
import re
import copy
import logging
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: read listed parameter file and write a new one for month range
def copy_parameter_files(base_dir, parameter_file, NEW_START, NEW_END, NEW_DREL,
    NEW_DSET, LMAX=None, MMAX=None, DEG1=None, SLR_C20=None, SLR_C30=None,
    COMMENTS=False, OVERWRITE=False, MODE=0o775):

    # read prior parameter file
    # Opening parameter file and assigning file ID object (fid)
    parameter_file = pathlib.Path(parameter_file).expanduser().absolute()
    with parameter_file.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()

    # get parameters for GRACE processing center and data release
    # for each line in the file will extract the parameter (name and value)
    parameters = {}
    for fileline in file_contents:
        # Splitting the input line between parameter name and value
        part = fileline.split()
        # filling the parameter definition variable
        parameters[part[0]] = part[1]

    # input and output GRACE processing center
    PROC = parameters['--center']
    # input and output GRACE release
    OLD_DREL = parameters['--release']
    if NEW_DREL is None:
        NEW_DREL = parameters['--release']
        drel_str = ''
    elif (NEW_DREL == OLD_DREL):
        drel_str = ''
    else:
        drel_str = '_{NEW_DREL}'

    # input and output GRACE product
    OLD_DSET = parameters['--product']
    if NEW_DSET is None:
        NEW_DSET = parameters['--product']
        dset_str = ''
    elif (NEW_DSET == OLD_DSET):
        dset_str = ''
    else:
        dset_str = f'_{NEW_DSET}'

    # if changing the truncation
    if LMAX is None:
        LMAX = np.int64(parameters['--lmax'])
    # If MMAX is None == LMAX
    if MMAX is not None:
        order_str = f'M{MMAX:d}'
    else:
        order_str = ''

    # if not entered: keep previous geocenter product in output parameter file
    if DEG1 is None:
        DEG1 = parameters['--geocenter']

    # if not entered: keep previous oblateness product in output parameter file
    if SLR_C20 is None:
        SLR_C20 = parameters['--slr-c20']

    # if not entered: keep previous C30 product in output parameter file
    if SLR_C30 is None:
        SLR_C30 = parameters['--slr-c30']

    # previous date range from parameter file
    OLD_START = np.int64(parameters['--start'])
    OLD_END = np.int64(parameters['--end'])
    OLD_MISSING = np.array(parameters['--missing'].split(','),dtype=np.int64)
    # get the latest GRACE months for output data set
    grace_months = gravtk.grace_find_months(base_dir, PROC, NEW_DREL, DSET=NEW_DSET)
    # default start and end months are the first and latest months
    if NEW_START is None:
        NEW_START = grace_months['start']
    if NEW_END is None:
        NEW_END = grace_months['end']
    # GRACE missing months (keeping any old unique missing months)
    NEW_MISSING = sorted(set(OLD_MISSING) | set(grace_months['missing']))

    # new output file for new months range
    mapping = []
    mapping.append(('parameters', f'parameters{drel_str}{dset_str}'))
    mapping.append((f'L{LMAX:d}', f'L{LMAX:d}{order_str}'))
    mapping.append((f'{OLD_START:03d}-{OLD_END:03d}', f'{NEW_START:03d}-{NEW_END:03d}'))
    FILE = copy.copy(parameter_file.name)
    for key, val in mapping:
        FILE = FILE.replace(key, val)
    # print old and new files
    output_file = parameter_file.with_name(FILE)
    logging.info(f'{str(parameter_file)} ->\n\t{str(output_file)}\n')
    # open file ID object for FILE
    fid = output_file.open(mode='w', encoding='utf8')
    # for each line in the original parameter file
    for fileline in file_contents:
        # Splitting the input line between parameter name and value
        part = fileline.split()
        comments = ' '.join(p for p in part[2:])
        # filling the parameter definition variables
        regex = r'(\-\-output\-directory|\-\-reconstruct\-file|\-\-remove\-file)'
        if bool(re.match(regex,part[0])):
            # add MMAX string for order 30 solutions
            # replace the old starting month with the new one
            # replace the old ending month with the new one
            # replace the old dataset string with the new dataset string
            mapping = []
            mapping.append(('L60', f'L60{order_str}'))
            mapping.append((f'{OLD_START:03d}-{OLD_END:03d}'
                f'{NEW_START:03d}-{NEW_END:03d}'))
            mapping.append((OLD_DREL, NEW_DREL))
            mapping.append((OLD_DSET, NEW_DSET))
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif bool(re.match(rf'#[\s]?{regex}', part[0])) and COMMENTS:
            # add MMAX string for order 30 solutions
            # replace the old starting month with the new one
            # replace the old ending month with the new one
            # replace the old dataset string with the new dataset string
            mapping = []
            mapping.append(('L60', f'L60{order_str}'))
            mapping.append((f'{OLD_START:03d}-{OLD_END:03d}',
                f'{NEW_START:03d}-{NEW_END:03d}'))
            mapping.append((OLD_DSET, NEW_DSET))
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif (part[0] == '--start'):
            # replace the old starting month with the new one
            mapping = [(f'{OLD_START:d}', f'{OLD_START:d}')]
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif (part[0] == '--end'):
            # replace the old ending month with the new one
            mapping = [(f'{OLD_END:d}', f'{NEW_END:d}')]
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif (part[0] == '--missing'):
            # update the missing months
            missing = ' '.join(f'{m:d}' for m in NEW_MISSING)
            print(f'{part[0]}\t{missing}\t{comments}', file=fid)

        elif (part[0] == '--release'):
            # replace the old release with the new one
            mapping = [(OLD_DREL, NEW_DREL)]
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif (part[0] == '--product'):
            # replace the old product with the new one
            mapping = [(OLD_DSET, NEW_DSET)]
            for key, val in mapping:
                fileline = fileline.replace(key, val)
            print(fileline, file=fid)

        elif (part[0] == '--geocenter'):
            # Geocenter correction to DEG1
            print(f'{part[0]}\t{DEG1}\t{comments}', file=fid)

        elif (part[0] == '--slr-c20'):
            # Oblateness correction
            print(f'{part[0]}\t{SLR_C20}\t{comments}', file=fid)

        elif (part[0] == '--slr-c30'):
            # C30 correction
            print(f'{part[0]}\t{SLR_C30}\t{comments}', file=fid)

        elif (part[0] == '--lmax'):
            # LMAX to new spherical harmonic degree
            print(f'{part[0]}\t{LMAX:d}\t{comments}', file=fid)

        elif (part[0] == '--mmax') and (MMAX is not None):
            # MMAX to new spherical harmonic order
            print(f'{part[0]}\t{MMAX:d}\t{comments}', file=fid)

        elif not bool(re.match(r'#', part[0])):
            # keep only uncommented lines
            print(fileline, file=fid)

        elif bool(re.match(r'#', part[0])) and COMMENTS:
            # keep commented lines if COMMENTS is set
            print(fileline, file=fid)

    # close the output parameter file
    fid.close()
    # change the permissions mode of the output parameter file
    output_file.chmod(mode=MODE)
    # remove the previous parameter file if OVERWRITE
    if OVERWRITE:
        parameter_file.unlink()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Copies a set of parameter files for the latest
            months with updated missing
            """
    )
    # command line parameters
    parser.add_argument('parameters',
        type=pathlib.Path, nargs='+',
        help='Parameter files containing specific variables for each analysis')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=None,
        help='Starting GRACE/GRACE-FO month for output file')
    parser.add_argument('--end','-E',
        type=int, default=None,
        help='Ending GRACE/GRACE-FO month for output file')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default=None,
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default=None,
        help='GRACE/GRACE-FO Data Product')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # geocenter, oblateness and low degree zonals
    parser.add_argument('--geocenter',
        type=str, metavar='DEG1', default=None,
        help='Geocenter Product to use in file')
    parser.add_argument('--slr-c20',
        type=str, metavar='C20', default=None,
        help='C20 Product to use in file')
    parser.add_argument('--slr-c30',
        type=str, metavar='C30', default=None,
        help='C30 Product to use in file')
    # keep commented lines
    parser.add_argument('--comments','-C',
        default=False, action='store_true',
        help='Keep commented lines in parameter file')
    # keep commented lines
    parser.add_argument('--overwrite','-O',
        default=False, action='store_true',
        help='Remove previous version of parameter file')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the output files (octal)
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

    # create logger for verbosity level
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run directly for file
    for f in args.parameters:
        copy_parameter_files(args.directory, f, args.start, args.end,
            args.release, args.product,
            LMAX=args.lmax,
            MMAX=args.mmax,
            DEG1=args.geocenter,
            SLR_C20=args.slr_c20,
            SLR_C30=args.slr_c30,
            COMMENTS=args.comments,
            OVERWRITE=args.overwrite,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
