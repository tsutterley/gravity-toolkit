#!/usr/bin/env python
u"""
make_grace_index.py
Written by Tyler Sutterley (12/2022)
Creates index files of GRACE/GRACE-FO Level-2 data

CALLING SEQUENCE:
    python make_grace_index.py --version 0 1

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases to run
    -p X, --product X: GRACE/GRACE-FO Data Products to run
    -v X, --version X: GRACE/GRACE-FO Level-2 Data Version to run
    -V, --verbose: Output information for each output file
    -M X, --mode X: Local permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: make the data product optional
    Written 08/2022
"""
from __future__ import print_function

import sys
import os
import logging
import argparse
import gravity_toolkit as gravtk

# PURPOSE: Creates index files of GRACE/GRACE-FO data
def make_grace_index(DIRECTORY, PROC=[], DREL=[], DSET=[],
    VERSION=[], MODE=None):

    # mission shortnames
    shortname = {'grace':'GRAC', 'grace-fo':'GRFO'}
    # GRACE/GRACE-FO level-2 spherical harmonic products
    logging.info('GRACE/GRACE-FO L2 Global Spherical Harmonics:')
    # for each processing center (CSR, GFZ, JPL)
    for pr in PROC:
        # for each data release (RL04, RL05, RL06)
        for rl in DREL:
            # for each level-2 product
            for ds in DSET:
                # local directory for exact data product
                local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                # check if local directory exists
                if not os.access(local_dir, os.F_OK):
                    continue
                # list of GRACE/GRACE-FO files for index
                grace_files = []
                # for each satellite mission (grace, grace-fo)
                for i,mi in enumerate(['grace','grace-fo']):
                    # print string of exact data product
                    logging.info(f'{mi} {pr}/{rl}/{ds}')
                    # regular expression operator for data product
                    rx = gravtk.utilities.compile_regex_pattern(pr, rl, ds,
                        mission=shortname[mi], version=VERSION[i])
                    # find local GRACE/GRACE-FO files to create index
                    granules = [f for f in os.listdir(local_dir) if rx.match(f)]
                    # extend list of GRACE/GRACE-FO files
                    grace_files.extend(granules)

                # outputting GRACE/GRACE-FO filenames to index
                index_file = os.path.join(local_dir,'index.txt')
                with open(index_file, mode='w', encoding='utf8') as fid:
                    for fi in sorted(grace_files):
                        print(fi, file=fid)
                # change permissions of index file
                os.chmod(index_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates index files of GRACE/GRACE-FO
            monthly Level-2 data
            """
    )
    # command line parameters
    # # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # GRACE/GRACE-FO processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'], choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO processing center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'], choices=['RL06'],
        help='GRACE/GRACE-FO data release')
    # GRACE/GRACE-FO data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str.upper, nargs='+',
        default=['GSM'], choices=['GAA','GAB','GAC','GAD','GSM'],
        help='GRACE/GRACE-FO Level-2 data product')
    # GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs=2,
        default=['0','1'], choices=['0','1','2','3'],
        help='GRACE/GRACE-FO Level-2 data version')
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of files created')
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

    # run program with parameters
    make_grace_index(args.directory, PROC=args.center,
        DREL=args.release, DSET=args.product,
        VERSION=args.version, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
