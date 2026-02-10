#!/usr/bin/env python
u"""
nominal_grace_date.py (05/2023)
Creates a GRACE/GRACE-FO date file using nominal calendar dates

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -r X, --release X: GRACE/GRACE-FO data release
    -Y X, --year X: start and end year to run

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 05/2023: use f-strings for formatting date file
        reduce output columns to fit with other GRACE date programs
        use pathlib to define and operate on paths
    Updated 03/2023: use f-strings for formatting output date lines
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 12/2020: using utilities from time module
    Updated 10/2020: use argparse to set command line parameters
    Updated 04/2020: start at GRACE/GRACE-FO month 004
    Written 03/2018
"""

from __future__ import print_function

import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

def nominal_grace_date(base_dir, DREL, RANGE=None, MODE=0o775):
    # output to AOD1B directory
    PROC, DSET = ('AOD1B', 'GSM')
    # directory of exact product
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    grace_dir = base_dir.joinpath(PROC, DREL, DSET)
    # recursively create directories if they do not exist
    grace_dir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output date file
    grace_date_file = grace_dir.joinpath(f'{PROC}_{DREL}_DATES.txt')
    fid = grace_date_file.open(mode='w', encoding='utf8')
    # date file header information
    a = ('Mid-date','Month','Start_Day','End_Day','Total_Days')
    print('{0} {1:>10} {2:>11} {3:>10} {4:>13}'.format(*a),file=fid)

    # for each year in the range
    count = 0
    for yr in range(RANGE[0], RANGE[1]+1):
        # cumulative days per month (check if year is a leap year)
        if ((np.int64(yr) % 4) == 0):
            # Leap Year
            cdays = [1,32,61,92,122,153,183,214,245,275,306,336,367]
        else:
            # Standard Year
            cdays = [1,32,60,91,121,152,182,213,244,274,305,335,366]
        # for each month
        for m in range(0,12):
            # calculate year decimal and GRACE month
            tdec, = gravtk.time.convert_calendar_decimal(yr,m+1)
            grace_month = gravtk.time.calendar_to_grace(yr,m+1)
            number_of_days = cdays[m+1] - cdays[m]
            print((f'{tdec:13.8f} {grace_month:03d} '
                f'{yr:8.0f} {cdays[m]:03.0f} '
                f'{yr:8.0f} {cdays[m+1]-1:03.0f} '
                f'{count:8.0f}'), file=fid)
            # add to count
            count += number_of_days
    # close the date file
    fid.close()
    # change the permissions mode of the output file
    grace_date_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a GRACE/GRACE-FO date file using nominal
            calendar dates
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'],
        help='GRACE/GRACE-FO Data Release')
    # start and end year to calculate nominal dates
    parser.add_argument('--year','-Y',
        metavar=('start','end'), type=int, nargs=2,
        default=[2002,2023],
        help='Year range to run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
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

    # run program for each release
    for DREL in args.release:
        nominal_grace_date(args.directory, DREL,
            RANGE=args.year, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
