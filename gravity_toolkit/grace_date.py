#!/usr/bin/env python
u"""
grace_date.py
Written by Tyler Sutterley (11/2022)
Contributions by Hugo Lecomte and Yara Mohajerani

Reads index file from podaac_grace_sync.py or gfz_isdc_grace_ftp.py
Parses dates of each GRACE/GRACE-FO file and assigns the month number
Creates an index of dates for GRACE/GRACE-FO files

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data

OPTIONS:
    PROC: Data processing center or satellite mission
        CSR: University of Texas Center for Space Research
        GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        JPL: Jet Propulsion Laboratory
        CNES: French Centre National D'Etudes Spatiales
        GRAZ: Institute of Geodesy from GRAZ University of Technology
        COSTG: Combination Service for Time-variable Gravity Fields
        Swarm: Time-variable gravity data from Swarm satellites
    DREL: GRACE/GRACE-FO/Swarm data release
    DSET: GRACE/GRACE-FO/Swarm dataset
        GAA: non-tidal atmospheric correction
        GAB: non-tidal oceanic correction
        GAC: combined non-tidal atmospheric and oceanic correction
        GAD: ocean bottom pressure product
        GSM: monthly static field product
    OUTPUT: create index of dates for GRACE/GRACE-FO data
    MODE: permissions mode of output files

OUTPUTS:
    dictionary of GRACE/GRACE-FO files indexed by month

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: raise exception if index file cannot be found
        use logging for debugging level verbose output
    Updated 08/2022: moved file parsing functions to time module
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 09/2021: adjust regular expression operators for Swarm and GRAZ
        use functions for converting to and from GRACE months
    Updated 07/2021: remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: using utilities from time module
        Add Swarm data time-variable gravity fields
    Updated 11/2020: added CNES RL04 & RL05 and GRAZ 2018 (monthly fields)
    Updated 10/2020: use argparse to set command line parameters
    Updated 07/2020: added function docstrings
    Updated 03/2020: for public release
    Updated 11/2018: updated regular expression pattern for RL06 GFZ
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: can read new GRACE file format for RL06 and GRACE-FO
    Updated 10/2017: added option OUTPUT to write the text file with GRACE dates
        now will output from the function the GRACE month and file name
        adjusted regular expression for extracting parameters from filename
    Updated 02/2017: added mode to modify the permissions of the output file
    Updated 05-06/2016: using __future__ print function, format date lines
    Updated 01/2016: minor clean up.  using enumerate for loop
    Updated 10/2015: added manual fix for month 161 (centered in 160)
    Updated 05/2015: additional 2 digits to date file to reduce the differences
        in dates if importing from date file or binary data file
    Updated 03/2015: added main definition for running from command line
        further generalization, calculates the Julian date of the GRACE date
    Updated 11/2014: output more decimal points for date
    Updated 09/2014: using regular expressions, general code updates
    Updated 03/2014: printing GRACE month with zero-padding
    Updated 02/2014: minor update to if statements
    Updated 01/2014: updated for CNES RL03 (monthly fields)
    Updated 10/2013: created a directory with both RL04 and RL05 (combining)
        made a slight edit for the drift rates
    Updated 09/2013: changed dating schemes for all products
        for CNES: Solution 1 == 001, versus 10-day from start of 2002
    Updated 05/2013: modified for use with GUI program
    Updated 07/2012: changed some variable names, and saving variables
        start_yr, start_day, end_yr, end_day
    Updated 07/2012: fix missing date for CSR RL05
        One of the months 119 registered as 118 as half of 119 is in 118 due to
        accelerometer issues during 119
    Updated 06/2012: fixes missing dates to be automatic
        Also added options to enter the processing center, the data release and
        the dataset from an external main level program
    Updated 04/2012: changes for RL05 data
"""
from __future__ import print_function

import os
import logging
import argparse
import numpy as np
import gravity_toolkit.time

# PURPOSE: parses GRACE/GRACE-FO data files and assigns month numbers
def grace_date(base_dir, PROC='', DREL='', DSET='', OUTPUT=True, MODE=0o775):
    """
    Reads index file containing GRACE/GRACE-FO/Swarm data files

    Parses dates of each GRACE/GRACE-FO file and assigns the month number

    Creates an index of dates for GRACE/GRACE-FO files

    Parameters
    ----------
    base_dir: str
        Working data directory
    PROC: str, default ''
        Data processing center or satellite mission

            - ``'CSR'``: University of Texas Center for Space Research
            - ``'GFZ'``: German Research Centre for Geosciences (GeoForschungsZentrum)
            - ``'JPL'``: Jet Propulsion Laboratory
            - ``'CNES'``: French Centre National D'Etudes Spatiales
            - ``'GRAZ'``: Institute of Geodesy from GRAZ University of Technology
            - ``'COSTG'``: Combination Service for Time-variable Gravity Fields
            - ``'Swarm'``: Time-variable gravity data from Swarm satellites
    DREL: str, default ''
        GRACE/GRACE-FO/Swarm data release
    DSET: str, default ''
        GRACE/GRACE-FO/Swarm dataset

            - ``'GAA'``: non-tidal atmospheric correction
            - ``'GAB'``: non-tidal oceanic correction
            - ``'GAC'``: combined non-tidal atmospheric and oceanic correction
            - ``'GAD'``: ocean bottom pressure product
            - ``'GSM'``: corrected monthly static gravity field product
    OUTPUT: bool, default True
        create index file of dates for GRACE/GRACE-FO data
    MODE: oct, default 0o775
        Permission mode of directories and files

    Returns
    -------
    output_files: dict
        dictionary of GRACE/GRACE-FO files indexed by month
    """

    #  Directory of exact product
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)
    # index file containing GRACE/GRACE-FO data filenames
    index_file = os.path.join(grace_dir, 'index.txt')
    # check that index file exists
    if not os.access(index_file, os.F_OK):
        raise FileNotFoundError(f'{index_file} not found')
    # log index file if debugging
    logging.debug(f'Reading index file: {index_file}')
    # read index file for GRACE/GRACE-FO filenames
    with open(index_file, mode='r', encoding='utf8') as f:
        input_files = f.read().splitlines()

    #  number of lines in input_files
    n_files = len(input_files)

    # define date variables
    start_yr = np.zeros((n_files))# year start date
    end_yr = np.zeros((n_files))# year end date
    start_day = np.zeros((n_files))# day number start date
    end_day = np.zeros((n_files))# day number end date
    mid_day = np.zeros((n_files))# mid-month day
    tot_days = np.zeros((n_files))# number of days since Jan 2002
    tdec = np.zeros((n_files))# tdec is the date in decimal form
    mon = np.zeros((n_files,),dtype=np.int64)# GRACE/GRACE-FO month number

    # for each data file
    for t,infile in enumerate(input_files):
        if PROC in ('GRAZ','Swarm',):
            # get date lists for the start and end of fields
            start_date,end_date = gravity_toolkit.time.parse_gfc_file(
                infile, PROC, DSET)
            # start and end year
            start_yr[t] = np.float64(start_date[0])
            end_yr[t] = np.float64(end_date[0])
            # number of days in each month for the calendar year
            dpm = gravity_toolkit.time.calendar_days(start_yr[t])
            # start and end day of the year
            start_day[t] = np.sum(dpm[:start_date[1]-1]) + start_date[2] + \
                start_date[3]/24. + start_date[4]/1440. + start_date[5]/86400.
            end_day[t] = np.sum(dpm[:end_date[1]-1]) + end_date[2] + \
                end_date[3]/24. + end_date[4]/1440. + end_date[5]/86400.
        else:
            # get date lists for the start and end of fields
            start_date,end_date = gravity_toolkit.time.parse_grace_file(infile)
            # start and end year
            start_yr[t] = np.float64(start_date[0])
            end_yr[t] = np.float64(end_date[0])
            # start and end day of the year
            start_day[t] = np.float64(start_date[1])
            end_day[t] = np.float64(end_date[1])

        # number of days in the starting year for leap and standard years
        dpy = gravity_toolkit.time.calendar_days(start_yr[t]).sum()
        # end date taking into account measurements taken on different years
        end_cyclic = (end_yr[t]-start_yr[t])*dpy + end_day[t]
        # calculate mid-month value
        mid_day[t] = np.mean([start_day[t], end_cyclic])

        # calculate Modified Julian Day from start_yr and mid_day
        MJD = gravity_toolkit.time.convert_calendar_dates(start_yr[t],
            1.0,mid_day[t],epoch=(1858,11,17,0,0,0))
        # convert from Modified Julian Days to calendar dates
        cal_date = gravity_toolkit.time.convert_julian(MJD+2400000.5)

        # Calculating the mid-month date in decimal form
        tdec[t] = start_yr[t] + mid_day[t]/dpy

        # Calculation of total days since start of campaign
        count = 0
        n_yrs = np.int64(start_yr[t]-2002)
        # for each of the GRACE years up to the file year
        for iyr in range(n_yrs):
            # year
            year = 2002 + iyr
            # add all days from prior years to count
            # number of days in year i (if leap year or standard year)
            count += gravity_toolkit.time.calendar_days(year).sum()

        # calculating the total number of days since 2002
        tot_days[t] = np.mean([count+start_day[t], count+end_cyclic])

        # Calculates the month number (or 10-day number for CNES RL01,RL02)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            mon[t] = np.round(1.0+(tot_days[t]-tot_days[0])/10.0)
        else:
            # calculate the GRACE/GRACE-FO month (Apr02 == 004)
            # https://grace.jpl.nasa.gov/data/grace-months/
            # Notes on special months (e.g. 119, 120) below
            mon[t] = gravity_toolkit.time.calendar_to_grace(
                cal_date['year'],cal_date['month'])

    # The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    # Accelerometer shutoffs make the relation between month number
    # and date more complicated as days from other months are used
    # For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    # For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    # For all: May 2015 (161) is centered in Apr 2015 (160)
    mon = gravity_toolkit.time.adjust_months(mon)

    # Output GRACE/GRACE-FO date ascii file
    if OUTPUT:
        date_file = f'{PROC}_{DREL}_DATES.txt'
        fid = open(os.path.join(grace_dir,date_file), mode='w', encoding='utf8')
        # date file header information
        args = ('Mid-date','Month','Start_Day','End_Day','Total_Days')
        print('{0} {1:>10} {2:>11} {3:>10} {4:>13}'.format(*args),file=fid)

    # create python dictionary mapping input file names with GRACE months
    grace_files = {}
    # for each data file
    for t, infile in enumerate(input_files):
        # add file to python dictionary mapped to GRACE/GRACE-FO month
        grace_files[mon[t]] = os.path.join(grace_dir,infile)
        # print to GRACE dates ascii file (NOTE: tot_days will be rounded)
        if OUTPUT:
            print(('{0:13.8f} {1:03d} {2:8.0f} {3:03.0f} {4:8.0f} {5:03.0f} '
                '{6:8.0f}').format(tdec[t],mon[t],start_yr[t],start_day[t],
                end_yr[t],end_day[t],tot_days[t]), file=fid)

    # close date file
    # set permissions level of output date file
    if OUTPUT:
        fid.close()
        os.chmod(os.path.join(grace_dir, date_file), MODE)

    # return the python dictionary that maps GRACE months with GRACE files
    return grace_files

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Parses dates of each GRACE/GRACE-FO file and
            assigns the month number.
            Creates an index of dates for GRACE/GRACE-FO files.
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # Data processing center or satellite mission
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'],
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str.upper, nargs='+',
        default=['GAC','GAD','GSM'],
        choices=['GAA','GAB','GAC','GAD','GSM'],
        help='GRACE/GRACE-FO Level-2 data product')
    # output GRACE/GRACE-FO ascii date file
    parser.add_argument('--output','-O',
        default=False, action='store_true',
        help='Overwrite existing data')
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

    # run GRACE/GRACE-FO date program
    for pr in args.center:
        for rl in args.release:
            for ds in args.product:
                grace_date(args.directory, PROC=pr, DREL=rl, DSET=ds,
                    OUTPUT=args.output, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
