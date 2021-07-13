#!/usr/bin/env python
u"""
grace_date.py
Written by Tyler Sutterley (02/2021)

Reads index file from podaac_grace_sync.py or gfz_isdc_grace_ftp.py
Parses dates of each GRACE/GRACE-FO file and assigns the month number
Creates an index of dates for GRACE/GRACE-FO files

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data

OPTIONS:
    PROC: GRACE data processing center (CSR/CNES/JPL/GFZ/GRAZ/COSTG/SWARM)
    DREL: GRACE/GRACE-FO Data Release (RL03 for CNES) (RL06 for CSR/GFZ/JPL)
    DSET: GRACE dataset (GAA/GAB/GAC/GAD/GSM)
        GAA is the non-tidal atmospheric correction
        GAB is the non-tidal oceanic correction
        GAC is the combined non-tidal atmospheric and oceanic correction
        GAD is the GRACE ocean bottom pressure product
        GSM is corrected monthly GRACE/GRACE-FO static field product
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
    Updated 02/2021: use adjust_months function to fix special months cases
    Updated 12/2020: Add SWARM data compilance
    Updated 12/2020: using utilities from time module
    Updated 11/2020: updated for CNES RL04 & RL05 and GRAZ 2018 (monthly fields)
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
import re
import argparse
import numpy as np
import gravity_toolkit.time

def grace_date(base_dir, PROC='', DREL='', DSET='', OUTPUT=True, MODE=0o775):
    """
    Reads index file from podaac_grace_sync.py or gfz_isdc_grace_ftp.py
    Parses dates of each GRACE/GRACE-FO file and assigns the month number
    Creates an index of dates for GRACE/GRACE-FO files

    Arguments
    ---------
    base_dir: working data directory

    Keyword arguments
    -----------------
    PROC: GRACE data processing center
        CSR: University of Texas Center for Space Research
        GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        JPL: Jet Propulsion Laboratory
        CNES: French Centre National D'Etudes Spatiales
        GRAZ: Institute of Geodesy from GRAZ University of Technology
        COSTG: International Combination Service for Time-variable Gravity Fields

        SWARM: gravity data from SWARM satellite
    DREL: GRACE/GRACE-FO data release
    DSET: GRACE/GRACE-FO dataset
        GAA: non-tidal atmospheric correction
        GAB: non-tidal oceanic correction
        GAC: combined non-tidal atmospheric and oceanic correction
        GAD: ocean bottom pressure product
        GSM: corrected monthly static gravity field product
    OUTPUT: create index file of dates for GRACE/GRACE-FO data
    MODE: Permission mode of directories and files

    Returns
    -------
    output_files: dictionary of GRACE/GRACE-FO files indexed by month
    """

    #--  Directory of exact product
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)
    #-- input index file containing GRACE data filenames
    with open(os.path.join(grace_dir, 'index.txt'),'r') as f:
        input_files = f.read().splitlines()

    #--  number of lines in input_files
    n_files = len(input_files)

    #-- define date variables
    start_yr = np.zeros((n_files))#-- year start date
    end_yr = np.zeros((n_files))#-- year end date
    start_day = np.zeros((n_files))#-- day number start date
    end_day = np.zeros((n_files))#-- day number end date
    mid_day = np.zeros((n_files))#-- mid-month day
    tot_days = np.zeros((n_files))#-- number of days since Jan 2002
    tdec = np.zeros((n_files))#-- tdec is the date in decimal form
    mon = np.zeros((n_files,),dtype=np.int)#-- GRACE/GRACE-FO month number

    if PROC in ('CSR', 'GFZ', 'JPL', 'CNES', 'COSTG'):
        #-- compile numerical expression operator for parameters from files
        #-- will work with previous releases and releases for GRACE-FO
        #-- UTCSR: The University of Texas at Austin Center for Space Research
        #-- EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
        #-- GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
        #-- JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
        #-- JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
        #-- GRGS: CNES Groupe de Recherche de Géodésie Spatiale
        regex_pattern = (r'(.*?)-2_(\d+)-(\d+)_(.*?)_({0})_(.*?)_(\d+)(.*?)'
           r'(\.gz|\.gfc|\.txt)?$').format(r'UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC|GRGS|COSTG')
    elif PROC == 'GRAZ':
        # -- GRAZ: Institute of Geodesy from GRAZ University of Technology
        regex_pattern = (r'(.*?)-({0})_(.*?)_(\d+)-(\d+)'
                         r'(\.gz|\.gfc|\.txt)').format(r'Grace_operational|Grace2018')
    elif PROC == 'SWARM':
        # -- SWARM: data from SWARM satellite
        regex_pattern = (r'({0})_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)'
                         r'(\.gz|\.gfc|\.txt)').format(r'SW')
    else:
        raise ValueError("Unknown PROC value:", PROC)

    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- for each data file
    for t, infile in enumerate(input_files):
        #-- extract parameters from input filename
        if PROC in ('CSR', 'GFZ', 'JPL', 'CNES', 'COSTG'):
            PFX,start_date,end_date,AUX,PRC,F1,DRL,F2,SFX = rx.findall(infile).pop()

            #-- find start date, end date and number of days
            start_yr[t] = np.float(start_date[:4])
            end_yr[t] = np.float(end_date[:4])
            start_day[t] = np.float(start_date[4:])
            end_day[t] = np.float(end_date[4:])

        elif PROC == 'GRAZ' or PROC == 'SWARM':
            if PROC == 'GRAZ':
                PFX,SAT,trunc,year,month,SFX = rx.findall(infile).pop()
                #-- find start year, end year
                start_yr[t] = np.float(year)
                end_yr[t] = np.float(year)
            elif PROC == 'SWARM':
                SAT, tmp, PROD, start_date, end_date, RL, SFX = rx.findall(os.path.basename(infile)).pop()

                start_yr[t] = int(start_date[:4])
                end_yr[t] = int(end_date[:4])
                month = int(start_date[4:6])

            #-- Calculation of total days since start of campaign
            #-- Get information on the current year (day per month and day per year)
            dpm = gravity_toolkit.time.dpm_count(start_yr[t])

            #-- find start day, end day
            start_day[t] = np.sum(dpm[:np.int(month) - 1]) + 1
            end_day[t] = np.sum(dpm[:np.int(month)])

        # -- number of days in the starting year for leap and standard years
        dpy = gravity_toolkit.time.calendar_days(start_yr[t]).sum()

        # -- end date taking into account measurements taken on different years
        end_cyclic = (end_yr[t] - start_yr[t]) * dpy + end_day[t]

        #-- Calculation of Mid-month value
        mid_day[t] = np.mean([start_day[t], end_cyclic])

        # -- calculate Modified Julian Day from start_yr and mid_day
        MJD = gravity_toolkit.time.convert_calendar_dates(start_yr[t],
                                                          1.0, mid_day[t], epoch=(1858, 11, 17, 0, 0, 0))
        # -- convert from Modified Julian Days to calendar dates
        cal_date = gravity_toolkit.time.convert_julian(MJD + 2400000.5)

        #-- Calculating the mid-month date in decimal form
        tdec[t] = start_yr[t] + mid_day[t] / dpy

        # -- Calculation of total days since start of campaign
        count = 0
        n_yrs = np.int(start_yr[t] - 2002)
        # -- for each of the GRACE years up to the file year
        for iyr in range(n_yrs):
            # -- year
            year = 2002 + iyr
            # -- add all days from prior years to count
            # -- number of days in year i (if leap year or standard year)
            count += gravity_toolkit.time.calendar_days(year).sum()

        # -- calculating the total number of days since 2002
        tot_days[t] = np.mean([count + start_day[t], count + end_cyclic])

        #-- Calculates the month number (or 10-day number for CNES RL01,RL02)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            mon[t] = np.round(1.0+(tot_days[t]-tot_days[0])/10.0)
        else:
            #-- calculate the GRACE/GRACE-FO month (Apr02 == 004)
            #-- https://grace.jpl.nasa.gov/data/grace-months/
            #-- Notes on special months (e.g. 119, 120) below
            mon[t] = 12*(cal_date['year']-2002) + cal_date['month']

    #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    #-- Accelerometer shutoffs make the relation between month number
    #-- and date more complicated as days from other months are used
    #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    #-- For all: May 2015 (161) is centered in Apr 2015 (160)
    mon = gravity_toolkit.time.adjust_months(mon)

    #-- Output GRACE/GRACE-FO date ascii file
    if OUTPUT:
        date_file = '{0}_{1}_DATES.txt'.format(PROC, DREL)
        fid = open(os.path.join(grace_dir,date_file), 'w')
        #-- date file header information
        args = ('Mid-date','Month','Start_Day','End_Day','Total_Days')
        print('{0} {1:>10} {2:>11} {3:>10} {4:>13}'.format(*args),file=fid)

    #-- create python dictionary mapping input file names with GRACE months
    grace_files = {}
    #-- for each data file
    for t, infile in enumerate(input_files):
        #-- add file to python dictionary mapped to GRACE/GRACE-FO month
        grace_files[mon[t]] = os.path.join(grace_dir,infile)

        #-- print to GRACE DATES ascii file (NOTE: tot_days will be rounded up)
        if OUTPUT:
            print(('{0:13.8f} {1:03d} {2:8.0f} {3:03.0f} {4:8.0f} {5:03.0f} '
                '{6:8.0f}').format(tdec[t],mon[t],start_yr[t],start_day[t],
                end_yr[t],end_day[t],tot_days[t]), file=fid)

    #-- close date file
    #-- set permissions level of output date file
    if OUTPUT:
        fid.close()
        os.chmod(os.path.join(grace_dir, date_file), MODE)

    #-- return the python dictionary that maps GRACE months with GRACE files
    return grace_files

#-- PURPOSE: program that calls grace_date() with set parameters
def main():
    #-- command line parameters
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Parses dates of each GRACE/GRACE-FO file and
            assigns the month number.
            Creates an index of dates for GRACE/GRACE-FO files.
            """
    )
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'],
        choices=['CSR','GFZ','JPL', 'CNES','GRAZ','SWARM', 'COSTG'],
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'],
        help='GRACE/GRACE-FO Data Release')
    #-- GRACE/GRACE-FO data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str.upper, nargs='+',
        default=['GAC','GAD','GSM'],
        choices=['GAA','GAB','GAC','GAD','GSM'],
        help='GRACE/GRACE-FO dealiasing product')
    #-- output GRACE/GRACE-FO ascii date file
    parser.add_argument('--output','-O',
        default=False, action='store_true',
        help='Overwrite existing data')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- run GRACE/GRACE-FO date program
    for pr in args.center:
        for rl in args.release:
            for ds in args.product:
                grace_date(args.directory, PROC=pr, DREL=rl, DSET=ds,
                    OUTPUT=args.output, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
