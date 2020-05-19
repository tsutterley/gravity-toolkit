#!/usr/bin/env python
u"""
grace_date.py
Written by Tyler Sutterley (03/2020)

Reads dates of each GRACE/GRACE-FO data file and assigns the month number
    reads the start and end date from the filename,
    calculates the mean date in decimal format (correcting for leap years)
Creates an index of dates for GRACE/GRACE-FO data if specified

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data

OPTIONS:
    PROC: GRACE data processing center (CSR/CNES/JPL/GFZ)
    DREL: GRACE data release (RL03 for CNES) (RL06 for CSR/GFZ/JPL)
    DSET: GRACE dataset (GAA/GAB/GAC/GAD/GSM)
        GAA is the non-tidal atmospheric correction
        GAB is the non-tidal oceanic correction
        GAC is the combined non-tidal atmospheric and oceanic correction
        GAD is the GRACE ocean bottom pressure product
        GSM is corrected monthly GRACE/GRACE-FO static field product
    OUTPUT: create index of dates for GRACE/GRACE-FO data

OUTPUTS:
    dictionary of files mapped by GRACE/GRACE-FO month

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    future: Compatibility layer between Python 2 and Python 3
        (https://python-future.org/)

PROGRAM DEPENDENCIES:
    convert_julian.py: converts a julian date into a calendar date

UPDATE HISTORY:
    Updated 03/2020 for public release
"""
from __future__ import print_function

import sys
import os
import re
import getopt
import numpy as np
from gravity_toolkit.convert_julian import convert_julian

def grace_date(base_dir, PROC='', DREL='', DSET='', OUTPUT=True, MODE=0o775):
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
    JD = np.zeros((n_files))#-- Julian date of mid-month
    tot_days = np.zeros((n_files))#-- number of days since Jan 2002
    tdec = np.zeros((n_files))#-- tdec is the date in decimal form
    mon = np.zeros((n_files,),dtype=np.int)#-- GRACE/GRACE-FO month number

    #-- compile numerical expression operator for parameters from files
    #-- will work with previous releases and releases for GRACE-FO
    #-- UTCSR: The University of Texas at Austin Center for Space Research
    #-- EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
    #-- GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
    #-- JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
    #-- JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    regex_pattern = ('(.*?)-2_(\d+)-(\d+)_(.*?)_({0})_(.*?)_(\d+)(.*?)'
        '(\.gz|\.gfc)?$').format('UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC')
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- Output GRACE date ascii file
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
        #-- extract parameters from input filename
        PFX,start_date,end_date,AUX,PRC,F1,DRL,F2,SFX = rx.findall(infile).pop()
        #-- find start date, end date and number of days
        start_yr[t] = np.float(start_date[:4])
        end_yr[t] = np.float(end_date[:4])
        start_day[t] = np.float(start_date[4:])
        end_day[t] = np.float(end_date[4:])
        #-- end_day (will be changed if the month crosses 2 years)
        end_plus = np.copy(end_day[t])

        #-- calculate mid-month date taking into account if measurements are
        #-- on different years
        if ((start_yr[t] % 4) == 0):#-- Leap Year (% = modulus)
            dpy = 366.0
        else:#-- Standard Year
            dpy = 365.0
        #-- For data that crosses years
        if (start_yr[t] != end_yr[t]):
            #-- end_yr - start_yr should be 1
            end_plus = (end_yr[t]-start_yr[t])*dpy + end_day[t]
        #-- Calculation of Mid-month value
        mid_day[t] = np.mean([start_day[t], end_plus])

        #-- Calculation of the Julian date from start_yr and mid_day
        JD[t] = np.float(367.0*start_yr[t] - \
            np.floor(7.0*(start_yr[t] + np.floor(10.0/12.0))/4.0) - \
            np.floor(3.0*(np.floor((start_yr[t] - 8.0/7.0)/100.0) + 1.0)/4.0) +\
            np.floor(275.0/9.0) + mid_day[t] + 1721028.5)
        #-- convert the julian date into calendar dates (hour, day, month, year)
        cal_date = convert_julian(JD[t])

        #-- Calculating the mid-month date in decimal form
        tdec[t] = start_yr[t] + mid_day[t]/dpy

        #-- Calculation of total days since start of campaign
        count = 0
        n_yrs = np.int(start_yr[t]-2002)
        #-- for each of the GRACE years up to the file year
        for iyr in range(n_yrs):
            #-- year i
            year = 2002 + iyr
            #-- number of days in year i (if leap year or standard year)
            if ((year % 4) == 0):
                #-- Leap Year
                dpm=[31,29,31,30,31,30,31,31,30,31,30,31]
            else:
                #-- Standard Year
                dpm=[31,28,31,30,31,30,31,31,30,31,30,31]
            #-- add all days from prior years to count
            count += np.sum(dpm)

        #-- calculating the total number of days since 2002
        tot_days[t] = np.mean([count+start_day[t], count+end_plus])

        #-- Calculates the month number (or 10-day number for CNES RL01,RL02)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            mon[t] = np.round(1.0+(tot_days[t]-tot_days[0])/10.0)
        else:
            #-- calculate the GRACE/GRACE-FO month (Apr02 == 004)
            #-- https://grace.jpl.nasa.gov/data/grace-months/
            #-- Notes on special months (e.g. 119, 120) below
            mon[t] = 12*(cal_date['year']-2002) + cal_date['month']

            #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
            #-- Accelerometer shutoffs make this relation between month number
            #-- and date more complicated as days from other months are used
            #-- For CSR and GFZ: Nov11 (month 119) is centered in Oct11 (118)
            #-- For JPL: Dec 2011 (month 120) is centered in Jan12 (121)
            #-- For all: May15 (month 161) is centered in Apr15 (160)
            if PROC in ('CSR','GFZ') and (mon[t] == mon[t-1]) and (mon[t-1] == 118):
                mon[t] = mon[t-1] + 1
            elif (mon[t] == mon[t-1]) and (mon[t-1] == 160):
                mon[t] = mon[t-1] + 1
            elif PROC in ('JPL') and (mon[t-1] == 119):
                mon[t] = mon[t-1] + 1

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

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' --directory=X\t\tGRACE/GRACE-FO working directory')
    print(' -C X, --center=X\tGRACE/GRACE-FO Processing Center (CSR,GFZ,JPL)')
    print(' -R X, --release=X\tGRACE/GRACE-FO data releases (RL04,RL05,RL06)')
    print(' -D X, --dataset=X\tGRACE/GRACE-FO dataset (GAC,GAD,GSM)')
    print(' -O, --output\t\tOutput GRACE/GRACE-FO ascii date file')
    print(' -M X, --mode=X\t\tPermissions mode of output files\n')

#-- PURPOSE: program that calls grace_date() with set parameters
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','center=','release=','dataset=',
        'output','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hC:R:D:OM:',long_options)

    #-- GRACE/GRACE-FO directory
    base_dir = os.getcwd()
    #-- GRACE/GRACE-FO Processing Centers to run
    PROC = ['CSR','GFZ','JPL']
    #-- Data release
    DREL = ['RL06']
    #-- Dataset
    DSET = ['GAC','GAD','GSM']
    #-- output GRACE/GRACE-FO ascii date file
    OUTPUT = False
    #-- permissions mode of output files (e.g. 0o775)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("-C","--center"):
            PROC = arg.upper().split(',')
        elif opt in ("-R","--release"):
            DREL = arg.upper().split(',')
        elif opt in ("-D","--dataset"):
            DSET = arg.upper().split(',')
        elif opt in ("-O","--output"):
            OUTPUT = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- run GRACE/GRACE-FO date program
    for pr in PROC:
        for rl in DREL:
            for ds in DSET:
                grace_date(base_dir, PROC=pr, DREL=rl, DSET=ds,
                    OUTPUT=OUTPUT, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
