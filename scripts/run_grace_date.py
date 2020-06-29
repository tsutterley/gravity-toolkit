#!/usr/bin/env python
u"""
run_grace_date.py
Written by Tyler Sutterley (05/2020)

Wrapper program for running GRACE date and months programs

Processing Centers: [CSR, GFZ, JPL]
Data Releases CSR, GFZ, JPL: RL06
Data Products: [GAA, GAB, GAC, GAD, GSM]
    CSR: (GAC, GAD, GSM) only

CALLING SEQUENCE:
    run_grace_date(base_dir, PROC, DREL, VERBOSE=False)

INPUTS:
    base_dir: working GRACE/GRACE-FO data directory
    PROC: GRACE processing centers to run
    DREL: GRACE data releases to run

OPTIONS:
    VERBOSE: Track program progress
    MODE: permissions mode of output date files

COMMAND LINE OPTIONS:
    -D X, --directory=X: GRACE/GRACE-FO working data directory
    -C X, --center=X: GRACE/GRACE-FO Processing Center (CSR,GFZ,JPL)
    -R X, --release=X: GRACE/GRACE-FO data releases (RL04,RL05,RL06)
    -V, --verbose: Track program progress
    -M X, --mode=X: Permissions mode of output date files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    grace_date.py: computes the dates of the GRACE datasets
    grace_months_index.py: creates a single file showing the GRACE dates

UPDATE HISTORY:
    Updated 05/2020: updated for public release
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2018: added GFZ RL06 parameters for LMAX 60
        using full release string (RL05 instead of 5)
    Updated 06/2018: using python3 compatible octal and input. RL06 updates.
        if running as program: will use inputs from getopt
    Updated 01/2017: added MODE to set file and directory permissions
    Updated 05/2016: using __future__ print function
    Updated 11/2015: switched to dictionaries only. cleaned up obsolete portions
    Updated 12/2014: added main definition for running from the command line
    Updated 02/2014: minor update to if statements
    Written 07/2012
"""
from __future__ import print_function

import sys
import os
import getopt
from gravity_toolkit.grace_date import grace_date
from gravity_toolkit.grace_months_index import grace_months_index

def run_grace_date(base_dir, PROC, DREL, VERBOSE=False, MODE=0o775):
    #-- allocate python dictionaries for each processing center
    DSET = {}
    VALID = {}
    #-- CSR RL04/5/6 at LMAX 60
    DSET['CSR'] = {'RL04':['GAC', 'GAD', 'GSM'], 'RL05':['GAC', 'GAD', 'GSM'],
        'RL06':['GAC', 'GAD', 'GSM']}
    VALID['CSR'] = ['RL04','RL05','RL06']
    #-- GFZ RL04/5 at LMAX 90
    #-- GFZ RL06 at LMAX 60
    DSET['GFZ'] = {'RL04':['GAA', 'GAB', 'GAC', 'GAD', 'GSM'],
        'RL05':['GAA', 'GAB', 'GAC', 'GAD', 'GSM'],
        'RL06':['GAA', 'GAB', 'GAC', 'GAD', 'GSM']}
    VALID['GFZ'] = ['RL04','RL05','RL06']
    #-- JPL RL04/5/6 at LMAX 60
    DSET['JPL'] = {'RL04':['GAA', 'GAB', 'GAC', 'GAD', 'GSM'],
        'RL05':['GAA', 'GAB', 'GAC', 'GAD', 'GSM'],
        'RL06':['GAA', 'GAB', 'GAC', 'GAD', 'GSM']}
    VALID['JPL'] = ['RL04','RL05','RL06']

    #-- for each processing center
    for p in PROC:
        #-- for each valid dataset release from the processing center
        drel = [r for r in DREL if r in VALID[p]]
        for r in DREL:
            #-- for number of data products
            for d in DSET[p][r]:
                print('GRACE Date Program: {0} {1} {2}'.format(p,r,d))
                #-- run program for processing center, data release and product
                grace_date(base_dir,PROC=p,DREL=r,DSET=d,OUTPUT=True,MODE=MODE)

    #-- run grace months program for data releases
    print('GRACE Months Program')
    grace_months_index(base_dir, DREL=DREL, MODE=MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tGRACE working data directory')
    print(' -C X, --center=X\tGRACE Processing Center (CSR,GFZ,JPL)')
    print(' -R X, --release=X\tGRACE data releases (RL04,RL05,RL06)')
    print(' -V, --verbose\t\t Track program progress')
    print(' -M X, --mode=X\t\tPermissions mode of output files\n')

#-- PURPOSE: program that calls run_grace_date() with set parameters
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','center=','release=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hC:R:VM:',long_options)

    #-- command line options
    #-- working data directory
    base_dir = os.getcwd()
    #-- GRACE Processing Centers to run
    PROC = ['CSR', 'GFZ', 'JPL']
    #-- Data release
    DREL = ['RL06']
    #-- output file information
    VERBOSE = False
    #-- permissions mode of output files (e.g. 0o775)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("-C","--center"):
            PROC = arg.upper().split(',')
        elif opt in ("-R","--release"):
            DREL = arg.split(',')
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- run GRACE preliminary program
    run_grace_date(base_dir, PROC, DREL, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
