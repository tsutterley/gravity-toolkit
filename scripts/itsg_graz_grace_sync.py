#!/usr/bin/env python
u"""
itsg_graz_grace_sync.py
Written by Tyler Sutterley (12/2022)
Syncs GRACE/GRACE-FO and auxiliary data from the ITSG GRAZ server

CALLING SEQUENCE:
    python itsg_graz_grace_sync.py --release Grace2018 --lmax 60

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -r X, --release X: GRAZ Data Releases to sync
        Grace2014
        Grace2016
        Grace2018
        Grace_operational
    --lmax X: Maximum degree and order of GRAZ products
        60
        96
        120
    -t X, --timeout X: Timeout in seconds for blocking operations
    --log: output log of files downloaded
    --list: print files to be transferred, but do not execute transfer
    --clobber: Overwrite existing data in transfer
    -M X, --mode X: permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Written 09/2021
"""
from __future__ import print_function

import sys
import os
import re
import time
import shutil
import logging
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: sync local GRACE/GRACE-FO files with ITSG GRAZ server
def itsg_graz_grace_sync(DIRECTORY, RELEASE=None, LMAX=None, TIMEOUT=0,
    LOG=False, LIST=False, MODE=0o775, CLOBBER=False):

    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: ITSG_GRAZ_GRACE_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'ITSG_GRAZ_GRACE_sync_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'ITSG GRAZ GRACE Sync Log ({today})')
        logging.info(f'Release: {RELEASE}')
        logging.info(f'LMAX: {LMAX:d}')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # ITSG GRAZ server
    HOST = ['http://ftp.tugraz.at','outgoing','ITSG','GRACE']
    # open connection with ITSG GRAZ server at remote directory
    release_directory = f'ITSG-{RELEASE}'
    # regular expression operators for ITSG data and models
    itsg_products = []
    itsg_products.append(r'atmosphere')
    itsg_products.append(r'dealiasing')
    itsg_products.append(r'oceanBottomPressure')
    itsg_products.append(r'ocean')
    itsg_products.append(r'Grace2014')
    itsg_products.append(r'Grace2016')
    itsg_products.append(r'Grace2018')
    itsg_products.append(r'Grace_operational')
    itsg_pattern = (r'(AOD1B_RL\d+|model|ITSG)[-_]({0})(_n\d+)?_'
        r'(\d+)-(\d+)(\.gfc)').format(r'|'.join(itsg_products))
    R1 = re.compile(itsg_pattern, re.VERBOSE | re.IGNORECASE)
    # local directory for release
    DREL = {}
    DREL['Grace2014'] = '2014'
    DREL['Grace2016'] = '2016'
    DREL['Grace2018'] = '2018'
    DREL['Grace_operational'] = '2018'
    # local dealiasing directories for each product
    DEALIASING = {}
    DEALIASING['atmosphere'] = 'GAA'
    DEALIASING['ocean'] = 'GAB'
    DEALIASING['dealiasing'] = 'GAC'
    DEALIASING['oceanBottomPressure'] = 'GAD'

    # sync ITSG GRAZ dealiasing products
    subdir = 'background' if (RELEASE == 'Grace2014') else 'monthly_background'
    REMOTE = [*HOST,release_directory,'monthly',subdir]
    files,mtimes = gravtk.utilities.http_list(REMOTE,
        timeout=TIMEOUT,pattern=R1,sort=True)
    # for each file on the remote directory
    for colname,remote_mtime in zip(files,mtimes):
        # extract parameters from input filename
        PFX,PRD,trunc,year,month,SFX = R1.findall(colname).pop()
        # local directory for output GRAZ data
        local_dir=os.path.join(DIRECTORY,'GRAZ',DREL[RELEASE],DEALIASING[PRD])
        # check if local directory exists and recursively create if not
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        # local and remote versions of the file
        local_file = os.path.join(local_dir,colname)
        remote_file = posixpath.join(*REMOTE,colname)
        # copy file from remote directory comparing modified dates
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER, MODE=MODE)

    # sync ITSG GRAZ data for truncation
    subdir = f'monthly_n{LMAX:d}'
    REMOTE = [*HOST,release_directory,'monthly',subdir]
    files,mtimes = gravtk.utilities.http_list(REMOTE,
        timeout=TIMEOUT,pattern=R1,sort=True)
    # local directory for output GRAZ data
    local_dir = os.path.join(DIRECTORY,'GRAZ',DREL[RELEASE],'GSM')
    # check if local directory exists and recursively create if not
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    # for each file on the remote directory
    for colname,remote_mtime in zip(files,mtimes):
        # local and remote versions of the file
        local_file = os.path.join(local_dir,colname)
        remote_file = posixpath.join(*REMOTE,colname)
        # copy file from remote directory comparing modified dates
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER, MODE=MODE)

    # create index file for GRACE/GRACE-FO L2 Spherical Harmonic Data
    # DATA PRODUCTS (GAC, GAD, GSM, GAA, GAB)
    for ds in ['GAA','GAB','GAC','GAD','GSM']:
        # local directory for exact data product
        local_dir = os.path.join(DIRECTORY,'GRAZ',DREL[RELEASE],ds)
        if not os.access(local_dir,os.F_OK):
            continue
        # find local GRACE files to create index
        grace_files=[fi for fi in os.listdir(local_dir) if R1.match(fi)]
        # outputting GRACE filenames to index
        index_file = os.path.join(local_dir,'index.txt')
        with open(index_file, mode='w', encoding='utf8') as fid:
            for fi in sorted(grace_files):
                print(fi, file=fid)
        # change permissions of index file
        os.chmod(index_file, MODE)

    # close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file,remote_mtime,local_file,
    TIMEOUT=0,LIST=False,CLOBBER=False,MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if os.access(local_file, os.F_OK):
        # check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        # if remote file is newer: overwrite the local file
        if (gravtk.utilities.even(remote_mtime) >
            gravtk.utilities.even(local_mtime)):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    # if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        # Printing files transferred
        logging.info(f'{remote_file} --> ')
        logging.info(f'\t{local_file}{OVERWRITE}\n')
        # if executing copy command (not only printing the files)
        if not LIST:
            # Create and submit request. There are a wide range of exceptions
            # that can be thrown here, including HTTPError and URLError.
            request = gravtk.utilities.urllib2.Request(remote_file)
            response = gravtk.utilities.urllib2.urlopen(request,
                timeout=TIMEOUT)
            # chunked transfer encoding size
            CHUNK = 16 * 1024
            # copy contents to local file using chunked transfer encoding
            # transfer should work properly with ascii and binary data formats
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO and auxiliary data from the
            ITSG GRAZ server
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # ITSG GRAZ releases
    choices = ['Grace2014','Grace2016','Grace2018','Grace_operational']
    parser.add_argument('--release','-r',
        type=str, nargs='+', metavar='DREL',
        default=['Grace2018','Grace_operational'],choices=choices,
        help='GRAZ Data Releases to sync')
    parser.add_argument('--lmax',
        type=int, default=60, choices=[60,96,120],
        help='Maximum degree and order of GRAZ products')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # ITSG_GRAZ_GRACE_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # check internet connection before attempting to run program
    HOST = posixpath.join('http://ftp.tugraz.at')
    if gravtk.utilities.check_connection(HOST):
        # for each ITSG GRAZ release
        for RELEASE in args.release:
            itsg_graz_grace_sync(args.directory, RELEASE=RELEASE,
                LMAX=args.lmax, TIMEOUT=args.timeout, LOG=args.log,
                LIST=args.list, CLOBBER=args.clobber, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
