#!/usr/bin/env python
u"""
podaac_cumulus.py
Written by Tyler Sutterley (04/2022)

Syncs GRACE/GRACE-FO data from NASA JPL PO.DAAC Cumulus AWS S3 bucket
S3 Cumulus syncs are only available in AWS instances in us-west-2

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

CALLING SEQUENCE:
    python podaac_cumulus.py --user <username>
    where <username> is your NASA Earthdata username

OUTPUTS:
    CSR RL06: GAC/GAD/GSM
    GFZ RL06: GAA/GAB/GAC/GAD/GSM
    JPL RL06: GAA/GAB/GAC/GAD/GSM
    GFZ RL06: Level-1b dealiasing solutions

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases to sync
    -v X, --version X: GRACE/GRACE-FO Level-2 Data Version to sync
    -a, --aod1b: sync GRACE/GRACE-FO Level-1B dealiasing products
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing data in transfer
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    boto3: Amazon Web Services (AWS) SDK for Python
        https://boto3.amazonaws.com/v1/documentation/api/latest/index.html
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 04/2022: added option for GRACE/GRACE-FO Level-2 data version
        refactor to always try syncing from both grace and grace-fo missions
        use granule identifiers from CMR query to build output file index
        use argparse descriptions within sphinx documentation
    Written 03/2022 with release of PO.DAAC Cumulus
"""
from __future__ import print_function

import sys
import os
import re
import time
import shutil
import logging
import argparse
import gravity_toolkit.utilities

#-- PURPOSE: create and compile regular expression operator to find GRACE files
def compile_regex_pattern(PROC, DREL, DSET, version='0'):
    if ((DSET == 'GSM') and (PROC == 'CSR') and (DREL in ('RL04','RL05'))):
        #-- CSR GSM: only monthly degree 60 products
        #-- not the longterm degree 180, degree 96 dataset or the
        #-- special order 30 datasets for the high-resonance months
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        regex_pattern=r'{0}-2_\d+-\d+_\d+_UTCSR_0060_000{1:d}(\.gz)?$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'CSR') and (DREL == 'RL06')):
        #-- CSR GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', release.zfill(2), version.zfill(2))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_UTCSR_{2}_{3}{4}(\.gz)?$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL04')):
        #-- GFZ RL04: only unconstrained solutions (not GK2 products)
        regex_pattern=r'{0}-2_\d+-\d+_\d+_EIGEN_G---_0004(\.gz)?$'.format(DSET)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL05')):
        #-- GFZ RL05: updated RL05a products which are less constrained to
        #-- the background model.  Allow regularized fields
        regex_unconst=r'{0}-2_\d+-\d+_\d+_EIGEN_G---_005a(\.gz)?$'.format(DSET)
        regex_regular=r'{0}-2_\d+-\d+_\d+_EIGEN_GK2-_005a(\.gz)?$'.format(DSET)
        regex_pattern=r'{0}|{1}'.format(regex_unconst,regex_regular)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL06')):
        #-- GFZ GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', release.zfill(2), version.zfill(2))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_GFZOP_{2}_{3}{4}(\.gz)?$' .format(*args)
    elif (PROC == 'JPL') and DREL in ('RL04','RL05'):
        #-- JPL: RL04a and RL05a products (denoted by 0001)
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        regex_pattern=r'{0}-2_\d+-\d+_\d+_JPLEM_0001_000{1:d}(\.gz)?$'.format(*args)
    elif ((DSET == 'GSM') and (PROC == 'JPL') and (DREL == 'RL06')):
        #-- JPL GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', release.zfill(2), version.zfill(2))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_JPLEM_{2}_{3}{4}(\.gz)?$' .format(*args)
    else:
        regex_pattern=r'{0}-2_([a-zA-Z0-9_\-]+)(\.gz)?$'.format(DSET)
    #-- return the compiled regular expression operator used to find files
    return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local GRACE/GRACE-FO files with JPL PO.DAAC AWS S3 bucket
def podaac_cumulus(client, DIRECTORY, PROC=[], DREL=[], VERSION=[],
    AOD1B=False, LOG=False, CLOBBER=False, MODE=None):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- PO.DAAC cumulus bucket
    bucket = 'podaac-ops-cumulus-protected'
    #-- datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: PODAAC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'PODAAC_sync_{0}.log'.format(today)
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('PO.DAAC Cumulus Sync Log ({0})'.format(today))
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)


    #-- GRACE/GRACE-FO AOD1B dealiasing products
    if AOD1B:
        logging.info('GRACE L1B Dealiasing Products:')
        #-- for each data release (RL04, RL05, RL06)
        for rl in DREL:
            #-- print string of exact data product
            logging.info('{0}/{1}/{2}'.format('GFZ','AOD1B',rl))
            #-- local directory for exact data product
            local_dir = os.path.join(DIRECTORY,'AOD1B',rl)
            #-- check if directory exists and recursively create if not
            if not os.path.exists(local_dir):
                os.makedirs(local_dir,MODE)
            #-- query CMR for dataset
            ids,urls,mtimes = gravity_toolkit.utilities.cmr(
                mission='grace', level='L1B', center='GFZ', release=rl,
                product='AOD1B', start_date='2002-01-01T00:00:00',
                provider='POCLOUD', endpoint='s3')
            #-- for each model id and url
            for id,url in zip(ids,urls):
                #-- retrieve GRACE/GRACE-FO files
                key = gravity_toolkit.utilities.s3_key(url)
                granule = gravity_toolkit.utilities.url_split(url)[-1]
                response = client.get_object(Bucket=bucket, Key=key)
                s3_pull_file(response, os.path.join(local_dir,granule),
                    CLOBBER=CLOBBER, MODE=MODE)

    #-- GRACE/GRACE-FO level-2 spherical harmonic products
    logging.info('GRACE/GRACE-FO L2 Global Spherical Harmonics:')
    #-- for each processing center (CSR, GFZ, JPL)
    for pr in PROC:
        #-- for each data release (RL04, RL05, RL06)
        for rl in DREL:
            #-- for each level-2 product (GAC, GAD, GSM, GAA, GAB)
            for ds in DSET[pr]:
                #-- local directory for exact data product
                local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                #-- check if directory exists and recursively create if not
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir,MODE)
                #-- list of GRACE/GRACE-FO files for index
                grace_files = []
                #-- for each satellite mission (grace, grace-fo)
                for i,mi in enumerate(['grace','grace-fo']):
                    #-- print string of exact data product
                    logging.info('{0} {1}/{2}/{3}'.format(mi, pr, rl, ds))
                    #-- query CMR for dataset
                    ids,urls,mtimes = gravity_toolkit.utilities.cmr(
                        mission=mi, center=pr, release=rl, product=ds,
                        version=VERSION[i], provider='POCLOUD', endpoint='s3')
                    #-- regular expression operator for data product
                    rx = compile_regex_pattern(pr, rl, ds, version=VERSION[i])
                    #-- for each model id and url
                    for id,url,mtime in zip(ids,urls,mtimes):
                        #-- retrieve GRACE/GRACE-FO files
                        key = gravity_toolkit.utilities.s3_key(url)
                        response = client.get_object(Bucket=bucket, Key=key)
                        granule = gravity_toolkit.utilities.url_split(url)[-1]
                        local_file = os.path.join(local_dir, granule)
                        s3_pull_file(response, mtime, local_file,
                            CLOBBER=CLOBBER, MODE=MODE)
                        #-- extend list of GRACE/GRACE-FO files with granule
                        grace_files.append(granule) if rx.match(granule) else None

                #-- outputting GRACE/GRACE-FO filenames to index
                with open(os.path.join(local_dir,'index.txt'),'w') as fid:
                    for fi in sorted(grace_files):
                        print('{0}'.format(fi), file=fid)
                #-- change permissions of index file
                os.chmod(os.path.join(local_dir,'index.txt'), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from AWS s3 bucket checking if file exists locally
#-- and if the remote file is newer than the local file
def s3_pull_file(response, remote_mtime, local_file, CLOBBER=False, MODE=0o775):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if os.access(local_file, os.F_OK):
        #-- check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (gravity_toolkit.utilities.even(remote_mtime) >
            gravity_toolkit.utilities.even(local_mtime)):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- Printing files transferred
        logging.info('{0}{1}'.format(local_file, OVERWRITE))
        #-- chunked transfer encoding size
        CHUNK = 16 * 1024
        #-- copy remote file contents to local file
        with open(local_file, 'wb') as f:
            shutil.copyfileobj(response['Body'], f, CHUNK)
        #-- keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        os.chmod(local_file, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO and auxiliary data from the
            NASA JPL PO.DAAC Cumulus AWS bucket.
            """
    )
    #-- command line parameters
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- GRACE/GRACE-FO processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'], choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO processing center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'], choices=['RL06'],
        help='GRACE/GRACE-FO data release')
    #-- GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs='+',
        default=['0','1'], choices=['0','1','2','3'],
        help='GRACE/GRACE-FO Level-2 data version')
    #-- GRACE/GRACE-FO dealiasing products
    parser.add_argument('--aod1b','-a',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Level-1B dealiasing products')
    #-- connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    #-- Output log file in form
    #-- PODAAC_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    #-- check internet connection before attempting to run program
    opener = gravity_toolkit.utilities.attempt_login(URS,
        username=args.user, password=args.password,
        netrc=args.netrc)

    #-- Create and submit request to create AWS session
    #-- There are a range of exceptions that can be thrown here
    #-- including HTTPError and URLError.
    HOST = 'https://archive.podaac.earthdata.nasa.gov/s3credentials'
    #-- get aws s3 client object
    client = gravity_toolkit.utilities.s3_client(HOST, args.timeout)
    #-- retrieve data objects from s3 client
    podaac_cumulus(client, args.directory, PROC=args.center,
        DREL=args.release, VERSION=args.version, AOD1B=args.aod1b,
        LOG=args.log, CLOBBER=args.clobber, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
