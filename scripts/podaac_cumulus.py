#!/usr/bin/env python
u"""
podaac_cumulus.py
Written by Tyler Sutterley (12/2022)

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
    -e X, --endpoint X: CMR url endpoint type
    -t X, --timeout X: Timeout in seconds for blocking operations
    --gzip, -G: Compress output GRACE/GRACE-FO Level-2 granules
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
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: added CMR queries for GRACE/GRACE-FO technical notes
        recursively create geocenter directory if not in file system
    Updated 08/2022: moved regular expression function to utilities
        Dynamically select newest version of granules for index
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
import gzip
import time
import shutil
import logging
import argparse
import gravity_toolkit as gravtk

# PURPOSE: sync local GRACE/GRACE-FO files with JPL PO.DAAC AWS S3 bucket
def podaac_cumulus(client, DIRECTORY, PROC=[], DREL=[], VERSION=[],
    AOD1B=False, ENDPOINT='s3', TIMEOUT=None, GZIP=False, LOG=False,
    CLOBBER=False, MODE=None):

    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    # mission shortnames
    shortname = {'grace':'GRAC', 'grace-fo':'GRFO'}
    # datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: PODAAC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        LOGFILE = f'PODAAC_sync_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'PO.DAAC Cumulus Sync Log ({today})')
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # Degree 1 (geocenter) coefficients
    logging.info('Degree 1 Coefficients:')
    # SLR C2,0 and C3,0 coefficients
    logging.info('C2,0 and C3,0 Coefficients:')
    # compile regular expression operator for remote files
    R1 = re.compile(r'TN-13_GEOC_(CSR|GFZ|JPL)_(.*?).txt', re.VERBOSE)
    R2 = re.compile(r'TN-(14)_C30_C20_GSFC_SLR.txt', re.VERBOSE)
    # check if geocenter directory exists and recursively create if not
    local_dir = os.path.join(DIRECTORY,'geocenter')
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    # current time stamp to use for local files
    mtime = time.time()
    # for each processing center (CSR, GFZ, JPL)
    for pr in PROC:
        # for each data release (RL04, RL05, RL06)
        for rl in DREL:
            # for each unique version of data to sync
            for version in set(VERSION):
                # query CMR for product metadata
                urls = gravtk.utilities.cmr_metadata(
                    mission='grace-fo', center=pr, release=rl,
                    version=version, provider='POCLOUD',
                    endpoint=ENDPOINT)

                # TN-13 JPL degree 1 files
                url, = [url for url in urls if R1.search(url)]
                granule = gravtk.utilities.url_split(url)[-1]
                local_file = os.path.join(DIRECTORY,'geocenter',granule)
                # access auxiliary data from endpoint
                if (ENDPOINT == 'data'):
                    http_pull_file(url, mtime, local_file,
                        TIMEOUT=TIMEOUT, CLOBBER=CLOBBER, MODE=MODE)
                elif (ENDPOINT == 's3'):
                    bucket = gravtk.utilities.s3_bucket(url)
                    key = gravtk.utilities.s3_key(url)
                    response = client.get_object(Bucket=bucket, Key=key)
                    s3_pull_file(response, mtime, local_file,
                        CLOBBER=CLOBBER, MODE=MODE)

                # TN-14 SLR C2,0 and C3,0 files
                url, = [url for url in urls if R2.search(url)]
                granule = gravtk.utilities.url_split(url)[-1]
                local_file = os.path.join(DIRECTORY,granule)
                # access auxiliary data from endpoint
                if (ENDPOINT == 'data'):
                    http_pull_file(url, mtime, local_file,
                        TIMEOUT=TIMEOUT, CLOBBER=CLOBBER, MODE=MODE)
                elif (ENDPOINT == 's3'):
                    bucket = gravtk.utilities.s3_bucket(url)
                    key = gravtk.utilities.s3_key(url)
                    response = client.get_object(Bucket=bucket, Key=key)
                    s3_pull_file(response, mtime, local_file,
                        CLOBBER=CLOBBER, MODE=MODE)

    # GRACE/GRACE-FO AOD1B dealiasing products
    if AOD1B:
        logging.info('GRACE L1B Dealiasing Products:')
        # for each data release (RL04, RL05, RL06)
        for rl in DREL:
            # print string of exact data product
            logging.info(f'GFZ/AOD1B/{rl}')
            # local directory for exact data product
            local_dir = os.path.join(DIRECTORY,'AOD1B',rl)
            # check if directory exists and recursively create if not
            if not os.path.exists(local_dir):
                os.makedirs(local_dir,MODE)
            # query CMR for dataset
            ids,urls,mtimes = gravtk.utilities.cmr(
                mission='grace', level='L1B', center='GFZ', release=rl,
                product='AOD1B', start_date='2002-01-01T00:00:00',
                provider='POCLOUD', endpoint=ENDPOINT)
            # for each model id and url
            for id,url,mtime in zip(ids,urls,mtimes):
                # retrieve GRACE/GRACE-FO files
                granule = gravtk.utilities.url_split(url)[-1]
                local_file = os.path.join(local_dir,granule)
                # access data from endpoint
                if (ENDPOINT == 'data'):
                    http_pull_file(url, mtime, local_file,
                        TIMEOUT=TIMEOUT, CLOBBER=CLOBBER, MODE=MODE)
                elif (ENDPOINT == 's3'):
                    bucket = gravtk.utilities.s3_bucket(url)
                    key = gravtk.utilities.s3_key(url)
                    response = client.get_object(Bucket=bucket, Key=key)
                    s3_pull_file(response, mtime, local_file,
                        CLOBBER=CLOBBER, MODE=MODE)

    # GRACE/GRACE-FO level-2 spherical harmonic products
    logging.info('GRACE/GRACE-FO L2 Global Spherical Harmonics:')
    # for each processing center (CSR, GFZ, JPL)
    for pr in PROC:
        # for each data release (RL04, RL05, RL06)
        for rl in DREL:
            # for each level-2 product (GAC, GAD, GSM, GAA, GAB)
            for ds in DSET[pr]:
                # local directory for exact data product
                local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                # check if directory exists and recursively create if not
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir,MODE)
                # list of GRACE/GRACE-FO files for index
                grace_files = []
                # for each satellite mission (grace, grace-fo)
                for i,mi in enumerate(['grace','grace-fo']):
                    # print string of exact data product
                    logging.info(f'{mi} {pr}/{rl}/{ds}')
                    # query CMR for dataset
                    ids,urls,mtimes = gravtk.utilities.cmr(
                        mission=mi, center=pr, release=rl, product=ds,
                        version=VERSION[i], provider='POCLOUD',
                        endpoint=ENDPOINT)
                    # regular expression operator for data product
                    rx = gravtk.utilities.compile_regex_pattern(
                        pr, rl, ds, mission=shortname[mi])
                    # for each model id and url
                    for id,url,mtime in zip(ids,urls,mtimes):
                        # retrieve GRACE/GRACE-FO files
                        granule = gravtk.utilities.url_split(url)[-1]
                        suffix = '.gz' if GZIP else ''
                        local_file = os.path.join(local_dir, f'{granule}{suffix}')
                        # access data from endpoint
                        if (ENDPOINT == 'data'):
                            http_pull_file(url, mtime, local_file,
                                GZIP=GZIP, TIMEOUT=TIMEOUT,
                                CLOBBER=CLOBBER, MODE=MODE)
                        elif (ENDPOINT == 's3'):
                            bucket = gravtk.utilities.s3_bucket(url)
                            key = gravtk.utilities.s3_key(url)
                            response = client.get_object(Bucket=bucket, Key=key)
                            s3_pull_file(response, mtime, local_file,
                                GZIP=GZIP, CLOBBER=CLOBBER, MODE=MODE)
                    # find local GRACE/GRACE-FO files to create index
                    granules = sorted([f for f in os.listdir(local_dir) if rx.match(f)])
                    # reduce list of GRACE/GRACE-FO files to unique dates
                    granules = gravtk.time.reduce_by_date(granules)
                    # extend list of GRACE/GRACE-FO files with granules
                    grace_files.extend(granules)

                # outputting GRACE/GRACE-FO filenames to index
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
def http_pull_file(remote_file, remote_mtime, local_file,
    GZIP=False, TIMEOUT=120, CLOBBER=False, MODE=0o775):
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
        logging.info(f'{remote_file} -->')
        logging.info(f'\t{local_file}{OVERWRITE}\n')
        # chunked transfer encoding size
        CHUNK = 16 * 1024
        # Create and submit request.
        # There are a range of exceptions that can be thrown here
        # including HTTPError and URLError.
        request = gravtk.utilities.urllib2.Request(remote_file)
        response = gravtk.utilities.urllib2.urlopen(request,
            timeout=TIMEOUT)
        # copy remote file contents to local file
        if GZIP:
            with gzip.GzipFile(local_file, 'wb', 9, None, remote_mtime) as f:
                shutil.copyfileobj(response, f)
        else:
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
        # keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        os.chmod(local_file, MODE)

# PURPOSE: pull file from AWS s3 bucket checking if file exists locally
# and if the remote file is newer than the local file
def s3_pull_file(response, remote_mtime, local_file,
    GZIP=False, CLOBBER=False, MODE=0o775):
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
        logging.info(f'{local_file}{OVERWRITE}')
        # chunked transfer encoding size
        CHUNK = 16 * 1024
        # copy remote file contents to local file
        if GZIP:
            with gzip.GzipFile(local_file, 'wb', 9, None, remote_mtime) as f:
                shutil.copyfileobj(response['Body'], f)
        else:
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response['Body'], f, CHUNK)
        # keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO and auxiliary data from the
            NASA JPL PO.DAAC Cumulus AWS bucket.
            """
    )
    # command line parameters
    # NASA Earthdata credentials
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
    # working data directory
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
    # GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs=2,
        default=['0','1'], choices=['0','1','2','3'],
        help='GRACE/GRACE-FO Level-2 data version')
    # GRACE/GRACE-FO dealiasing products
    parser.add_argument('--aod1b','-a',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Level-1B dealiasing products')
    # CMR endpoint type
    parser.add_argument('--endpoint','-e',
        type=str, default='s3', choices=['s3','data'],
        help='CMR url endpoint type')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # output compressed files
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Compress output GRACE/GRACE-FO Level-2 granules')
    # Output log file in form
    # PODAAC_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
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

    # NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    # check internet connection before attempting to run program
    opener = gravtk.utilities.attempt_login(URS,
        username=args.user, password=args.password,
        netrc=args.netrc)

    # Create and submit request to create AWS session
    # There are a range of exceptions that can be thrown here
    # including HTTPError and URLError.
    HOST = 'https://archive.podaac.earthdata.nasa.gov/s3credentials'
    # get aws s3 client object
    client = gravtk.utilities.s3_client(HOST, args.timeout)

    # retrieve data objects from s3 client
    podaac_cumulus(client, args.directory, PROC=args.center,
        DREL=args.release, VERSION=args.version, AOD1B=args.aod1b,
        ENDPOINT=args.endpoint, TIMEOUT=args.timeout,
        GZIP=args.gzip, LOG=args.log, CLOBBER=args.clobber,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
