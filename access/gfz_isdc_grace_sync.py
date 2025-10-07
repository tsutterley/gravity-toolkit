#!/usr/bin/env python
u"""
gfz_isdc_grace_sync.py
Written by Tyler Sutterley (10/2025)
Syncs GRACE/GRACE-FO data from the GFZ Information System and Data Center (ISDC)

CALLING SEQUENCE:
    python gfz_isdc_grace_sync.py

OUTPUTS:
    CSR RL06: GAC/GAD/GSM
    GFZ RL06: GAA/GAB/GAC/GAD/GSM/AOD1b
    JPL RL06: GAA/GAB/GAC/GAD/GSM

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO data releases to sync
    -v X, --version X: GRACE/GRACE-FO Level-2 data version to sync
    -n, --newsletters: sync GRACE/GRACE-FO newsletters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing data in transfer
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2025: switch to https as ftp server is being retired
    Updated 09/2023: don't restrict version number to a set list
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: fix version check for mission
    Updated 08/2022: moved regular expression function to utilities
        Dynamically select newest version of granules for index
    Updated 04/2022: added option for GRACE/GRACE-FO Level-2 data version
        sync GRACE/GRACE-FO technical notes and newsletters
        refactor to always try syncing from both grace and grace-fo missions
        use argparse descriptions within sphinx documentation
    Updated 03/2022: update regular expression pattern for finding files
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 01/2021: using utilities module to list files from ftp
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 03/2020 for public release
    Updated 03/2020: new GFZ ISDC ftp server website
    Updated 09/2019: checksum option to not overwrite existing data files
        added GRACE Follow-On data sync
    Written 08/2018
"""
from __future__ import print_function

import sys
import os
import re
import ssl
import copy
import time
import shutil
import logging
import pathlib
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: sync local GRACE/GRACE-FO files with GFZ ISDC server
def gfz_isdc_grace_sync(DIRECTORY, PROC=[], DREL=[], VERSION=[],
    NEWSLETTERS=False, TIMEOUT=None, LOG=False, LIST=False,
    CLOBBER=False, MODE=None):

    # check if directory exists and recursively create if not
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # GFZ ISDC https host
    HOST = 'https://isdc-data.gfz.de/'
    # mission shortnames
    shortname = {'grace':'GRAC', 'grace-fo':'GRFO'}
    # datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: GFZ_ISDC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = DIRECTORY.joinpath(f'GFZ_ISDC_sync_{today}.log')
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'GFZ ISDC Sync Log ({today})')
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # Degree 1 (geocenter) coefficients
    logging.info('Degree 1 Coefficients:')
    local_dir = DIRECTORY.joinpath('geocenter')
    # check if geocenter directory exists and recursively create if not
    local_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
    # TN-13 JPL degree 1 files
    # compile regular expression operator for remote files
    R1 = re.compile(r'TN-13_GEOC_(CSR|GFZ|JPL)_(.*?).txt$', re.VERBOSE)
    # get filenames from remote directory
    remote_files,remote_mtimes = http_list(
        [HOST,'grace-fo','DOCUMENTS','TECHNICAL_NOTES'],
        timeout=TIMEOUT, pattern=R1, sort=True)
    # for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        # extract filename from regex object
        remote_path = [HOST,'grace-fo','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = local_dir.joinpath(fi)
        http_pull_file(remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, MODE=MODE)

    # SLR C2,0 coefficients
    logging.info('C2,0 Coefficients:')
    # compile regular expression operator for remote files
    R1 = re.compile(r'TN-(05|07|11)_C20_SLR_RL(.*?).txt$', re.VERBOSE)
    # get filenames from remote directory
    remote_files,remote_mtimes = http_list(
        [HOST,'grace','DOCUMENTS','TECHNICAL_NOTES'],
        timeout=TIMEOUT, pattern=R1, sort=True)
    # for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        # extract filename from regex object
        remote_path = [HOST,'grace','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = DIRECTORY.joinpath(re.sub(r'(_RL.*?).txt','.txt',fi))
        http_pull_file(remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, MODE=MODE)

    # SLR C3,0 coefficients
    logging.info('C3,0 Coefficients:')
    # compile regular expression operator for remote files
    R1 = re.compile(r'TN-(14)_C30_C20_SLR_GSFC.txt$', re.VERBOSE)
    # get filenames from remote directory
    remote_files,remote_mtimes = http_list(
        [HOST,'grace-fo','DOCUMENTS','TECHNICAL_NOTES'],
        timeout=TIMEOUT, pattern=R1, sort=True)
    # for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        # extract filename from regex object
        remote_path = [HOST,'grace-fo','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = DIRECTORY.joinpath(re.sub(r'(SLR_GSFC)','GSFC_SLR',fi))
        http_pull_file(remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, MODE=MODE)

    # TN-08 GAE, TN-09 GAF and TN-10 GAG ECMWF atmosphere correction products
    logging.info('TN-08 GAE, TN-09 GAF and TN-10 GAG products:')
    ECMWF_files = []
    ECMWF_files.append('TN-08_GAE-2_2006032-2010031_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-09_GAF-2_2010032-2015131_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-10_GAG-2_2015132-2099001_0000_EIGEN_G---_0005.gz')
    # compile regular expression operator for remote files
    R1 = re.compile(r'({0}|{1}|{2})'.format(*ECMWF_files), re.VERBOSE)
    # get filenames from remote directory
    remote_files,remote_mtimes = http_list(
        [HOST,'grace','DOCUMENTS','TECHNICAL_NOTES'],
        timeout=TIMEOUT, pattern=R1, sort=True)
    # for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        # extract filename from regex object
        remote_path = [HOST,'grace','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = DIRECTORY.joinpath(fi)
        http_pull_file(remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, MODE=MODE)

    # GRACE and GRACE-FO newsletters
    if NEWSLETTERS:
        # local newsletter directory (place GRACE and GRACE-FO together)
        local_dir = DIRECTORY.joinpath('newsletters')
        # check if newsletters directory exists and recursively create if not
        local_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
        # for each satellite mission (grace, grace-fo)
        for i,mi in enumerate(['grace','grace-fo']):
            logging.info(f'{mi} Newsletters:')
            # compile regular expression operator for remote files
            NAME = mi.upper().replace('-','_')
            R1 = re.compile(rf'{NAME}_SDS_NL_(\d+).pdf', re.VERBOSE)
            # find years for GRACE/GRACE-FO newsletters
            years,_  = http_list([HOST,mi,'DOCUMENTS','NEWSLETTER'],
                timeout=TIMEOUT, pattern=r'\d+',
                sort=True)
            # for each year of GRACE/GRACE-FO newsletters
            for Y in years:
                # find GRACE/GRACE-FO newsletters
                remote_files,remote_mtimes = http_list(
                    [HOST,mi,'DOCUMENTS','NEWSLETTER',Y],
                    timeout=TIMEOUT, pattern=R1,
                    sort=True)
                # for each file on the remote server
                for fi,remote_mtime in zip(remote_files,remote_mtimes):
                    # extract filename from regex object
                    remote_path = [HOST,mi,'DOCUMENTS','NEWSLETTER',Y,fi]
                    local_file = local_dir.joinpath(fi)
                    http_pull_file(remote_path, remote_mtime,
                        local_file, TIMEOUT=TIMEOUT, LIST=LIST,
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
                local_dir = DIRECTORY.joinpath(pr, rl, ds)
                # check if directory exists and recursively create if not
                local_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
                # list of GRACE/GRACE-FO files for index
                grace_files = []
                # for each satellite mission (grace, grace-fo)
                for i,mi in enumerate(['grace','grace-fo']):
                    # modifiers for intermediate data releases
                    if (int(VERSION[i]) > 0):
                        drel_str = f'{rl}.{VERSION[i]}'
                    else:
                        drel_str = copy.copy(rl)
                    # print string of exact data product
                    logging.info(f'{mi}/{pr}/{drel_str}/{ds}')
                    # compile the regular expression operator to find files
                    R1 = re.compile(rf'({ds}-(.*?)(gz|txt|dif))')
                    # get filenames from remote directory
                    remote_files,remote_mtimes = http_list(
                        [HOST,mi,'Level-2',pr,drel_str], timeout=TIMEOUT,
                        pattern=R1, sort=True)
                    for fi,remote_mtime in zip(remote_files,remote_mtimes):
                        # extract filename from regex object
                        remote_path = [HOST,mi,'Level-2',pr,drel_str,fi]
                        local_file = local_dir.joinpath(fi)
                        http_pull_file(remote_path, remote_mtime,
                            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
                            CLOBBER=CLOBBER, MODE=MODE)
                    # regular expression operator for data product
                    rx = gravtk.utilities.compile_regex_pattern(
                        pr, rl, ds, mission=shortname[mi])
                    # find local GRACE/GRACE-FO files to create index
                    granules = sorted([f.name for f in local_dir.iterdir()
                        if rx.match(f.name)])
                    # reduce list of GRACE/GRACE-FO files to unique dates
                    granules = gravtk.time.reduce_by_date(granules)
                    # extend list of GRACE/GRACE-FO files with granules
                    grace_files.extend(granules)

                # outputting GRACE/GRACE-FO filenames to index
                index_file = local_dir.joinpath('index.txt')
                with index_file.open(mode='w', encoding='utf8') as fid:
                    for fi in sorted(grace_files):
                        print(fi, file=fid)
                # change permissions of index file
                index_file.chmod(mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

# PURPOSE: list a directory on the GFZ https server
def http_list(
        HOST: str | list,
        timeout: int | None = None,
        context: ssl.SSLContext = gravtk.utilities._default_ssl_context,
        pattern: str | re.Pattern = '',
        sort: bool = False
    ):
    """
    List a directory on the GFZ https Server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = gravtk.utilities.url_split(HOST)
    # regular expression pattern for finding files and modification times
    parser = r'\<a\shref=.*?\>(.*?)\<\/a\>\s+(\d{4}-\d{2}-\d{2}\s+\d{2}\:\d{2})'
    rx = re.compile(parser, re.VERBOSE)
    # try listing from http
    try:
        # Create and submit request.
        request = gravtk.utilities.urllib2.Request(posixpath.join(*HOST))
        response = gravtk.utilities.urllib2.urlopen(request,
            timeout=timeout, context=context)
    except Exception as exc:
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    # read the directory listing
    contents = response.readlines()
    # read and parse request for files (column names and modified times)
    lines = [l for l in contents if rx.search(l.decode('utf-8'))]
    # column names and last modified times
    colnames = [None]*len(lines)
    collastmod = [None]*len(lines)
    for i, l in enumerate(lines):
        colnames[i], lastmod = rx.findall(l.decode('utf-8')).pop()
        # get the Unix timestamp value for a modification time
        collastmod[i] = gravtk.utilities.get_unix_time(lastmod,
            format='%Y-%m-%d %H:%M')
    # reduce using regular expression pattern
    if pattern:
        i = [i for i,f in enumerate(colnames) if re.search(pattern, f)]
        # reduce list of column names and last modified times
        colnames = [colnames[indice] for indice in i]
        collastmod = [collastmod[indice] for indice in i]
    # sort the list
    if sort:
        i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
        # sort list of column names and last modified times
        colnames = [colnames[indice] for indice in i]
        collastmod = [collastmod[indice] for indice in i]
    # return the list of column names and last modified times
    return (colnames, collastmod)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_path, remote_mtime, local_file,
    TIMEOUT=0, LIST=False, CLOBBER=False, MODE=0o775):
    # verify inputs for remote http host
    if isinstance(remote_path, str):
        remote_path = gravtk.utilities.url_split(remote_path)
    # construct remote file path
    remote_file = posixpath.join(*remote_path)
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    local_file = pathlib.Path(local_file).expanduser().absolute()
    if local_file.exists():
        # check last modification time of local file
        local_mtime = local_file.stat().st_mtime
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
        logging.info(f'\t{str(local_file)}{OVERWRITE}\n')
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
            with local_file.open(mode='wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (local_file.stat().st_atime, remote_mtime))
            local_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO data from the GFZ
            Information System and Data Center (ISDC)
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'], choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO processing center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'], choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    # GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs=2,
        default=['0','1'],
        help='GRACE/GRACE-FO Level-2 data version')
    # GRACE/GRACE-FO newsletters
    parser.add_argument('--newsletters','-n',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Newsletters')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # GFZ_ISDC_sync_2002-04-01.log
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

    # GFZ ISDC https host
    HOST = 'https://isdc-data.gfz.de/'
    # check internet connection before attempting to run program
    if gravtk.utilities.check_connection(HOST):
        gfz_isdc_grace_sync(args.directory, PROC=args.center,
            DREL=args.release, VERSION=args.version,
            NEWSLETTERS=args.newsletters, TIMEOUT=args.timeout,
            LIST=args.list, LOG=args.log, CLOBBER=args.clobber,
            MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

# run main program
if __name__ == '__main__':
    main()
