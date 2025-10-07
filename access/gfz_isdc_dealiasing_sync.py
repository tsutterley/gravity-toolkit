#!/usr/bin/env python
u"""
gfz_isdc_dealiasing_sync.py
Written by Tyler Sutterley (10/2025)
Syncs GRACE Level-1b dealiasing products from the GFZ Information
    System and Data Center (ISDC)

Optionally outputs as monthly tar files

CALLING SEQUENCE:
    python gfz_isdc_dealiasing_sync.py --year=2015 --release=RL06 --tar

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -r X, --release X: GRACE/GRACE-FO Data Releases to run (RL05,RL06)
    -Y X, --year X: Years of data to sync
    -m X, --month X: Months of data to sync
    -T, --tar: Output data as monthly tar files (.tar.gz or .tgz)
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing data in transfers
    -M X, --mode X: permissions mode of files synced

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2025: switch to https as ftp server is being retired
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: increase default year range to sync
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: added option to sync only specific months
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 01/2021: using utilities module to list files from ftp
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 03/2020: new GFZ ISDC ftp server website
    Updated 06/2019: different suffix with GRACE/GRACE-FO release 6
    Updated 03/2018: made tar file creation optional with --tar
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import ssl
import time
import shutil
import logging
import pathlib
import tarfile
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: syncs GRACE Level-1b dealiasing products from the GFZ data server
# and optionally outputs as monthly tar files
def gfz_isdc_dealiasing_sync(base_dir, DREL, YEAR=None, MONTHS=None, TAR=False,
    TIMEOUT=None, LOG=False, CLOBBER=False, MODE=None):
    # check if directory exists and recursively create if not
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    grace_dir = base_dir.joinpath('AOD1B',DREL)
    grace_dir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: GFZ_AOD1B_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = base_dir.joinpath(f'GFZ_AOD1B_sync_{today}.log')
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'GFZ AOD1b Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # GFZ ISDC https host
    HOST = 'https://isdc-data.gfz.de/'

    # compile regular expression operator for years to sync
    if YEAR is None:
        regex_years = r'\d{4}'
    else:
        regex_years = r'|'.join(rf'{y:d}' for y in YEAR)
    # compile regular expression operator for years to sync
    R1 = re.compile(rf'({regex_years})', re.VERBOSE)
    # suffix for each data release
    SUFFIX = dict(RL04='tar.gz', RL05='tar.gz', RL06='tgz')

    # find remote yearly directories for DREL
    YRS,_ = http_list([HOST,'grace','Level-1B', 'GFZ','AOD',DREL],
        timeout=TIMEOUT, basename=True, pattern=R1, sort=True)
    # for each year
    for Y in YRS:
        # for each month of interest
        for M in MONTHS:
            # output tar file for year and month
            args = (Y, M, DREL.replace('RL',''), SUFFIX[DREL])
            FILE = 'AOD1B_{0}-{1:02d}_{2}.{3}'.format(*args)
            # check if output tar file exists (if TAR)
            local_tar_file = grace_dir.joinpath(FILE)
            TEST = not local_tar_file.exists()
            # compile regular expressions operators for file dates
            # will extract year and month and calendar day from the ascii file
            regex_pattern = r'AOD1B_({0})-({1:02d})-(\d+)_X_\d+.asc.gz$'
            R2 = re.compile(regex_pattern.format(Y,M), re.VERBOSE)
            remote_files,remote_mtimes = http_list(
                [HOST,'grace','Level-1B','GFZ','AOD',DREL,Y],
                timeout=TIMEOUT, basename=True, pattern=R2, sort=True)
            file_count = len(remote_files)
            # if compressing into monthly tar files
            if TAR and (file_count > 0) and (TEST or CLOBBER):
                # copy each gzip file and store within monthly tar files
                tar = tarfile.open(name=local_tar_file, mode='w:gz')
                for fi,remote_mtime in zip(remote_files,remote_mtimes):
                    # remote version of each input file
                    remote = [HOST,'grace','Level-1B','GFZ','AOD',DREL,Y,fi]
                    logging.info(posixpath.join(*remote))
                    # retrieve bytes from remote file
                    remote_buffer = gravtk.utilities.from_sync(remote,
                        timeout=TIMEOUT)
                    # add file to tar
                    tar_info = tarfile.TarInfo(name=fi)
                    tar_info.mtime = remote_mtime
                    tar_info.size = remote_buffer.getbuffer().nbytes
                    tar.addfile(tarinfo=tar_info, fileobj=remote_buffer)
                # close tar file and set permissions level to MODE
                tar.close()
                logging.info(f' --> {local_tar_file}\n')
                local_tar_file.chmod(mode=MODE)
            elif (file_count > 0) and not TAR:
                # copy each gzip file and keep as individual daily files
                for fi,remote_mtime in zip(remote_files,remote_mtimes):
                    # remote and local version of each input file
                    remote = [HOST,'grace','Level-1B','GFZ','AOD',DREL,Y,fi]
                    local_file = grace_dir.joinpath(fi)
                    http_pull_file(remote,remote_mtime,local_file,
                        CLOBBER=CLOBBER, MODE=MODE)

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
        description="""Syncs GRACE Level-1b dealiasing products from
            the GFZ Information System and Data Center (ISDC)
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'], choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    # years to download
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,2023),
        help='Years of data to sync')
    # months to download
    parser.add_argument('--month','-m',
        type=int, nargs='+', default=range(1,13),
        help='Months of data to sync')
    # output dealiasing files as monthly tar files
    parser.add_argument('--tar','-T',
        default=False, action='store_true',
        help='Output data as monthly tar files')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # GFZ_AOD1B_sync_2002-04-01.log
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

    # GFZ ISDC https host
    HOST = 'https://isdc-data.gfz.de/'
    # check internet connection before attempting to run program
    if gravtk.utilities.check_connection(HOST):
        for DREL in args.release:
            gfz_isdc_dealiasing_sync(args.directory, DREL=DREL,
                YEAR=args.year, MONTHS=args.month, TAR=args.tar,
                TIMEOUT=args.timeout, LOG=args.log,
                CLOBBER=args.clobber, MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

# run main program
if __name__ == '__main__':
    main()
