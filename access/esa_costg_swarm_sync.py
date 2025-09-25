#!/usr/bin/env python
u"""
esa_costg_swarm_sync.py
Written by Tyler Sutterley (05/2023)
Syncs Swarm gravity field products from the ESA Swarm Science Server
    https://earth.esa.int/eogateway/missions/swarm/data
    https://www.esa.int/Applications/Observing_the_Earth/Swarm

CALLING SEQUENCE:
    python esa_costg_swarm_sync.py

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -r X, --release X: Data release to sync
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -L, --list: print files to be transferred, but do not execute transfer
    -C, --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode=X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Written 09/2021
"""
from __future__ import print_function

import sys
import re
import os
import io
import json
import time
import shutil
import logging
import pathlib
import argparse
import posixpath
import lxml.etree
import gravity_toolkit as gravtk

# PURPOSE: sync local Swarm files with ESA server
def esa_costg_swarm_sync(DIRECTORY, RELEASE=None, TIMEOUT=None, LOG=False,
    LIST=False, CLOBBER=False, CHECKSUM=False, MODE=0o775):

    # check if directory exists and recursively create if not
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    # local directory for exact data product
    local_dir = DIRECTORY.joinpath('Swarm',RELEASE,'GSM')
    local_dir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: ESA_Swarm_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = DIRECTORY.joinpath(f'ESA_Swarm_sync_{today}.log')
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'ESA Swarm Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # Swarm Science Server url
    # using the JSON api protocols to retrieve files
    # static site is no longer available
    HOST = 'https://swarm-diss.eo.esa.int'
    # compile xml parsers for lxml
    XMLparser = lxml.etree.XMLParser()
    # create "opener" (OpenerDirector instance)
    gravtk.utilities.build_opener(None, None,
        authorization_header=False, urs=HOST)
    # All calls to urllib2.urlopen will now use handler
    # Make sure not to include the protocol in with the URL, or
    # HTTPPasswordMgrWithDefaultRealm will be confused.

    # compile regular expression operator for files
    swarm_data = r'(SW)_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)(\.gfc|\.ZIP)'
    R1 = re.compile(swarm_data, re.VERBOSE)

    # create combined list of filenames and last modified times
    colnames = []
    collastmod = []
    # position, maximum number of files to list, flag to check if done
    pos,maxfiles,prevmax = (0,500,500)
    # iterate to get a compiled list of files
    # will iterate until there are no more files to add to the lists
    while (maxfiles == prevmax):
        # set previous flag to maximum
        prevmax = maxfiles
        # open connection with Swarm science server at remote directory
        # to list maxfiles number of files at position
        parameters = gravtk.utilities.urlencode({'maxfiles':prevmax,
            'pos':pos,'file':posixpath.join('swarm','Level2longterm','EGF')})
        url=posixpath.join(HOST,f'?do=list&{parameters}')
        request = gravtk.utilities.urllib2.Request(url=url)
        response = gravtk.utilities.urllib2.urlopen(request,
            timeout=TIMEOUT)
        table = json.loads(response.read().decode())
        # extend lists with new files
        colnames.extend([t['name'] for t in table['results']])
        collastmod.extend([t['mtime'] for t in table['results']])
        # update maximum number of files
        maxfiles = len(table['results'])
        # update position
        pos += maxfiles

    # find lines of valid files
    valid_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
    # write each file to an index
    index_file = local_dir.joinpath(local_dir,'index.txt')
    fid = index_file.open(mode='w', encoding='utf8')
    # for each data and header file
    for i in valid_lines:
        # remote and local versions of the file
        parameters = gravtk.utilities.urlencode({'file':
            posixpath.join('swarm','Level2longterm','EGF',colnames[i])})
        remote_file = posixpath.join(HOST,
            f'?do=download&{parameters}')
        local_file = local_dir.joinpath(colnames[i])
        # check that file is not in file system unless overwriting
        http_pull_file(remote_file, collastmod[i], local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)
        # output Swarm filenames to index
        print(colnames[i], file=fid)
    # change permissions of index file
    index_file.chmod(mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file, TIMEOUT=120,
    LIST=False, CLOBBER=False, CHECKSUM=False, MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    local_file = pathlib.Path(local_file).expanduser().absolute()
    if CHECKSUM and local_file.exists():
        # generate checksum hash for local file
        # open the local_file in binary read mode
        local_hash = gravtk.utilities.get_hash(local_file)
        # Create and submit request.
        # There are a wide range of exceptions that can be thrown here
        # including HTTPError and URLError.
        req = gravtk.utilities.urllib2.Request(remote_file)
        resp = gravtk.utilities.urllib2.urlopen(req,timeout=TIMEOUT)
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(resp.read())
        remote_buffer.seek(0)
        # generate checksum hash for remote file
        remote_hash = gravtk.utilities.get_hash(remote_buffer)
        # compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = f' (checksums: {local_hash} {remote_hash})'
    elif local_file.exists():
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
            # chunked transfer encoding size
            CHUNK = 16 * 1024
            # copy bytes or transfer file
            if CHECKSUM and local_file.exists():
                # store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with local_file.open(mode='wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
            else:
                # Create and submit request.
                # There are a range of exceptions that can be thrown here
                # including HTTPError and URLError.
                request = gravtk.utilities.urllib2.Request(remote_file)
                response = gravtk.utilities.urllib2.urlopen(request,
                    timeout=TIMEOUT)
                # copy remote file contents to local file
                with local_file.open(mode='wb') as f:
                    shutil.copyfileobj(response, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(local_file, (local_file.stat().st_atime, remote_mtime))
            local_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs Swarm gravity field products from the
            ESA Swarm Science Server
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # data release
    parser.add_argument('--release','-r',
        type=str, default='RL01', choices=['RL01'],
        help='Data release to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # ESA_Swarm_sync_2002-04-01.log
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
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
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
    args = parser.parse_args()

    # check internet connection before attempting to run program
    HOST = 'https://swarm-diss.eo.esa.int'
    if gravtk.utilities.check_connection(HOST):
        esa_costg_swarm_sync(args.directory, RELEASE=args.release,
            TIMEOUT=args.timeout, LOG=args.log, LIST=args.list,
            CLOBBER=args.clobber, CHECKSUM=args.checksum, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
