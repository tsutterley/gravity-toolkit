#!/usr/bin/env python
u"""
gfz_icgem_costg_ftp.py
Written by Tyler Sutterley (12/2022)
Syncs GRACE/GRACE-FO/Swarm COST-G data from the GFZ International
    Centre for Global Earth Models (ICGEM)

CALLING SEQUENCE:
    python gfz_icgem_costg_ftp.py --mission Grace Grace-FO Swarm

OUTPUTS:
    Grace: COST-G GRACE spherical harmonic fields
    Grace-FO: COST-G GRACE-FO spherical harmonic fields
    Swarm: COST-G Swarm spherical harmonic fields

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -m X, --mission X: Sync GRACE, GRACE Follow-On (GRACE-FO) and Swarm data
         Grace
         Grace-FO
         Swarm
    -r X, --release X: Data release to sync
    -t X, --timeout X: Timeout in seconds for blocking operations
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

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
import ftplib
import shutil
import hashlib
import logging
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: create and compile regular expression operator to find files
def compile_regex_pattern(MISSION, DSET):
    if ((DSET == 'GSM') and (MISSION == 'Swarm')):
        # regular expression operators for Swarm data
        regex=r'(SW)_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)(\.gfc|\.ZIP)'
    elif ((DSET != 'GSM') and (MISSION == 'Swarm')):
        regex=r'(GAA|GAB|GAC|GAD)_Swarm_(\d+)_(\d{2})_(\d{4})(\.gfc|\.ZIP)'
    else:
        regex=rf'{DSET}-2_(.*?)\.gfc$'
    # return the compiled regular expression operator used to find files
    return re.compile(regex, re.VERBOSE)

# PURPOSE: sync local GRACE/GRACE-FO/Swarm files with GFZ ICGEM server
def gfz_icgem_costg_ftp(DIRECTORY, MISSION=[], RELEASE=None, TIMEOUT=None,
    LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=None):

    # connect and login to GFZ ICGEM ftp server
    ftp = ftplib.FTP('icgem.gfz-potsdam.de', timeout=TIMEOUT)
    ftp.login()

    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    # dealiasing datasets for each mission
    DSET = {}
    DSET['Grace'] = ['GAC','GSM']
    DSET['Grace-FO'] = ['GSM']
    DSET['Swarm'] = ['GAA','GAB','GAC','GAD','GSM']
    # local subdirectory for data
    LOCAL = {}
    LOCAL['Grace'] = 'COSTG'
    LOCAL['Grace-FO'] = 'COSTG'
    LOCAL['Swarm'] = 'Swarm'

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: GFZ_ICGEM_COST-G_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'GFZ_ICGEM_COST-G_sync_{today}.log'
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info(f'GFZ ICGEM COST-G Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # find files for a particular mission
    logging.info(f'{MISSION} Spherical Harmonics:')

    # Sync gravity field dealiasing products
    for ds in DSET[MISSION]:
        # print string of exact data product
        logging.info(f'{MISSION}/{RELEASE}/{ds}')
        # local directory for exact data product
        local_dir = os.path.join(DIRECTORY,LOCAL[MISSION],RELEASE,ds)
        # check if directory exists and recursively create if not
        if not os.path.exists(local_dir):
            os.makedirs(local_dir,MODE)
        # compile the regular expression operator to find files
        R1 = compile_regex_pattern(MISSION, ds)
        # set the remote path to download files
        if ds in ('GAA','GAB','GAC','GAD') and (MISSION == 'Swarm'):
            remote_path = [ftp.host,'02_COST-G',MISSION,'GAX_products',ds]
        elif ds in ('GAA','GAB','GAC','GAD') and (MISSION != 'Swarm'):
            remote_path = [ftp.host,'02_COST-G',MISSION,'GAX_products']
        elif (MISSION == 'Swarm'):
            remote_path = [ftp.host,'02_COST-G',MISSION,'40x40']
        elif (MISSION == 'Grace'):
            remote_path = [ftp.host,'02_COST-G',MISSION,'unfiltered']
        elif (MISSION == 'Grace-FO'):
            remote_path = [ftp.host,'02_COST-G',MISSION]
        # get filenames from remote directory
        remote_files,remote_mtimes = gravtk.utilities.ftp_list(
            remote_path, timeout=TIMEOUT, basename=True, pattern=R1,
            sort=True)
        # download the file from the ftp server
        for fi,remote_mtime in zip(remote_files,remote_mtimes):
            # remote and local versions of the file
            remote_path.append(fi)
            local_file = os.path.join(local_dir,fi)
            ftp_mirror_file(ftp, remote_path, remote_mtime,
                local_file, TIMEOUT=TIMEOUT, LIST=LIST,
                CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)
            # remove the file from the remote path list
            remote_path.remove(fi)
        # find local GRACE/GRACE-FO/Swarm files to create index
        grace_files=[fi for fi in os.listdir(local_dir) if R1.match(fi)]
        # write each file to an index
        index_file = os.path.join(local_dir,'index.txt')
        with open(index_file, mode='w', encoding='utf8') as fid:
            # output GRACE/GRACE-FO/Swarm filenames to index
            for fi in sorted(grace_files):
                print(fi, file=fid)
        # change permissions of index file
        os.chmod(index_file, MODE)

    # close the ftp connection
    ftp.quit()
    # close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def ftp_mirror_file(ftp,remote_path,remote_mtime,local_file,
    TIMEOUT=None,LIST=False,CLOBBER=False,CHECKSUM=False,MODE=0o775):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        # generate checksum hash for local file
        # open the local_file in binary read mode
        with open(local_file, 'rb') as local_buffer:
            local_hash = hashlib.md5(local_buffer.read()).hexdigest()
        # copy remote file contents to bytesIO object
        remote_buffer = gravtk.utilities.from_ftp(remote_path,
            timeout=TIMEOUT)
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = f' (checksums: {local_hash} {remote_hash})'
    elif os.access(local_file, os.F_OK):
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
        remote_ftp_url = posixpath.join('ftp://',*remote_path)
        logging.info(f'{remote_ftp_url} -->')
        logging.info(f'\t{local_file}{OVERWRITE}\n')
        # if executing copy command (not only printing the files)
        if not LIST:
            # copy file from ftp server or from bytesIO object
            if CHECKSUM and os.access(local_file, os.F_OK):
                # store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, 16 * 1024)
            else:
                # path to remote file
                remote_file = posixpath.join(*remote_path[1:])
                # copy remote file contents to local file
                with open(local_file, 'wb') as f:
                    ftp.retrbinary(f'RETR {remote_file}', f.write)
            # keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO/Swarm COST-G data from the
            GFZ International Centre for Global Earth Models (ICGEM)
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # mission (GRACE, GRACE Follow-On or Swarm)
    choices = ['Grace','Grace-FO','Swarm']
    parser.add_argument('--mission','-m',
        type=str, nargs='+',
        default=['Grace','Grace-FO','Swarm'], choices=choices,
        help='Mission to sync between GRACE, GRACE-FO and Swarm')
    # data release
    parser.add_argument('--release','-r',
        type=str, default='RL01', choices=['RL01'],
        help='Data release to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # GFZ_ICGEM_COST-G_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
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
    HOST = 'icgem.gfz-potsdam.de'
    if gravtk.utilities.check_ftp_connection(HOST):
        for m in args.mission:
            gfz_icgem_costg_ftp(args.directory, MISSION=m,
                RELEASE=args.release, TIMEOUT=args.timeout,
                LIST=args.list, LOG=args.log, CLOBBER=args.clobber,
                CHECKSUM=args.checksum, MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

# run main program
if __name__ == '__main__':
    main()
