#!/usr/bin/env python
u"""
cnes_grace_sync.py
Written by Tyler Sutterley (12/2022)

CNES/GRGS GRACE data download program for gravity field products
    https://grace.obs-mip.fr/

Downloads the tar file containing the CNES/GRGS GRACE data for a given release
Iterates through files to determine if any are not in the local file system
For any file to sync: copies from the tar file into a separate gzipped text file
    following the data storage structure used with other GRACE products
Creates an index file for each data product

CALLING SEQUENCE:
    python cnes_grace_sync.py --release RL04 RL05

OUTPUTS:
    CNES RL01: GAC/GSM
    CNES RL02: GAA/GAB/GSM
    CNES RL03: GAA/GAB/GSM
    CNES RL04: GSM
    CNES RL05: GAA/GAB/GSM

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory: working data directory
    -r X, --release X: CNES/GRGS data releases to sync
    -t X, --timeout X: Timeout in seconds for blocking operations
    -C, --clobber: overwrite existing data in transfer
    -M X, --mode X: Local permissions mode of the directories and files synced
    -l, --log: output log of files downloaded

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 12/2020: use argparse to set command line parameters
        use gravity toolkit utilities to download tar files
        update remote link paths for new structure
    Updated 05/2019: set urllib version based on python version
    Updated 10/2018: added RL04 of the GRGS/CNES product
    Updated 07/2018: python3 compatible method for extracting last modified date
    Updated 06/2018: using python3 compatible octal, input and urllib
    Updated 05/2017: using os.makedirs to recursively create directories
    Updated 04/2017: Sort members from tar file when iterating.  Comments update
        changed from using --rl01 and --rl02 options to setting with --release
        minor changes to check_connection function to parallel other programs
    Updated 02/2017: do not extract tar files to temp, extract contents of files
    Updated 01/2017: v3 of the RL03 dataset: solves the problems identified at
        the poles in RL03-v1 and the problem in C21/S21 identified in RL03-v2
        between January 2003 and December 2012. Added MODE to set permissions
    Updated 10/2016: added --directory option
    Updated 09/2016: v2 of the RL03 dataset (improvements for high-latitudes)
        compress output files with gzip to be similar to other GRACE data
        placed copy file as gzip portion within function (gzip_copy_file)
    Updated 06/2016: added clobber option and date check for individual files
        (will only copy over newer or overwritten files if clobber is not set)
    Updated 05-06/2016: using __future__ print function
    Updated 03/2016: using getopt to set parameters, whether or not to output a
        log file, added new help module.  Updated for latest CNES hosts
    Updated 08/2015: changed sys.exit to raise RuntimeError
    Updated 01/2015: added internet connectivity check
        added main definition for parameters
    Updated 10/2014: updated for RL03
        Individual files no longer available in GRACE format
        for RL01 and RL02.  Using tar files to copy all files to directory
        removing RL01 as option
    Updated 02/2014: quick code update for if statements
    Updated 10/2013: interruptions less frequent but still occur.
        Adding option to create sync logs.  Will output both a file sync log
        and a log of the messages from the wget command
        Updated filepaths with path.join to standardize for different OS
        Deleted commented lftp code (wget is much better at syncing over http)
        May consider switching PO.DAAC sync program to wget to test over ftp
    Updated 09/2013: added subprocess wait commands to prevent interruptions
        in the system call
    Updated 09/2013: switched back to wget as wget was much faster
        at syncing over http.  lftp code is commented and notes remain above.
        non-verbose flag added to output only the files retrieved and errors
    Updated 07/2013: switched to use lftp instead of wget
        Requires less setup as lftp is used with PO.DAAC program
    Updated 05/2013: converted to python
    Updated 03/2013: sync each file with wget versus downloading archive file
        added functionality for RL01 and RL03 (future release)
    Written 07/2012
"""
from __future__ import print_function

import sys
import os
import re
import copy
import time
import gzip
import struct
import shutil
import logging
import tarfile
import argparse
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: sync local GRACE/GRACE-FO files with CNES server
def cnes_grace_sync(DIRECTORY, DREL=[], TIMEOUT=None, LOG=False,
    CLOBBER=False, MODE=None):
    # remote CNES/GRGS host directory
    HOST = ['http://gravitegrace.get.obs-mip.fr','grgs.obs-mip.fr','data']
    # check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    # create dictionaries for dataset and host directory
    DSET = {}
    DSET['RL01'] = ['GSM', 'GAC']
    DSET['RL02'] = ['GSM', 'GAA', 'GAB']
    DSET['RL03'] = ['GSM', 'GAA', 'GAB']
    DSET['RL04'] = ['GSM']
    DSET['RL05'] = ['GSM', 'GAA', 'GAB']

    # remote path to tar files on CNES servers
    REMOTE = dict(RL01={},RL02={},RL03={},RL04={},RL05={})
    # RL01: GSM and GAC
    REMOTE['RL01']['GSM'] = ['RL01','variable','archives']
    REMOTE['RL01']['GAC'] = ['RL01','variable','archives']
    # RL02: GSM, GAA and GAB
    REMOTE['RL02']['GSM'] = ['RL02','variable','archives']
    REMOTE['RL02']['GAA'] = ['RL02','variable','archives']
    REMOTE['RL02']['GAB'] = ['RL02','variable','archives']
    # RL03: GSM, GAA and GAB
    REMOTE['RL03']['GSM'] = ['RL03-v3','archives']
    REMOTE['RL03']['GAA'] = ['RL03','variable','archives']
    REMOTE['RL03']['GAB'] = ['RL03','variable','archives']
    # RL04: GSM
    REMOTE['RL04']['GSM'] = ['RL04-v1','archives']
    # RL05: GSM, GAA, GAB for GRACE/GRACE-FO
    REMOTE['RL05']['GSM'] = ['RL05','archives']
    REMOTE['RL05']['GAA'] = ['RL05','archives']
    REMOTE['RL05']['GAB'] = ['RL05','archives']

    # tar file names for each dataset
    TAR = dict(RL01={},RL02={},RL03={},RL04={},RL05={})
    # RL01: GSM and GAC
    TAR['RL01']['GSM'] = ['GRGS.SH_models.GRACEFORMAT.RL01.tar.gz']
    TAR['RL01']['GAC'] = ['GRGS.dealiasing.RL01.tar.gz']
    # RL02: GSM, GAA and GAB
    TAR['RL02']['GSM'] = ['GRGS.SH_models.GRACEFORMAT.all.tar.gz']
    TAR['RL02']['GAA'] = ['GRGS.dealiasing.GRACEFORMAT.all.tar.gz']
    TAR['RL02']['GAB'] = ['GRGS.dealiasing.GRACEFORMAT.all.tar.gz']
    # RL03: GSM, GAA and GAB
    TAR['RL03']['GSM'] = ['CNES-GRGS.RL03-v3.monthly.coeff.tar.gz']
    TAR['RL03']['GAA'] = ['GRGS.RL03.dealiasing.monthly.tar.gz']
    TAR['RL03']['GAB'] = ['GRGS.RL03.dealiasing.monthly.tar.gz']
    # RL04: GSM
    # TAR['RL04']['GSM'] = ['CNES.RL04-v1.monthly.OLD_IERS2010_MEAN_POLE_CONVENTION.tar.gz']
    TAR['RL04']['GSM'] = ['CNES.RL04-v1.monthly.NEW_IERS2010_MEAN_POLE_CONVENTION.tar.gz']
    # RL05: GSM, GAA and GAB
    TAR['RL05']['GSM'] = ['CNES-GRGS.RL05.GRACE.monthly.tar.gz',
        'CNES-GRGS.RL05.GRACE-FO.monthly.tar.gz']
    TAR['RL05']['GAA'] = ['CNES-GRGS.RL05.monthly.dealiasing.tar.gz']
    TAR['RL05']['GAB'] = ['CNES-GRGS.RL05.monthly.dealiasing.tar.gz']

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # output to log file
        # format: CNES_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'CNES_sync_{today}.log'
        fid1 = open(os.path.join(DIRECTORY,LOGFILE), mode='w', encoding='utf8')
        logging.basicConfig(stream=fid1,level=logging.INFO)
        logging.info(f'CNES Sync Log ({today})')
    else:
        # standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    # DATA RELEASES (RL01, RL02, RL03, RL04)
    # RL01 and RL02 are no longer updated as default
    for rl in DREL:
        # datasets (GSM, GAA, GAB)
        for ds in DSET[rl]:
            logging.info(f'CNES/{rl}/{ds}')
            # specific GRACE directory
            local_dir = os.path.join(DIRECTORY, 'CNES', rl, ds)
            # check if GRACE directory exists and recursively create if not
            os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
            # retrieve each tar file from CNES
            for t in TAR[rl][ds]:
                remote_tar_path = copy.copy(HOST)
                remote_tar_path.extend(REMOTE[rl][ds])
                remote_tar_path.append(t)
                # local copy of CNES data tar file
                local_file = os.path.join(DIRECTORY, 'CNES', rl, t)
                MD5 = gravtk.utilities.get_hash(local_file)
                # copy remote tar file to local if new or updated
                gravtk.utilities.from_http(remote_tar_path,
                    local=local_file, timeout=TIMEOUT, hash=MD5, chunk=16384,
                    verbose=True, fid=fid1, mode=MODE)
                # Create and submit request to get modification time of file
                remote_file = posixpath.join(*remote_tar_path)
                request = gravtk.utilities.urllib2.Request(remote_file)
                response = gravtk.utilities.urllib2.urlopen(request,
                    timeout=TIMEOUT)
                # change modification time to remote
                time_string = response.headers['last-modified']
                remote_mtime = gravtk.utilities.get_unix_time(time_string,
                    format='%a, %d %b %Y %H:%M:%S %Z')
                # keep remote modification time of file and local access time
                os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))

                # open file with tarfile (read)
                tar = tarfile.open(name=local_file, mode='r:gz')

                # copy files from the tar file into the data directory
                member_list=[m for m in tar.getmembers() if re.search(ds,m.name)]
                # for each member of the dataset within the tar file
                for member in member_list:
                    # local gzipped version of the file
                    fi = os.path.basename(member.name)
                    local_file = os.path.join(local_dir, f'{fi}.gz')
                    gzip_copy_file(tar, member, local_file, CLOBBER, MODE)
                # close the tar file
                tar.close()

            # find GRACE files and sort by date
            grace_files = [fi for fi in os.listdir(local_dir) if re.search(ds,fi)]
            # outputting GRACE filenames to index
            index_file = os.path.join(local_dir, 'index.txt')
            with open(index_file, mode='w', encoding='utf8') as fid:
                for fi in sorted(grace_files):
                    print(fi, file=fid)
            # change permissions of index file
            os.chmod(index_file, MODE)

    # close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: copy file from tar file checking if file exists locally
# and if the original file is newer than the local file
def gzip_copy_file(tar, member, local_file, CLOBBER, MODE):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # last modification time of file within tar file
    file1_mtime = member.mtime
    # check if output compressed file exists in local directory
    if os.access(local_file, os.F_OK):
        # check last modification time of output gzipped file
        with gzip.open(local_file, 'rb') as fileID:
            fileobj = fileID.fileobj
            fileobj.seek(4)
            # extract little endian 4 bit unsigned integer
            file2_mtime, = struct.unpack("<I", fileobj.read(4))
        # if remote file is newer: overwrite the local file
        if (file1_mtime > file2_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    # if file does not exist, is to be overwritten, or CLOBBERed
    if TEST or CLOBBER:
        # Printing files copied from tar file to new compressed file
        logging.info(f'{tar.name}/{member.name} --> ')
        logging.info(f'\t{local_file}{OVERWRITE}\n')
        # extract file contents to new compressed file
        f_in = tar.extractfile(member)
        with gzip.GzipFile(local_file, 'wb', 9, None, file1_mtime) as f_out:
            shutil.copyfileobj(f_in, f_out)
        f_in.close()
        # keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, file1_mtime))
        os.chmod(local_file, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""CNES/GRGS GRACE data download program for
            gravity field products
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL05'], choices=['RL01','RL02','RL03','RL04','RL05'],
        help='GRACE/GRACE-FO data release')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # CNES_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
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
    HOST = 'http://gravitegrace.get.obs-mip.fr'
    if gravtk.utilities.check_connection(HOST):
        cnes_grace_sync(args.directory, DREL=args.release,
            TIMEOUT=args.timeout, LOG=args.log,
            CLOBBER=args.clobber, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
