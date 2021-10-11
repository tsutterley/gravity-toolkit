#!/usr/bin/env python
u"""
gfz_isdc_grace_ftp.py
Written by Tyler Sutterley (10/2021)
Syncs GRACE/GRACE-FO data from the GFZ Information System and Data Center (ISDC)
Syncs CSR/GFZ/JPL files for RL06 GAA/GAB/GAC/GAD/GSM
    GAA and GAB are GFZ/JPL only

CALLING SEQUENCE:
    python gfz_isdc_grace_ftp.py --mission grace grace-fo

OUTPUTS:
    CSR RL06: GAC/GAD/GSM
    GFZ RL06: GAA/GAB/GAC/GAD/GSM/AOD1b
    JPL RL06: GAA/GAB/GAC/GAD/GSM

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -m X, --mission X: Sync GRACE (grace) or GRACE Follow-On (grace-fo) data
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO data releases to sync (RL05,RL06)
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
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 01/2021: using utilities module to list files from ftp
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 03/2020 for public release
    Updated 03/2020: new GFZ ISDC ftp server website
    Updated 09/2019: added checksum option to not overwrite existing data files
        added GRACE Follow-On data sync
    Written 08/2018
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
import gravity_toolkit.utilities

#-- PURPOSE: create and compile regular expression operator to find GRACE files
def compile_regex_pattern(PROC, DREL, DSET):
    if ((DSET == 'GSM') and (PROC == 'CSR') and (DREL in ('RL04','RL05'))):
        #-- CSR GSM: only monthly degree 60 products
        #-- not the longterm degree 180, degree 96 dataset or the
        #-- special order 30 datasets for the high-resonance months
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        regex_pattern=r'{0}-2_\d+-\d+_\d+_UTCSR_0060_000{1:d}.gz$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'CSR') and (DREL == 'RL06')):
        #-- CSR GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_UTCSR_{2}_0{3:d}00.gz$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL04')):
        #-- GFZ RL04: only unconstrained solutions (not GK2 products)
        regex_pattern=r'{0}-2_\d+-\d+_\d+_EIGEN_G---_0004.gz$'.format(DSET)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL05')):
        #-- GFZ RL05: updated RL05a products which are less constrained to
        #-- the background model.  Allow regularized fields
        regex_unconst=r'{0}-2_\d+-\d+_\d+_EIGEN_G---_005a.gz$'.format(DSET)
        regex_regular=r'{0}-2_\d+-\d+_\d+_EIGEN_GK2-_005a.gz$'.format(DSET)
        regex_pattern=r'{0}|{1}'.format(regex_unconst,regex_regular)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL06')):
        #-- GFZ GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_GFZOP_{2}_0{3:d}00.gz$' .format(*args)
    elif (PROC == 'JPL') and DREL in ('RL04','RL05'):
        #-- JPL: RL04a and RL05a products (denoted by 0001)
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        regex_pattern=r'{0}-2_\d+-\d+_\d+_JPLEM_0001_000{1:d}.gz$'.format(*args)
    elif ((DSET == 'GSM') and (PROC == 'JPL') and (DREL == 'RL06')):
        #-- JPL GSM RL06: only monthly degree 60 products
        release, = re.findall(r'\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern=r'{0}-2_\d+-\d+_{1}_JPLEM_{2}_0{3:d}00.gz$' .format(*args)
    else:
        regex_pattern=r'{0}-2_(.*?).gz$'.format(DSET)
    #-- return the compiled regular expression operator used to find files
    return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local GRACE/GRACE-FO files with GFZ ISDC server
def gfz_isdc_grace_ftp(DIRECTORY, PROC, DREL=[], MISSION=[], TIMEOUT=None,
    LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=None):

    #-- connect and login to GFZ ISDC ftp server
    ftp = ftplib.FTP('isdcftp.gfz-potsdam.de', timeout=TIMEOUT)
    ftp.login()

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    #-- datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- output to log file
        #-- format: GFZ_ISDC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'GFZ_ISDC_sync_{0}.log'.format(today)
        logging.basicConfig(file=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('GFZ ISDC Sync Log ({0})'.format(today))
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- GRACE/GRACE-FO DATA
    for mi in MISSION:
        #-- PROCESSING CENTERS (CSR, GFZ, JPL)
        logging.info('{0} L2 Global Spherical Harmonics:'.format(mi.upper()))
        for pr in PROC:
            #-- DATA RELEASES (RL04, RL05, RL06)
            for rl in DREL:
                #-- modifiers for intermediate data releases
                if (pr == 'JPL') and (rl == 'RL04'):
                    #-- JPL RELEASE 4 = RL4.1
                    drel_str = '{0}.1'.format(rl)
                elif (pr == 'JPL') and (rl == 'RL05'):
                    #-- JPL RELEASE 5 = RL05.1 (11/2014)
                    drel_str = '{0}.1'.format(rl)
                else:
                    drel_str = rl
                #-- DATA PRODUCTS (GAC GAD GSM GAA GAB)
                for ds in DSET[pr]:
                    #-- print string of exact data product
                    logging.info('{0}/{1}/{2}'.format(pr, drel_str, ds))
                    #-- local directory for exact data product
                    local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                    #-- check if directory exists and recursively create if not
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir,MODE)
                    #-- compile the regular expression operator to find files
                    R1 = re.compile(r'({0}-(.*?)(gz|txt|dif))'.format(ds))
                    #-- get filenames from remote directory
                    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
                        [ftp.host,mi,'Level-2',pr,drel_str], timeout=TIMEOUT,
                        basename=True, pattern=R1, sort=True)
                    for fi,remote_mtime in zip(remote_files,remote_mtimes):
                        #-- extract filename from regex object
                        remote_path = [ftp.host,mi,'Level-2',pr,drel_str,fi]
                        local_file = os.path.join(local_dir,fi)
                        ftp_mirror_file(ftp, remote_path, remote_mtime,
                            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
                            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

                    #-- Create an index file for each GRACE/GRACE-FO product
                    #-- finding all dataset files *.gz in directory
                    rx = compile_regex_pattern(pr, rl, ds)
                    #-- find local GRACE/GRACE-FO files to create index
                    files = [fi for fi in os.listdir(local_dir) if rx.match(fi)]
                    #-- outputting GRACE/GRACE-FO filenames to index
                    with open(os.path.join(local_dir,'index.txt'),'w') as fid:
                        for fi in sorted(files):
                            print('{0}'.format(fi), file=fid)
                    #-- change permissions of index file
                    os.chmod(os.path.join(local_dir,'index.txt'), MODE)

    #-- close the ftp connection
    ftp.quit()
    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def ftp_mirror_file(ftp,remote_path,remote_mtime,local_file,
    TIMEOUT=None,LIST=False,CLOBBER=False,CHECKSUM=False,MODE=0o775):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(local_file, 'rb') as local_buffer:
            local_hash = hashlib.md5(local_buffer.read()).hexdigest()
        #-- copy remote file contents to bytesIO object
        remote_buffer = gravity_toolkit.utilities.from_ftp(remote_path,
            timeout=TIMEOUT)
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = ' (checksums: {0} {1})'.format(local_hash,remote_hash)
    elif os.access(local_file, os.F_OK):
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
        arg=(posixpath.join('ftp://',*remote_path),local_file,OVERWRITE)
        logging.info('{0} -->\n\t{1}{2}\n'.format(*arg))
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- copy file from ftp server or from bytesIO object
            if CHECKSUM and os.access(local_file, os.F_OK):
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, 16 * 1024)
            else:
                #-- path to remote file
                remote_file = posixpath.join(*remote_path[1:])
                #-- copy remote file contents to local file
                with open(local_file, 'wb') as f:
                    ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- Main program that calls gfz_isdc_grace_ftp()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO data from the GFZ
            Information System and Data Center (ISDC)
            """
    )
    #-- command line parameters
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- mission (GRACE or GRACE Follow-On)
    parser.add_argument('--mission','-m',
        type=str, nargs='+',
        default=['grace','grace-fo'], choices=['grace','grace-fo'],
        help='Mission to sync between GRACE and GRACE-FO')
    #-- GRACE/GRACE-FO processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+',
        default=['CSR','GFZ','JPL'], choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO processing center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, nargs='+',
        default=['RL06'], choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    #-- connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    #-- Output log file in form
    #-- GFZ_ISDC_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args,_ = parser.parse_known_args()

    #-- check internet connection before attempting to run program
    HOST = 'isdcftp.gfz-potsdam.de'
    if gravity_toolkit.utilities.check_ftp_connection(HOST):
        gfz_isdc_grace_ftp(args.directory, args.center, DREL=args.release,
            MISSION=args.mission, TIMEOUT=args.timeout, LIST=args.list,
            LOG=args.log, CLOBBER=args.clobber, CHECKSUM=args.checksum,
            MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

#-- run main program
if __name__ == '__main__':
    main()
