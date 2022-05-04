#!/usr/bin/env python
u"""
gfz_isdc_grace_ftp.py
Written by Tyler Sutterley (04/2022)
Syncs GRACE/GRACE-FO data from the GFZ Information System and Data Center (ISDC)
Syncs CSR/GFZ/JPL files for RL06 GAA/GAB/GAC/GAD/GSM
    GAA and GAB are GFZ/JPL only
Gets the latest technical note (TN) files
Gets the monthly GRACE/GRACE-FO newsletters

CALLING SEQUENCE:
    python gfz_isdc_grace_ftp.py

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
import copy
import time
import ftplib
import shutil
import hashlib
import logging
import argparse
import posixpath
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

#-- PURPOSE: sync local GRACE/GRACE-FO files with GFZ ISDC server
def gfz_isdc_grace_ftp(DIRECTORY, PROC=[], DREL=[], VERSION=[],
    NEWSLETTERS=False, TIMEOUT=None, LOG=False, LIST=False,
    CLOBBER=False, CHECKSUM=False, MODE=None):

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
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('GFZ ISDC Sync Log ({0})'.format(today))
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- Degree 1 (geocenter) coefficients
    logging.info('Degree 1 Coefficients:')
    local_dir = os.path.join(DIRECTORY,'geocenter')
    #-- check if geocenter directory exists and recursively create if not
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    #-- TN-13 JPL degree 1 files
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-13_GEOC_(CSR|GFZ|JPL)_(.*?).txt$', re.VERBOSE)
    #-- get filenames from remote directory
    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
        [ftp.host,'grace-fo','DOCUMENTS','TECHNICAL_NOTES'], timeout=TIMEOUT,
        basename=True, pattern=R1, sort=True)
    #-- for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        #-- extract filename from regex object
        remote_path = [ftp.host,'grace-fo','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = os.path.join(local_dir,fi)
        ftp_mirror_file(ftp, remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

    #-- SLR C2,0 coefficients
    logging.info('C2,0 Coefficients:')
    local_dir = os.path.expanduser(DIRECTORY)
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-(05|07|11)_C20_SLR_RL(.*?).txt$', re.VERBOSE)
    #-- get filenames from remote directory
    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
        [ftp.host,'grace','DOCUMENTS','TECHNICAL_NOTES'], timeout=TIMEOUT,
        basename=True, pattern=R1, sort=True)
    #-- for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        #-- extract filename from regex object
        remote_path = [ftp.host,'grace','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = os.path.join(local_dir,re.sub(r'(_RL.*?).txt','.txt',fi))
        ftp_mirror_file(ftp, remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

    #-- SLR C3,0 coefficients
    logging.info('C3,0 Coefficients:')
    local_dir = os.path.expanduser(DIRECTORY)
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-(14)_C30_C20_SLR_GSFC.txt$', re.VERBOSE)
    #-- get filenames from remote directory
    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
        [ftp.host,'grace-fo','DOCUMENTS','TECHNICAL_NOTES'], timeout=TIMEOUT,
        basename=True, pattern=R1, sort=True)
    #-- for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        #-- extract filename from regex object
        remote_path = [ftp.host,'grace-fo','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = os.path.join(local_dir,re.sub(r'(SLR_GSFC)','GSFC_SLR',fi))
        ftp_mirror_file(ftp, remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

    #-- TN-08 GAE, TN-09 GAF and TN-10 GAG ECMWF atmosphere correction products
    logging.info('TN-08 GAE, TN-09 GAF and TN-10 GAG products:')
    local_dir = os.path.expanduser(DIRECTORY)
    ECMWF_files = []
    ECMWF_files.append('TN-08_GAE-2_2006032-2010031_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-09_GAF-2_2010032-2015131_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-10_GAG-2_2015132-2099001_0000_EIGEN_G---_0005.gz')
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'({0}|{1}|{2})'.format(*ECMWF_files), re.VERBOSE)
    #-- get filenames from remote directory
    remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
        [ftp.host,'grace','DOCUMENTS','TECHNICAL_NOTES'], timeout=TIMEOUT,
        basename=True, pattern=R1, sort=True)
    #-- for each file on the remote server
    for fi,remote_mtime in zip(remote_files,remote_mtimes):
        #-- extract filename from regex object
        remote_path = [ftp.host,'grace','DOCUMENTS','TECHNICAL_NOTES',fi]
        local_file = os.path.join(local_dir,fi)
        ftp_mirror_file(ftp, remote_path, remote_mtime,
            local_file, TIMEOUT=TIMEOUT, LIST=LIST,
            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

    #-- GRACE and GRACE-FO newsletters
    if NEWSLETTERS:
        #-- local newsletter directory (place GRACE and GRACE-FO together)
        local_dir = os.path.join(DIRECTORY,'newsletters')
        #-- check if newsletters directory exists and recursively create if not
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- for each satellite mission (grace, grace-fo)
        for i,mi in enumerate(['grace','grace-fo']):
            logging.info('{0} Newsletters:'.format(mi))
            #-- compile regular expression operator for remote files
            NAME = mi.upper().replace('-','_')
            R1 = re.compile(r'{0}_SDS_NL_(\d+).pdf'.format(NAME), re.VERBOSE)
            #-- find years for GRACE/GRACE-FO newsletters
            years,_  = gravity_toolkit.utilities.ftp_list(
                [ftp.host,mi,'DOCUMENTS','NEWSLETTER'], timeout=TIMEOUT,
                basename=True, pattern=r'\d+', sort=True)
            #-- for each year of GRACE/GRACE-FO newsletters
            for Y in years:
                #-- find GRACE/GRACE-FO newsletters
                remote_files,remote_mtimes = gravity_toolkit.utilities.ftp_list(
                    [ftp.host,mi,'DOCUMENTS','NEWSLETTER',Y], timeout=TIMEOUT,
                    basename=True, pattern=R1, sort=True)
                #-- for each file on the remote server
                for fi,remote_mtime in zip(remote_files,remote_mtimes):
                    #-- extract filename from regex object
                    remote_path = [ftp.host,mi,'DOCUMENTS','NEWSLETTER',Y,fi]
                    local_file = os.path.join(local_dir,fi)
                    ftp_mirror_file(ftp, remote_path, remote_mtime,
                        local_file, TIMEOUT=TIMEOUT, LIST=LIST,
                        CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

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
                    #-- modifiers for intermediate data releases
                    if (int(VERSION) > 0):
                        drel_str = '{0}.{1}'.format(rl,VERSION[i])
                    else:
                        drel_str = copy.copy(rl)
                    #-- print string of exact data product
                    logging.info('{0}/{1}/{2}/{3}'.format(mi, pr, drel_str, ds))
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
                    #-- regular expression operator for data product
                    rx = compile_regex_pattern(pr, rl, ds, version=VERSION[i])
                    #-- find local GRACE/GRACE-FO files to create index
                    files = [fi for fi in os.listdir(local_dir) if rx.match(fi)]
                    #-- extend list of GRACE/GRACE-FO files
                    grace_files.extend(files)

                #-- outputting GRACE/GRACE-FO filenames to index
                with open(os.path.join(local_dir,'index.txt'),'w') as fid:
                    for fi in sorted(grace_files):
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

#-- PURPOSE: create argument parser
def arguments():
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
    #-- GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs='+',
        default=['0','1'], choices=['0','1','2','3'],
        help='GRACE/GRACE-FO Level-2 data version')
    #-- GRACE/GRACE-FO newsletters
    parser.add_argument('--newsletters','-n',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Newsletters')
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
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- check internet connection before attempting to run program
    HOST = 'isdcftp.gfz-potsdam.de'
    if gravity_toolkit.utilities.check_ftp_connection(HOST):
        gfz_isdc_grace_ftp(args.directory, PROC=args.center,
            DREL=args.release, VERSION=args.version,
            NEWSLETTERS=args.newsletters, TIMEOUT=args.timeout,
            LIST=args.list, LOG=args.log, CLOBBER=args.clobber,
            CHECKSUM=args.checksum, MODE=args.mode)
    else:
        raise RuntimeError('Check internet connection')

#-- run main program
if __name__ == '__main__':
    main()
