#!/usr/bin/env python
u"""
gfz_isdc_dealiasing_ftp.py
Written by Tyler Sutterley (08/2020)
Syncs GRACE Level-1b dealiasing products from the GFZ Information
    System and Data Center (ISDC)
Optionally outputs as monthly tar files

CALLING SEQUENCE:
    python gfz_isdc_dealiasing_ftp.py --year=2015 --release=RL06 --tar

COMMAND LINE OPTIONS:
    -D X, --directory=X: working data directory
    -R X, --release=X: GRACE data releases to run (RL05,RL06)
    -Y X, --year=X: Years to sync separated by commas
    -T, --tar: Output data as monthly tar files (.tar.gz or .tgz)
    -M X, --mode=X: permissions mode of files synced
    --log: output log of files downloaded
    --clobber: Overwrite existing data in transfers

UPDATE HISTORY:
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
import getopt
import ftplib
import tarfile
import posixpath
import calendar, time

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to the GFZ ftp host
    try:
        f = ftplib.FTP('isdcftp.gfz-potsdam.de')
        f.login()
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: syncs GRACE Level-1b dealiasing products from the GFZ data server
#-- and optionally outputs as monthly tar files
def gfz_isdc_dealiasing_ftp(base_dir, DREL, YEAR=None, TAR=False, LOG=False,
    CLOBBER=False, MODE=None):
    #-- output data directory
    grace_dir = os.path.join(base_dir,'AOD1B',DREL)
    os.makedirs(grace_dir) if not os.access(grace_dir,os.F_OK) else None
    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- output to log file
        #-- format: PODAAC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'GFZ_AOD1B_sync_{0}.log'.format(today)
        fid1 = open(os.path.join(base_dir,LOGFILE),'w')
        print('GFZ AOD1b Sync Log ({0})'.format(today), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- connect and login to GFZ ftp server
    ftp = ftplib.FTP('isdcftp.gfz-potsdam.de')
    ftp.login()

    #-- compile regular expression operator for years to sync
    if YEAR is None:
        regex_years = r'\d{4}'
    else:
        regex_years = r'|'.join(r'{0:d}'.format(y) for y in YEAR)
    #-- compile regular expression operator for years to sync
    R1 = re.compile(r'({0})'.format(regex_years), re.VERBOSE)
    #-- remote subdirectory for DREL on GFZ data server
    remote_sub=posixpath.join(posixpath.sep,'grace','Level-1B','GFZ','AOD',DREL)
    #-- suffix for each data release
    SUFFIX = dict(RL04='tar.gz',RL05='tar.gz',RL06='tgz')

    #-- find remote yearly directories for DREL
    YRS = [R1.search(Y).group(0) for Y in ftp.nlst(remote_sub) if R1.search(Y)]
    for Y in sorted(YRS):
        #-- full path to remote directory
        remote_dir = posixpath.join(remote_sub,Y)
        #-- run for each month
        for M in range(1,13):
            #-- output tar file for year and month
            args = (Y, M, DREL.replace('RL',''), SUFFIX[DREL])
            FILE = 'AOD1B_{0}-{1:02d}_{2}.{3}'.format(*args)
            #-- check if output tar file exists (if TAR)
            TEST = not os.access(os.path.join(grace_dir,FILE), os.F_OK)
            #-- compile regular expressions operators for file dates
            #-- will extract year and month and calendar day from the ascii file
            regex_pattern = r'AOD1B_({0})-({1:02d})-(\d+)_X_\d+.asc.gz$'
            R2 = re.compile(regex_pattern.format(Y,M), re.VERBOSE)
            remote_files = [R2.search(fi).group(0) for fi in
                ftp.nlst(remote_dir) if R2.search(fi)]
            file_count = len(remote_files)
            if TAR and (file_count > 0) and (TEST or CLOBBER):
                #-- copy each gzip file and store within monthly tar files
                tar = tarfile.open(name=os.path.join(grace_dir,FILE),mode='w:gz')
                for fi in sorted(remote_files):
                    #-- remote and local version of each input file
                    remote_file = posixpath.join(remote_dir,fi)
                    local_file = os.path.join(grace_dir,fi)
                    ftp_mirror_file(fid1,ftp,remote_file,local_file,CLOBBER,MODE)
                    #-- add file to tar and then remove local gzip file
                    tar.add(local_file, arcname=fi)
                    os.remove(local_file)
                #-- close tar file and set permissions level to MODE
                tar.close()
                print(' --> {0}\n'.format(os.path.join(grace_dir,FILE),file=fid1))
                os.chmod(os.path.join(grace_dir,FILE), MODE)
            elif (file_count > 0) and not TAR:
                #-- copy each gzip file and keep as individual daily files
                for fi in sorted(remote_files):
                    #-- remote and local version of each input file
                    remote_file = posixpath.join(remote_dir,fi)
                    local_file = os.path.join(grace_dir,fi)
                    ftp_mirror_file(fid1,ftp,remote_file,local_file,CLOBBER,MODE)
    #-- close the ftp connection
    ftp.quit()
    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(base_dir,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def ftp_mirror_file(fid, ftp, remote_file, local_file, CLOBBER, MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- get last modified date of remote file and convert into unix time
    mdtm = ftp.sendcmd('MDTM {0}'.format(remote_file))
    remote_mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
    #-- check if local version of file exists
    if os.access(local_file, os.F_OK):
        #-- check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- Printing files transferred
        print('{0}{1}/{2} --> '.format('ftp://',ftp.host,remote_file),file=fid)
        print('\t{0}{1}\n'.format(local_file,OVERWRITE),file=fid)
        #-- copy remote file contents to local file
        with open(local_file, 'wb') as f:
            ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
        #-- keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking local data directory')
    print(' -R X, --release=X\tGRACE data releases to run (RL05,RL06)')
    print(' -Y X, --year=X\t\tYears to sync separated by commas')
    print(' -T, --tar\t\tOutput data as monthly tar files')
    print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
    print(' -C, --clobber\t\tOverwrite existing data in transfer')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'GFZ_AOD1B_sync_{0}.log'.format(today)
    print('    Log file format: {0}\n'.format(LOGFILE))

#-- Main program that calls gfz_isdc_dealiasing_ftp()
def main():
    #-- Read the system arguments listed after the program
    short_options = 'hD:R:Y:TCM:l'
    long_options = ['help','directory=','release=','year=','tar',
        'log','mode=','clobber']
    optlist,arglist = getopt.getopt(sys.argv[1:],short_options,long_options)

    #-- command line parameters
    base_dir = os.getcwd()
    DREL = 'RL06'
    YEAR = None
    TAR = False
    LOG = False
    #-- permissions mode of the local directories and files
    MODE = 0o775
    CLOBBER = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("--release"):
            DREL = arg
        elif opt in ("-Y","--year"):
            YEAR = [int(Y) for Y in arg.split(',')]
        elif opt in ("-T","--tar"):
            TAR = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-C","--clobber"):
            CLOBBER = True

    #-- check internet connection before attempting to run program
    if check_connection():
        gfz_isdc_dealiasing_ftp(base_dir, DREL=DREL, YEAR=YEAR,
            TAR=TAR, LOG=LOG, MODE=MODE, CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
    main()
