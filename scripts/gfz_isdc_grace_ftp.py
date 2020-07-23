#!/usr/bin/env python
u"""
gfz_isdc_grace_ftp.py
Written by Tyler Sutterley (03/2020)
Syncs GRACE/GRACE-FO data from the GFZ Information System and Data Center (ISDC)
Syncs CSR/GFZ/JPL files for RL06 GAA/GAB/GAC/GAD/GSM
    GAA and GAB are GFZ/JPL only

CALLING SEQUENCE:
    gfz_isdc_grace_ftp(DIRECTORY, PROC)

OUTPUTS:
    CSR RL06: GAC/GAD/GSM
    GFZ RL06: GAA/GAB/GAC/GAD/GSM/AOD1b
    JPL RL06: GAA/GAB/GAC/GAD/GSM

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory=X: working data directory
    --mission: Sync GRACE (grace) or GRACE Follow-On (grace-fo) data
    -C X, --center=X: GRACE/GRACE-FO Processing Center
    -R X, --release=X: GRACE/GRACE-FO data releases to sync (RL05,RL06)
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
    --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode=X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

UPDATE HISTORY:
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
import io
import getopt
import ftplib
import shutil
import hashlib
import posixpath
import calendar, time

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to ftp host for GFZ ISDC
    try:
        f = ftplib.FTP('isdcftp.gfz-potsdam.de')
        f.login()
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: create and compile regular expression operator to find GRACE files
def compile_regex_pattern(PROC, DREL, DSET):
    if ((DSET == 'GSM') and (PROC == 'CSR') and (DREL in ('RL04','RL05'))):
        #-- CSR GSM: only monthly degree 60 products
        #-- not the longterm degree 180, degree 96 dataset or the
        #-- special order 30 datasets for the high-resonance months
        release, = re.findall('\d+', DREL)
        args = (DSET, int(release))
        regex_pattern='{0}-2_\d+-\d+_\d+_UTCSR_0060_000{1:d}.gz$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'CSR') and (DREL == 'RL06')):
        #-- CSR GSM RL06: only monthly degree 60 products
        release, = re.findall('\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern='{0}-2_\d+-\d+_{1}_UTCSR_{2}_0{3:d}00.gz$' .format(*args)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL04')):
        #-- GFZ RL04: only unconstrained solutions (not GK2 products)
        regex_pattern='{0}-2_\d+-\d+_\d+_EIGEN_G---_0004.gz$'.format(DSET)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL05')):
        #-- GFZ RL05: updated RL05a products which are less constrained to
        #-- the background model.  Allow regularized fields
        regex_unconst='{0}-2_\d+-\d+_\d+_EIGEN_G---_005a.gz$'.format(DSET)
        regex_regular='{0}-2_\d+-\d+_\d+_EIGEN_GK2-_005a.gz$'.format(DSET)
        regex_pattern='{0}|{1}'.format(regex_unconst,regex_regular)
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL06')):
        #-- GFZ GSM RL06: only monthly degree 60 products
        release, = re.findall('\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern='{0}-2_\d+-\d+_{1}_GFZOP_{2}_0{3:d}00.gz$' .format(*args)
    elif (PROC == 'JPL') and DREL in ('RL04','RL05'):
        #-- JPL: RL04a and RL05a products (denoted by 0001)
        release, = re.findall('\d+', DREL)
        args = (DSET, int(release))
        regex_pattern='{0}-2_\d+-\d+_\d+_JPLEM_0001_000{1:d}.gz$'.format(*args)
    elif ((DSET == 'GSM') and (PROC == 'JPL') and (DREL == 'RL06')):
        #-- JPL GSM RL06: only monthly degree 60 products
        release, = re.findall('\d+', DREL)
        args = (DSET, '(GRAC|GRFO)', 'BA01', int(release))
        regex_pattern='{0}-2_\d+-\d+_{1}_JPLEM_{2}_0{3:d}00.gz$' .format(*args)
    else:
        regex_pattern='{0}-2_(.*?).gz$'.format(DSET)
    #-- return the compiled regular expression operator used to find files
    return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local GRACE/GRACE-FO files with GFZ ISDC server
def gfz_isdc_grace_ftp(DIRECTORY, PROC, DREL=[], MISSION=['grace','grace-fo'],
    LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=None):

    #-- connect and login to GFZ ISDC ftp server
    ftp = ftplib.FTP('isdcftp.gfz-potsdam.de')
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
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('GFZ ISDC Sync Log ({0})'.format(today), file=fid1)
        print('CENTERS={0}'.format(','.join(PROC)), file=fid1)
        print('RELEASES={0}'.format(','.join(DREL)), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- GRACE/GRACE-FO DATA
    for mi in MISSION:
        #-- PROCESSING CENTERS (CSR, GFZ, JPL)
        print('{0} L2 Global Spherical Harmonics:'.format(mi.upper()),file=fid1)
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
                #-- remote directory for data release
                remote_dir = posixpath.join(mi, 'Level-2', pr, drel_str)
                #-- DATA PRODUCTS (GAC GAD GSM GAA GAB)
                for ds in DSET[pr]:
                    #-- print string of exact data product
                    print('{0}/{1}/{2}'.format(pr, drel_str, ds), file=fid1)
                    #-- local directory for exact data product
                    local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                    #-- check if directory exists and recursively create if not
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir,MODE)
                    #-- compile the regular expression operator to find files
                    R1 = re.compile('({0}-(.*?)(gz|txt|dif))'.format(ds))
                    #-- get filenames from remote directory
                    lines = [fi for fi in ftp.nlst(remote_dir) if R1.search(fi)]
                    for line in sorted(lines):
                        #-- extract filename from regex object
                        fi = R1.search(line).group(0)
                        remote_file = posixpath.join(remote_dir,fi)
                        local_file = os.path.join(local_dir,fi)
                        ftp_mirror_file(fid1, ftp, remote_file, local_file,
                            LIST, CLOBBER, CHECKSUM, MODE)

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
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def ftp_mirror_file(fid,ftp,remote_file,local_file,LIST,CLOBBER,CHECKSUM,MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- get last modified date of remote file and convert into unix time
    mdtm = ftp.sendcmd('MDTM {0}'.format(remote_file))
    remote_mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
    #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(local_file, 'rb') as local_buffer:
            local_hash = hashlib.md5(local_buffer.read()).hexdigest()
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary('RETR {0}'.format(remote_file), remote_buffer.write)
        remote_buffer.seek(0)
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
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- Printing files transferred
        arg=(posixpath.join('ftp://',ftp.host,remote_file),local_file,OVERWRITE)
        fid.write('{0} -->\n\t{1}{2}\n\n'.format(*arg))
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- copy file from ftp server or from bytesIO object
            if CHECKSUM and os.access(local_file, os.F_OK):
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, 16 * 1024)
            else:
                #-- copy remote file contents to local file
                with open(local_file, 'wb') as f:
                    ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking Data Directory')
    print(' --mission\t\tSync GRACE (grace) or GRACE Follow-On (grace-fo) data')
    print(' -C X, --center=X\tGRACE Processing Center (CSR,GFZ,JPL)')
    print(' -R X, --release=X\tGRACE data releases to sync (RL04,RL05)')
    print(' -L, --list\t\tOnly print files that are to be transferred')
    print(' --clobber\t\tOverwrite existing data in transfer')
    print(' --checksum\t\tCompare hashes to check if overwriting existing data')
    print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'GFZ_ISDC_sync_{0}.log'.format(today)
    print('    Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls gfz_isdc_grace_ftp()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','list','log','mission=','center=',
        'release=','clobber','checksum','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:lLC:R:M:',long_options)

    #-- command line parameters
    DIRECTORY = os.getcwd()
    LIST = False
    LOG = False
    CLOBBER = False
    #-- mission (GRACE or GRACE Follow-On)
    MISSION = ['grace','grace-fo']
    #-- GRACE/GRACE-FO Processing Centers to run
    PROC = ['CSR', 'GFZ', 'JPL']
    #-- Data release
    DREL = ['RL06']
    #-- Use hash for determining whether or not to overwrite
    CHECKSUM = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("-L","--list"):
            LIST = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("--clobber"):
            CLOBBER = True
        elif opt in ("--mission"):
            MISSION = arg.lower().split(',')
        elif opt in ("-C","--center"):
            PROC = arg.upper().split(',')
        elif opt in ("-R","--release"):
            DREL = arg.upper().split(',')
        elif opt in ("--checksum"):
            CHECKSUM = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- check internet connection before attempting to run program
    if check_connection():
        gfz_isdc_grace_ftp(DIRECTORY, PROC, DREL=DREL, MISSION=MISSION,
            LIST=LIST, LOG=LOG, CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)
    else:
        raise RuntimeError('Check internet connection')

#-- run main program
if __name__ == '__main__':
    main()
