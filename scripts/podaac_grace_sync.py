#!/usr/bin/env python
u"""
podaac_grace_sync.py
Written by Tyler Sutterley (04/2022)

Syncs GRACE/GRACE-FO and auxiliary data from the NASA JPL PO.DAAC Drive Server
Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM
    GAA and GAB are GFZ/JPL only
Syncs GFZ AOD1b files for RL04/RL05/RL06
Gets the latest technical note (TN) files
Gets the monthly GRACE/GRACE-FO newsletters

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add PO.DAAC Drive OPS to NASA Earthdata Applications and get WebDAV Password
https://podaac-tools.jpl.nasa.gov/drive

CALLING SEQUENCE:
    python podaac_grace_sync.py --user <username>
    where <username> is your NASA Earthdata username

OUTPUTS:
    CSR RL04/RL05/RL06: GAC/GAD/GSM
    GFZ RL04/RL05/RL06: GAA/GAB/GAC/GAD/GSM
    JPL RL04/RL05/RL06: GAA/GAB/GAC/GAD/GSM
    Tellus degree one coefficients (TN-13)
    Technical notes for satellite laser ranging coefficients
    Technical notes for Release-05 atmospheric corrections
    GFZ RL04/RL05/RL06: Level-1b dealiasing solutions
    Monthly GRACE/GRACE-FO newsletters

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --webdav X: WebDAV password for JPL PO.DAAC Drive Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases to sync
    -v X, --version X: GRACE/GRACE-FO Level-2 Data Version to sync
    -a, --aod1b: sync GRACE/GRACE-FO Level-1B dealiasing products
    -n, --newsletters: sync GRACE/GRACE-FO newsletters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
    -C, --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 04/2022: added option for GRACE/GRACE-FO Level-2 data version
        refactor to always try syncing from both grace and grace-fo missions
        use granule identifiers from CMR query to build output file index
        use argparse descriptions within sphinx documentation
    Updated 03/2022: update regular expression pattern for finding files
        use CMR queries for finding GRACE/GRACE-FO level-2 product urls
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 12/2020: generalized podaac_list() by renaming to drive_list()
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
        moved urllib opener to utilities. add credential check
        moved urllib directory listing to utilities
    Updated 07/2020: add back snippets to sync Level-1b dealiasing products
    Updated 06/2020: increased timeout to 2 minutes
    Updated 05/2020: simplified PO.DAAC Drive login
        added netrc option for alternative authentication method
    Updated 03/2020 for public release.  Set default release to RL06
    Updated 12/2019: GSFC TN-14 oblateness and C30 file in gracefo documents
        convert last modified time in a function for a given format
    Updated 10/2019: add GRACE-FO newsletters
    Updated 09/2019: added ssl context to urlopen headers
        added checksum option to not overwrite existing data files
    Updated 07/2019: added GSFC TN-14 oblateness and C30 files
    Updated 06/2019: added JPL TN-13 geocenter files (CSR/GFZ/JPL)
        added GRACE-FO spherical harmonic coefficient files.  RL06 AOD1b *.tgz
    Updated 04/2019: new podaac drive website https://podaac-tools.jpl.nasa.gov
    Updated 12/2018: decode authorization header for python3 compatibility
    Updated 11/2018: encode base64 strings for python3 compatibility
    Updated 08/2018: regular expression pattern for GFZ Release 6
        new lxml expressions to match PO.DAAC Drive website updates
    Updated 06/2018: using python3 compatible octal, input and urllib
        updated regular expression pattern for JPL Release 6
    Updated 05/2018: updated for CSR Release 6
    Updated 03/2018: updated header for log file output
    Updated 11/2017: increased urllib2.urlopen timeout to 20.  PROC as option
    Updated 08/2017: use raw_input() to enter NASA Earthdata credentials rather
        than exiting with error
    Updated 05/2017: exception if NASA Earthdata credentials weren't entered
        using os.makedirs to recursively create directories
        using getpass to enter server password securely (remove --password)
    Updated 04/2017: slight modification to AOD1b regular expression
        changed from using --rl04 option to setting with --release
        parsing HTML with lxml libraries (extract names and last modified dates)
        Authorization header in urllib2 opener. changes to check_connection()
    Updated 01-02/2017: converted to https (urllib2) for NASA Earthdata servers
        GRACE newsletters no longer standard, use --newsletters option
    Updated 01/2017: added --mode to set file and directory permissions. Cleanup
    Updated 09/2016: rewritten to use libftp within python (no external calls)
        previous lftp version renamed podaac_grace_lftp.py
    Updated 06/2016: added --list option for a dry-run (do not transfer files)
        absolute import of shutil package
    Updated 05/2016: using __future__ print function
    Updated 03/2016: using getopt to set parameters, whether or not to output a
        log file, added new help module
    Updated 02/2016: added option for CSR weekly 5X5 harmonics
    Updated 01/2016: parallel downloading for lftp mirror portions (8 files)
    Updated 12/2015: added ECMWF GAG "jump" correction. New TN-09 GAF correction
        using clobber flag for mget portions (rather than initial delete)
    Updated 08/2015: changed sys.exit to raise RuntimeError
        Sync the TN-08 and TN-09 GAE and GAF products
    Updated 05/2015: updated for Jan/Feb 2015 months (don't include LMAXx30)
        improved regular expression usage to create indices
    Updated 03/2015: update for JPL RL05.1 (see L2-JPL_ProcStds_v5.1.pdf)
    Updated 01/2015: added internet connectivity check
    Updated 11/2014: added main definition for parameters
    Updated 09/2014: updates to code structure.  Removed globs
    Updated 02/2014: minor update to if statements
    Updated 12/2013: Updated for GFZ RL05a GSM products
        These products are less constrained to the background model
        and are denoted with *_005a.gz
        Final directory is starting working directory (using getcwd)
    Updated 10/2013: Adding option to create sync logs.
        Updated filepaths with path.join to standardize for different OS
    Updated 09/2013: added subprocess wait commands to prevent interruptions
        in the system call
    Updated 05/2013: added RETIRE flag for RL04 as podaac directory moved
        added sorted to glob for linux computers
        glob parallels ls -U and has to be sorted
    Updated 03/2013: added argument tags for AOD1B, RL04, ECCO, and GLDAS
    Updated 01/2013: downloads daily AOD1B products (GFZ RL04 and RL05)
    Updated 11/2012: Turned off loop for RL04 as it is no longer updated and
    Updated 07/2012: adjusted the output directory for generalization
    Updated 06/2012: added in new SLR dataset for release 5 TN-07_C20_SLR.txt
    Updated 06/2012: index file created only lists the wanted *.gz files
    Updated 04/2012: added loops to run through the different products
    Written 04/2012
"""
from __future__ import print_function

import sys
import os
import re
import io
import time
import netrc
import shutil
import getpass
import logging
import argparse
import builtins
import posixpath
import lxml.etree
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

#-- PURPOSE: sync local GRACE/GRACE-FO files with JPL PO.DAAC drive server
def podaac_grace_sync(DIRECTORY, PROC=[], DREL=[], VERSION=[],
    AOD1B=False, NEWSLETTERS=False, TIMEOUT=None, LOG=False, LIST=False,
    CLOBBER=False, CHECKSUM=False, MODE=None):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- remote https server for GRACE data
    HOST = 'https://podaac-tools.jpl.nasa.gov'
    #-- RL04/RL05 have been moved on PO.DAAC to the retired directory
    retired = {}
    retired['RL04'] = 'retired'
    retired['RL05'] = 'retired'
    retired['RL06'] = ''
    #-- datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    #-- remote subdirectories for newsletters (note capital for grace-fo)
    newsletter_sub = {}
    newsletter_sub['grace'] = ['grace','docs','newsletters']
    newsletter_sub['grace-fo'] = ['gracefo','docs','Newsletters']
    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: PODAAC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'PODAAC_sync_{0}.log'.format(today)
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('PO.DAAC Sync Log ({0})'.format(today))
        logging.info('CENTERS={0}'.format(','.join(PROC)))
        logging.info('RELEASES={0}'.format(','.join(DREL)))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- Degree 1 (geocenter) coefficients
    logging.info('Degree 1 Coefficients:')
    PATH = [HOST,'drive','files','allData','tellus','L2','degree_1']
    remote_dir = posixpath.join(*PATH)
    local_dir = os.path.join(DIRECTORY,'geocenter')
    #-- check if geocenter directory exists and recursively create if not
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    #-- TN-13 JPL degree 1 files
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-13_GEOC_(CSR|GFZ|JPL)_(.*?).txt', re.VERBOSE)
    #-- open connection with PO.DAAC drive server at remote directory
    files,mtimes = gravity_toolkit.utilities.drive_list(PATH,
        timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
    #-- for each file on the remote server
    for colname,remote_mtime in zip(files,mtimes):
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colname)
        local_file = os.path.join(local_dir,colname)
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)

    #-- SLR C2,0 coefficients
    logging.info('C2,0 Coefficients:')
    PATH = [HOST,'drive','files','allData','grace','docs']
    remote_dir = posixpath.join(*PATH)
    local_dir = os.path.expanduser(DIRECTORY)
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-(05|07|11)_C20_SLR.txt', re.VERBOSE)
    #-- open connection with PO.DAAC drive server at remote directory
    files,mtimes = gravity_toolkit.utilities.drive_list(PATH,
        timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
    #-- for each file on the remote server
    for colname,remote_mtime in zip(files,mtimes):
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colname)
        local_file = os.path.join(local_dir,colname)
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)

    #-- SLR C3,0 coefficients
    logging.info('C3,0 Coefficients:')
    PATH = [HOST,'drive','files','allData','gracefo','docs']
    remote_dir = posixpath.join(*PATH)
    local_dir = os.path.expanduser(DIRECTORY)
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'TN-(14)_C30_C20_GSFC_SLR.txt', re.VERBOSE)
    #-- open connection with PO.DAAC drive server at remote directory
    files,mtimes = gravity_toolkit.utilities.drive_list(PATH,
        timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
    #-- for each file on the remote server
    for colname,remote_mtime in zip(files,mtimes):
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colname)
        local_file = os.path.join(local_dir,colname)
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)

    #-- TN-08 GAE, TN-09 GAF and TN-10 GAG ECMWF atmosphere correction products
    logging.info('TN-08 GAE, TN-09 GAF and TN-10 GAG products:')
    PATH = [HOST,'drive','files','allData','grace','docs']
    remote_dir = posixpath.join(*PATH)
    local_dir = os.path.expanduser(DIRECTORY)
    ECMWF_files = []
    ECMWF_files.append('TN-08_GAE-2_2006032-2010031_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-09_GAF-2_2010032-2015131_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-10_GAG-2_2015132-2099001_0000_EIGEN_G---_0005.gz')
    #-- compile regular expression operator for remote files
    R1 = re.compile(r'({0}|{1}|{2})'.format(*ECMWF_files), re.VERBOSE)
    #-- open connection with PO.DAAC drive server at remote directory
    files,mtimes = gravity_toolkit.utilities.drive_list(PATH,
        timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
    #-- for each file on the remote server
    for colname,remote_mtime in zip(files,mtimes):
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colname)
        local_file = os.path.join(local_dir,colname)
        http_pull_file(remote_file, remote_mtime, local_file,
            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
            CHECKSUM=CHECKSUM, MODE=MODE)

    #-- GRACE and GRACE-FO newsletters
    if NEWSLETTERS:
        #-- local newsletter directory (place GRACE and GRACE-FO together)
        local_dir = os.path.join(DIRECTORY,'newsletters')
        #-- check if newsletters directory exists and recursively create if not
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- for each satellite mission (grace, grace-fo)
        for i,mi in enumerate(['grace','grace-fo']):
            logging.info('{0} Newsletters:'.format(mi))
            PATH = [HOST,'drive','files','allData',*newsletter_sub[mi]]
            remote_dir = posixpath.join(*PATH)
            #-- compile regular expression operator for remote files
            NAME = mi.upper().replace('-','_')
            R1 = re.compile(r'{0}_SDS_NL_(\d+).pdf'.format(NAME), re.VERBOSE)
            #-- open connection with PO.DAAC drive server at remote directory
            files,mtimes = gravity_toolkit.utilities.drive_list(PATH,
                timeout=TIMEOUT,build=False,parser=parser,pattern=R1,sort=True)
            #-- for each file on the remote server
            for colname,remote_mtime in zip(files,mtimes):
                #-- remote and local versions of the file
                remote_file = posixpath.join(remote_dir,colname)
                local_file = os.path.join(local_dir,colname)
                http_pull_file(remote_file, remote_mtime, local_file,
                    TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                    CHECKSUM=CHECKSUM, MODE=MODE)

    #-- GRACE/GRACE-FO AOD1B dealiasing products
    if AOD1B:
        logging.info('GRACE L1B Dealiasing Products:')
        #-- for each data release (RL04, RL05, RL06)
        for rl in DREL:
            #-- print string of exact data product
            logging.info('{0}/{1}/{2}/{3}'.format('L1B','GFZ','AOD1B',rl))
            #-- remote and local directory for exact data product
            local_dir = os.path.join(DIRECTORY,'AOD1B',rl)
            #-- check if AOD1B directory exists and recursively create if not
            os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
            #-- query CMR for dataset
            ids,urls,mtimes = gravity_toolkit.utilities.cmr(
                mission='grace', level='L1B', center='GFZ', release=rl,
                product='AOD1B', start_date='2002-01-01T00:00:00',
                provider='PODAAC', endpoint='data')
            #-- for each id, url and modification time
            for id,url,mtime in zip(ids,urls,mtimes):
                #-- retrieve GRACE/GRACE-FO files
                granule = gravity_toolkit.utilities.url_split(url)[-1]
                http_pull_file(url, mtime, os.path.join(local_dir,granule),
                    TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                    CHECKSUM=CHECKSUM, MODE=MODE)

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
                    #-- print string of exact data product
                    logging.info('{0} {1}/{2}/{3}'.format(mi, pr, rl, ds))
                    #-- query CMR for dataset
                    ids,urls,mtimes = gravity_toolkit.utilities.cmr(
                        mission=mi, center=pr, release=rl, product=ds,
                        version=VERSION[i], provider='PODAAC', endpoint='data')
                    #-- regular expression operator for data product
                    rx = compile_regex_pattern(pr, rl, ds, version=VERSION[i])
                    #-- for each id, url and modification time
                    for id,url,mtime in zip(ids,urls,mtimes):
                        #-- retrieve GRACE/GRACE-FO files
                        granule = gravity_toolkit.utilities.url_split(url)[-1]
                        http_pull_file(url, mtime, os.path.join(local_dir,granule),
                            TIMEOUT=TIMEOUT, LIST=LIST, CLOBBER=CLOBBER,
                            CHECKSUM=CHECKSUM, MODE=MODE)
                        #-- extend list of GRACE/GRACE-FO files with granule
                        grace_files.append(granule) if rx.match(granule) else None

                #-- outputting GRACE/GRACE-FO filenames to index
                with open(os.path.join(local_dir,'index.txt'),'w') as fid:
                    for fi in sorted(grace_files):
                        print('{0}'.format(fi), file=fid)
                #-- change permissions of index file
                os.chmod(os.path.join(local_dir,'index.txt'), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file, TIMEOUT=120,
    LIST=False, CLOBBER=False, CHECKSUM=False, MODE=0o775):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        local_hash = gravity_toolkit.utilities.get_hash(local_file)
        #-- Create and submit request.
        #-- There are a wide range of exceptions that can be thrown here
        #-- including HTTPError and URLError.
        req=gravity_toolkit.utilities.urllib2.Request(remote_file)
        resp=gravity_toolkit.utilities.urllib2.urlopen(req,timeout=TIMEOUT)
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(resp.read())
        remote_buffer.seek(0)
        #-- generate checksum hash for remote file
        remote_hash = gravity_toolkit.utilities.get_hash(remote_buffer)
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
        logging.info('{0} --> '.format(remote_file))
        logging.info('\t{0}{1}\n'.format(local_file,OVERWRITE))
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- copy bytes or transfer file
            if CHECKSUM and os.access(local_file, os.F_OK):
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
            else:
                #-- Create and submit request.
                #-- There are a range of exceptions that can be thrown here
                #-- including HTTPError and URLError.
                request = gravity_toolkit.utilities.urllib2.Request(remote_file)
                response = gravity_toolkit.utilities.urllib2.urlopen(request,
                    timeout=TIMEOUT)
                #-- copy remote file contents to local file
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(response, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Syncs GRACE/GRACE-FO and auxiliary data from the
            NASA JPL PO.DAAC Drive Server.
            Syncs GRACE/GRACE-FO Level-1b dealiasing products (AOD1B).
            Gets the latest technical note (TN) files.
            Gets the monthly GRACE/GRACE-FO newsletters.
            """
    )
    #-- command line parameters
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--webdav','-W',
        type=str, default=os.environ.get('PODAAC_PASSWORD'),
        help='WebDAV Password for JPL PO.DAAC Drive Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
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
        default=['RL06'], choices=['RL06'],
        help='GRACE/GRACE-FO data release')
    #-- GRACE/GRACE-FO data version
    parser.add_argument('--version','-v',
        metavar='VERSION', type=str, nargs='+',
        default=['0','1'], choices=['0','1','2','3'],
        help='GRACE/GRACE-FO Level-2 data version')
    #-- GRACE/GRACE-FO dealiasing products
    parser.add_argument('--aod1b','-a',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Level-1B dealiasing products')
    #-- GRACE/GRACE-FO newsletters
    parser.add_argument('--newsletters','-n',
        default=False, action='store_true',
        help='Sync GRACE/GRACE-FO Newsletters')
    #-- connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    #-- Output log file in form
    #-- PODAAC_sync_2002-04-01.log
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

    #-- JPL PO.DAAC drive hostname
    HOST = 'podaac-tools.jpl.nasa.gov'
    #-- get NASA Earthdata and JPL PO.DAAC drive credentials
    try:
        args.user,_,args.webdav = netrc.netrc(args.netrc).authenticators(HOST)
    except:
        #-- check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = 'Username for {0}: '.format(HOST)
            args.user = builtins.input(prompt)
        #-- enter WebDAV password securely from command-line
        if not args.webdav:
            prompt = 'Password for {0}@{1}: '.format(args.user,HOST)
            args.webdav = getpass.getpass(prompt)

    #-- build a urllib opener for PO.DAAC Drive
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user,args.webdav)

    #-- check internet connection before attempting to run program
    #-- check JPL PO.DAAC Drive credentials before attempting to run program
    DRIVE = 'https://{0}/drive/files'.format(HOST)
    if gravity_toolkit.utilities.check_credentials(DRIVE):
        podaac_grace_sync(args.directory, PROC=args.center,
            DREL=args.release, VERSION=args.version,
            NEWSLETTERS=args.newsletters, AOD1B=args.aod1b,
            TIMEOUT=args.timeout, LIST=args.list, LOG=args.log,
            CLOBBER=args.clobber, CHECKSUM=args.checksum,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
