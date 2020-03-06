#!/usr/bin/env python
u"""
podaac_grace_sync.py
Written by Tyler Sutterley (03/2020)

Syncs GRACE/GRACE-FO and auxilary data from the NASA JPL PO.DAAC Drive Server
Syncs CSR/GFZ/JPL files for RL04/RL05 GAA/GAB/GAC/GAD/GSM
    GAA and GAB are GFZ/JPL only
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
    podaac_grace_sync(<directory>, PROC, USER=<username>)
        or
    python podaac_grace_sync.py --user=<username>
    where <username> is your NASA Earthdata username

OUTPUTS:
    CSR RL04/RL05/RL06: GAC/GAD/GSM
    GFZ RL04/RL05: GAA/GAB/GAC/GAD/GSM
    JPL RL04/RL05: GAA/GAB/GAC/GAD/GSM
    Tellus degree one coefficients
    Technical notes for satellite laser ranging coefficients
    Technical notes for Release-05 atmospheric corrections
    Monthly GRACE newsletters

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user=X: username for NASA Earthdata Login
    -D X, --directory=X: working data directory
    -C X, --center=X: GRACE Processing Center
    -R X, --release=X: GRACE data releases to sync (RL05,RL06)
    --newsletters: sync GRACE newsletters
    -L, --list: print files to be transferred, but do not execute transfer
    -l, --log: output log of files downloaded
    --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode=X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        http://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

UPDATE HISTORY:
    Updated 03/2020 for public release
"""
from __future__ import print_function

import sys
import os
import re
import io
import ssl
import getopt
import shutil
import base64
import hashlib
import getpass
import builtins
import posixpath
import lxml.etree
import calendar, time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    import urllib2
else:
    from http.cookiejar import CookieJar
    import urllib.request as urllib2

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to https host for PO.DAAC
    try:
        HOST = posixpath.join('https://podaac-tools.jpl.nasa.gov','drive')
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
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

#-- PURPOSE: sync local GRACE/GRACE-FO files with JPL PO.DAAC drive server
def podaac_grace_sync(DIRECTORY, PROC, USER=None, PASSWORD=None, DREL=[],
    NEWSLETTERS=False, LOG=False, LIST=False, CLOBBER=False, CHECKSUM=False,
    MODE=None):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
    #-- RL04 has been moved on PO.DAAC to the retired directory
    remote_sub = {}
    remote_sub['RL04'] = 'retired'
    remote_sub['RL05'] = 'retired'
    remote_sub['RL06'] = ''
    #-- datasets for each processing center
    DSET = {}
    DSET['CSR'] = ['GAC', 'GAD', 'GSM']
    DSET['GFZ'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    DSET['JPL'] = ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']
    #-- remote subdirectories for newsletters (note capital for gracefo)
    newsletter_sub = {}
    newsletter_sub['grace'] = ['grace','docs','newsletters']
    newsletter_sub['gracefo'] = ['gracefo','docs','Newsletters']

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: PODAAC_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'PODAAC_sync_{0}.log'.format(today)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('PO.DAAC Sync Log ({0})'.format(today), file=fid1)
        print('CENTERS={0}'.format(','.join(PROC)), file=fid1)
        print('RELEASES={0}'.format(','.join(DREL)), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- https://docs.python.org/3/howto/urllib2.html#id5
    #-- create a password manager
    password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    #-- Add the username for NASA Earthdata Login system and password for WebDAV
    password_mgr.add_password(None, 'https://urs.earthdata.nasa.gov',
        USER, PASSWORD)
    #-- Encode username/password for request authorization headers
    base64_string = base64.b64encode('{0}:{1}'.format(USER, PASSWORD).encode())
    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()
    #-- Create cookie jar for storing cookies. This is used to store and return
    #-- the session cookie given to use by the data server (otherwise will just
    #-- keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(
        urllib2.HTTPBasicAuthHandler(password_mgr),
        urllib2.HTTPSHandler(context=ssl.SSLContext()),
        urllib2.HTTPCookieProcessor(cookie_jar))
    #-- add Authorization header to opener
    authorization_header = "Basic {0}".format(base64_string.decode())
    opener.addheaders = [("Authorization", authorization_header)]
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.

    #-- remote https server for GRACE data
    HOST = posixpath.join('https://podaac-tools.jpl.nasa.gov','drive','files')

    #-- DEGREE 1 COEFFICIENTS
    print('Degree 1 Coefficients:', file=fid1)
    remote_dir = posixpath.join(HOST,'allData','tellus','L2','degree_1')
    local_dir = os.path.join(DIRECTORY,'geocenter')
    #-- check if geocenter directory exists and recursively create if not
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    #-- open connection with PO.DAAC drive server at remote directory
    req = urllib2.Request(url=remote_dir)
    #-- read and parse request for files (find column names and modified dates)
    tree = lxml.etree.parse(urllib2.urlopen(req, timeout=20), parser)
    colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
    collastmod = tree.xpath('//tr/td[3]/text()')
    #-- TN-13 JPL degree 1 files
    #-- compile regular expression operator for remote files
    R1 = re.compile('TN-13_GEOC_(CSR|GFZ|JPL)_(.*?).txt', re.VERBOSE)
    remote_file_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
    #-- for each file on the remote server
    for i in remote_file_lines:
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colnames[i])
        local_file = os.path.join(local_dir,colnames[i])
        #-- get last modified date of file and convert into unix time
        remote_mtime = get_mtime(collastmod[i])
        http_pull_file(fid1, remote_file, remote_mtime, local_file,
            LIST, CLOBBER, CHECKSUM, MODE)
    #-- close request
    req = None

    #-- SLR C2,0 COEFFICIENTS
    print('C2,0 Coefficients:', file=fid1)
    remote_dir = posixpath.join(HOST,'allData','grace','docs')
    local_dir = os.path.expanduser(DIRECTORY)
    #-- open connection with PO.DAAC drive server at remote directory
    req = urllib2.Request(url=remote_dir)
    #-- read and parse request for files (find column names and modified dates)
    tree = lxml.etree.parse(urllib2.urlopen(req, timeout=20), parser)
    colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
    collastmod = tree.xpath('//tr/td[3]/text()')
    #-- compile regular expression operator for remote files
    R1 = re.compile('TN-(05|07|11)_C20_SLR.txt', re.VERBOSE)
    remote_file_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
    #-- for each file on the remote server
    for i in remote_file_lines:
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colnames[i])
        local_file = os.path.join(local_dir,colnames[i])
        #-- get last modified date of file and convert into unix time
        remote_mtime = get_mtime(collastmod[i])
        http_pull_file(fid1, remote_file, remote_mtime, local_file,
            LIST, CLOBBER, CHECKSUM, MODE)
    #-- close request
    req = None

    #-- SLR C3,0 COEFFICIENTS
    print('C3,0 Coefficients:', file=fid1)
    remote_dir = posixpath.join(HOST,'allData','gracefo','docs')
    local_dir = os.path.expanduser(DIRECTORY)
    #-- open connection with PO.DAAC drive server at remote directory
    req = urllib2.Request(url=remote_dir)
    #-- read and parse request for files (find column names and modified dates)
    tree = lxml.etree.parse(urllib2.urlopen(req, timeout=20), parser)
    colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
    collastmod = tree.xpath('//tr/td[3]/text()')
    #-- compile regular expression operator for remote files
    R1 = re.compile('TN-(14)_C30_C20_GSFC_SLR.txt', re.VERBOSE)
    remote_file_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
    #-- for each file on the remote server
    for i in remote_file_lines:
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colnames[i])
        local_file = os.path.join(local_dir,colnames[i])
        #-- get last modified date of file and convert into unix time
        remote_mtime = get_mtime(collastmod[i])
        http_pull_file(fid1, remote_file, remote_mtime, local_file,
            LIST, CLOBBER, CHECKSUM, MODE)
    #-- close request
    req = None

    #-- TN-08 GAE, TN-09 GAF and TN-10 GAG ECMWF atmosphere correction products
    print('TN-08 GAE, TN-09 GAF and TN-10 GAG products:', file=fid1)
    remote_dir = posixpath.join(HOST,'allData','grace','docs')
    local_dir = os.path.expanduser(DIRECTORY)
    ECMWF_files = []
    ECMWF_files.append('TN-08_GAE-2_2006032-2010031_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-09_GAF-2_2010032-2015131_0000_EIGEN_G---_0005.gz')
    ECMWF_files.append('TN-10_GAG-2_2015132-2099001_0000_EIGEN_G---_0005.gz')
    #-- open connection with PO.DAAC drive server at remote directory
    req = urllib2.Request(url=remote_dir)
    #-- read and parse request for files (find column names and modified dates)
    tree = lxml.etree.parse(urllib2.urlopen(req, timeout=20), parser)
    colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
    collastmod = tree.xpath('//tr/td[3]/text()')
    #-- compile regular expression operator for remote files
    R1 = re.compile('({0}|{1}|{2})'.format(*ECMWF_files), re.VERBOSE)
    remote_file_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
    #-- for each file on the remote server
    for i in remote_file_lines:
        #-- remote and local versions of the file
        remote_file = posixpath.join(remote_dir,colnames[i])
        local_file = os.path.join(local_dir,colnames[i])
        #-- get last modified date of file and convert into unix time
        remote_mtime = get_mtime(collastmod[i])
        http_pull_file(fid1, remote_file, remote_mtime, local_file,
            LIST, CLOBBER, CHECKSUM, MODE)
    #-- close request
    req = None

    #-- GRACE and GRACE-FO Newsletters
    if NEWSLETTERS:
        #-- local newsletter directory (place GRACE and GRACE-FO together)
        local_dir = os.path.join(DIRECTORY,'newsletters')
        #-- check if newsletters directory exists and recursively create if not
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- for each mission
        for MISSION,NAME in zip(['grace','gracefo'],['GRACE','GRACE_FO']):
            print('{0} Newsletters:'.format(NAME.replace('_','-')), file=fid1)
            #-- open connection with PO.DAAC drive server at remote directory
            remote_dir = posixpath.join(HOST,'allData',*newsletter_sub[MISSION])
            req = urllib2.Request(url=remote_dir)
            #-- read and parse request for files (find names and modified dates)
            tree = lxml.etree.parse(urllib2.urlopen(req, timeout=20), parser)
            colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
            collastmod = tree.xpath('//tr/td[3]/text()')
            #-- compile regular expression operator for remote files
            R1 = re.compile('{0}_SDS_NL_(\d+).pdf'.format(NAME), re.VERBOSE)
            remote_file_lines = [i for i,f in enumerate(colnames) if R1.match(f)]
            #-- for each file on the remote server
            for i in remote_file_lines:
                #-- remote and local versions of the file
                remote_file = posixpath.join(remote_dir,colnames[i])
                local_file = os.path.join(local_dir,colnames[i])
                #-- get last modified date of file and convert into unix time
                remote_mtime = get_mtime(collastmod[i])
                http_pull_file(fid1, remote_file, remote_mtime, local_file,
                    LIST, CLOBBER, CHECKSUM, MODE)
            #-- close request
            req = None

    #-- GRACE DATA
    #-- PROCESSING CENTERS (CSR, GFZ, JPL)
    print('GRACE L2 Global Spherical Harmonics:', file=fid1)
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
            remote_dir = posixpath.join(HOST,'allData','grace',remote_sub[rl],
                'L2', pr, drel_str)
            #-- open connection with PO.DAAC drive server at remote directory
            req = urllib2.Request(url=remote_dir)
            #-- read and parse request for files (names and modified dates)
            tree = lxml.etree.parse(urllib2.urlopen(req,timeout=20), parser)
            colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
            collastmod = tree.xpath('//tr/td[3]/text()')
            #-- DATA PRODUCTS (GAC, GAD, GSM, GAA, GAB)
            for ds in DSET[pr]:
                #-- print string of exact data product
                print('GRACE {0}/{1}/{2}'.format(pr, drel_str, ds), file=fid1)
                #-- local directory for exact data product
                local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                #-- check if directory exists and recursively create if not
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir,MODE)
                #-- compile regular expression operator to find GRACE files
                R1 = re.compile('({0}-(.*?)(gz|txt|dif))'.format(ds),re.VERBOSE)
                line = [i for i,f in enumerate(colnames) if R1.match(f)]
                #-- for each file on the remote server
                for i in line:
                    #-- remote and local versions of the file
                    remote_file = posixpath.join(remote_dir,colnames[i])
                    local_file = os.path.join(local_dir,colnames[i])
                    #-- get last modified date and convert into unix time
                    remote_mtime = get_mtime(collastmod[i])
                    http_pull_file(fid1, remote_file, remote_mtime, local_file,
                        LIST, CLOBBER, CHECKSUM, MODE)

    #-- GRACE-FO DATA
    #-- PROCESSING CENTERS (CSR, GFZ, JPL)
    #-- GRACE-FO data are stored separately for each year
    print('GRACE-FO L2 Global Spherical Harmonics:', file=fid1)
    for pr in PROC:
        #-- DATA RELEASES (RL06)
        valid_gracefo_releases = [d for d in DREL if d not in ('RL04','RL05')]
        for rl in valid_gracefo_releases:
            #-- remote directory for data release
            remote_dir = posixpath.join(HOST,'allData','gracefo',remote_sub[rl],
                'L2', pr, rl)
            #-- open connection with PO.DAAC drive server at remote directory
            req = urllib2.Request(url=remote_dir)
            #-- read and parse request for files (names and modified dates)
            tree = lxml.etree.parse(urllib2.urlopen(req,timeout=20), parser)
            colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
            collastmod = tree.xpath('//tr/td[3]/text()')
            years = [d for i,d in enumerate(colnames) if re.match(r'\d{4}',d)]
            for yr in years:
                #-- open connection with PO.DAAC drive server at remote directory
                req = urllib2.Request(url=posixpath.join(remote_dir,yr))
                #-- read and parse request for files (names and modified dates)
                tree = lxml.etree.parse(urllib2.urlopen(req,timeout=20), parser)
                colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
                collastmod = tree.xpath('//tr/td[3]/text()')
                #-- DATA PRODUCTS (GAC, GAD, GSM, GAA, GAB)
                for ds in DSET[pr]:
                    #-- print string of exact data product
                    args = (pr, rl, ds, yr)
                    print('GRACE-FO {0}/{1}/{2}/{3}'.format(*args), file=fid1)
                    #-- local directory for exact data product
                    local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                    #-- check if directory exists and recursively create if not
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir,MODE)
                    #-- compile regular expression operator to find GRACE files
                    R1 = re.compile('({0}-(.*?)(gz|txt|dif))'.format(ds))
                    line = [i for i,f in enumerate(colnames) if R1.match(f)]
                    #-- for each file on the remote server
                    for i in line:
                        #-- remote and local versions of the file
                        remote_file = posixpath.join(remote_dir,yr,colnames[i])
                        local_file = os.path.join(local_dir,colnames[i])
                        #-- get last modified date and convert into unix time
                        remote_mtime = get_mtime(collastmod[i])
                        http_pull_file(fid1, remote_file, remote_mtime,
                            local_file, LIST, CLOBBER, CHECKSUM, MODE)

    #-- create index file for GRACE/GRACE-FO L2 Spherical Harmonic Data
    #-- PROCESSING CENTERS (CSR, GFZ, JPL)
    for pr in PROC:
        #-- DATA RELEASES (RL04, RL05, RL06)
        for rl in DREL:
            #-- DATA PRODUCTS (GAC, GAD, GSM, GAA, GAB)
            for ds in DSET[pr]:
                #-- local directory for exact data product
                local_dir = os.path.join(DIRECTORY, pr, rl, ds)
                #-- Create an index file for each GRACE product
                #-- finding all dataset files *.gz in directory
                rx = compile_regex_pattern(pr, rl, ds)
                #-- find local GRACE files to create index
                grace_files=[fi for fi in os.listdir(local_dir) if rx.match(fi)]
                #-- outputting GRACE filenames to index
                with open(os.path.join(local_dir,'index.txt'),'w') as fid:
                    for fi in sorted(grace_files):
                        print('{0}'.format(fi), file=fid)
                #-- change permissions of index file
                os.chmod(os.path.join(local_dir,'index.txt'), MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(fid, remote_file, remote_mtime, local_file, LIST, CLOBBER,
    CHECKSUM, MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(local_file, 'rb') as local_buffer:
            local_hash = hashlib.md5(local_buffer.read()).hexdigest()
        #-- Create and submit request.
        #-- There are a wide range of exceptions that can be thrown here
        #-- including HTTPError and URLError.
        request = urllib2.Request(remote_file)
        response = urllib2.urlopen(request, timeout=20)
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(response.read())
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
        print('{0} --> '.format(remote_file), file=fid)
        print('\t{0}{1}\n'.format(local_file,OVERWRITE), file=fid)
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
                #-- There are a wide range of exceptions that can be thrown here
                #-- including HTTPError and URLError.
                request = urllib2.Request(remote_file)
                response = urllib2.urlopen(request, timeout=20)
                #-- copy contents to local file using chunked transfer encoding
                #-- transfer should work properly with ascii and binary formats
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(response, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- PURPOSE: returns the Unix timestamp value for a modification time
def get_mtime(collastmod, FORMAT='%Y-%m-%d %H:%M:%S'):
    lastmodtime = time.strptime(collastmod.rstrip(), FORMAT)
    return calendar.timegm(lastmodtime)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -U X, --user=X\t\tUsername for NASA Earthdata Login')
    print(' -D X, --directory=X\tWorking Data Directory')
    print(' -C X, --center=X\tGRACE Processing Center (CSR,GFZ,JPL)')
    print(' -R X, --release=X\tGRACE data releases to sync (RL05,RL06)')
    print(' --newsletters\t\tSync GRACE Newsletters')
    print(' -L, --list\t\tOnly print files that are to be transferred')
    print(' --clobber\t\tOverwrite existing data in transfer')
    print(' --checksum\t\tCompare hashes to check if overwriting existing data')
    print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'PODAAC_sync_{0}.log'.format(today)
    print('    Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls podaac_grace_sync()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','user=','directory=','list','log','center=',
        'release=','newsletters','clobber','checksum','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hU:D:lLC:R:M:',long_options)

    #-- command line parameters
    USER = ''
    #-- working data directory
    DIRECTORY = os.getcwd()
    LIST = False
    LOG = False
    CLOBBER = False
    #-- GRACE Processing Centers to run
    PROC = ['CSR', 'GFZ', 'JPL']
    #-- Data release
    DREL = ['RL05']
    NEWSLETTERS = False
    #-- Use hash for determining whether or not to overwrite
    CHECKSUM = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-U","--user"):
            USER = arg
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("-L","--list"):
            LIST = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("--clobber"):
            CLOBBER = True
        elif opt in ("-C","--center"):
            PROC = arg.upper().split(',')
        elif opt in ("-R","--release"):
            DREL = arg.upper().split(',')
        elif (opt in '--newsletters'):
            NEWSLETTERS = True
        elif opt in ("--checksum"):
            CHECKSUM = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- JPL PO.DAAC drive hostname
    HOST = 'podaac-tools.jpl.nasa.gov'
    #-- check that NASA Earthdata credentials were entered
    if not USER:
        USER = builtins.input('Username for {0}: '.format(HOST))
    #-- enter password securely from command-line
    PASSWORD = getpass.getpass('WebDAV Password for {0}@{1}: '.format(USER,HOST))

    #-- check internet connection before attempting to run program
    if check_connection():
        podaac_grace_sync(DIRECTORY, PROC, USER=USER, PASSWORD=PASSWORD,
            DREL=DREL, NEWSLETTERS=NEWSLETTERS, LIST=LIST, LOG=LOG,
            CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
