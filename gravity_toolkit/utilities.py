"""
utilities.py
Written by Tyler Sutterley (12/2020)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python (https://pypi.python.org/pypi/lxml)

UPDATE HISTORY:
    Updated 12/2020: added ICGEM list for static models
        added figshare geocenter download for Sutterley and Velicogna files
        added download for satellite laser ranging (SLR) files from UTCSR
        added file object keyword for downloads if verbose printing to file
        renamed podaac_list() and from_podaac() to drive_list() and from_drive()
        added username and password to ftp functions. added ftp connection check
    Updated 09/2020: copy from http and https to bytesIO object in chunks
        use netrc credentials if not entered from PO.DAAC functions
        generalize build opener function for different Earthdata instances
    Updated 08/2020: add PO.DAAC Drive opener, login and download functions
    Written 08/2020
"""
from __future__ import print_function, division

import sys
import os
import re
import io
import ssl
import json
import netrc
import ftplib
import shutil
import base64
import socket
import inspect
import hashlib
import posixpath
import lxml.etree
import calendar,time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
    import urllib.request as urllib2

def get_data_path(relpath):
    """
    Get the absolute path within a package from a relative path

    Arguments
    ---------
    relpath: relative path
    """
    #-- current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    if isinstance(relpath,list):
        #-- use *splat operator to extract from list
        return os.path.join(filepath,*relpath)
    elif isinstance(relpath,str):
        return os.path.join(filepath,relpath)

#-- PURPOSE: get the MD5 hash value of a file
def get_hash(local):
    """
    Get the MD5 hash value from a local file or BytesIO object

    Arguments
    ---------
    local: BytesIO object or path to file
    """
    #-- check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        return hashlib.md5(local.getvalue()).hexdigest()
    elif os.access(os.path.expanduser(local),os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(os.path.expanduser(local), 'rb') as local_buffer:
            return hashlib.md5(local_buffer.read()).hexdigest()
    else:
        return ''

#-- PURPOSE: recursively split a url path
def url_split(s):
    """
    Recursively split a url path into a list

    Arguments
    ---------
    s: url string
    """
    head, tail = posixpath.split(s)
    if head in ('http:','https:'):
        return s,
    elif head in ('', posixpath.sep):
        return tail,
    return url_split(head) + (tail,)

#-- PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S'):
    """
    Get the Unix timestamp value for a formatted date string

    Arguments
    ---------
    time_string: formatted time string to parse

    Keyword arguments
    -----------------
    format: format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        return None
    else:
        return calendar.timegm(parsed_time)

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(value):
    """
    Rounds a number to an even number less than or equal to original

    Arguments
    ---------
    value: number to be rounded
    """
    return 2*int(value//2)

#-- PURPOSE: make a copy of a file with all system information
def copy(source, destination, verbose=False, move=False):
    """
    Copy or move a file with all system information

    Arguments
    ---------
    source: source file
    destination: copied destination file

    Keyword arguments
    -----------------
    verbose: print file transfer information
    move: remove the source file
    """
    source = os.path.abspath(os.path.expanduser(source))
    destination = os.path.abspath(os.path.expanduser(destination))
    print('{0} -->\n\t{1}'.format(source,destination)) if verbose else None
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    if move:
        os.remove(source)

#-- PURPOSE: check ftp connection
def check_ftp_connection(HOST,username=None,password=None):
    """
    Check internet connection with ftp host

    Arguments
    ---------
    HOST: remote ftp host

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    """
    #-- attempt to connect to ftp host
    try:
        f = ftplib.FTP(HOST)
        f.login(username, password)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        return True

#-- PURPOSE: list a directory on a ftp host
def ftp_list(HOST,username=None,password=None,timeout=None,
    basename=False,pattern=None,sort=False):
    """
    List a directory on a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    timeout: timeout in seconds for blocking operations
    basename: return the file or directory basename instead of the full path
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    output: list of items in a directory
    mtimes: list of last modification times for items in the directory
    """
    #-- try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    else:
        ftp.login(username,password)
        #-- list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        #-- get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        #-- iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                #-- try sending modification time command
                mdtm = ftp.sendcmd('MDTM {0}'.format(f))
            except ftplib.error_perm:
                #-- directories will return with an error
                pass
            else:
                #-- convert the modification time into unix time
                mtimes[i] = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        #-- reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            #-- reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            #-- sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- close the ftp connection
        ftp.close()
        #-- return the list of items and last modified times
        return (output,mtimes)

#-- PURPOSE: download a file from a ftp host
def from_ftp(HOST,username=None,password=None,timeout=None,local=None,
    hash='',chunk=8192,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    timeout: timeout in seconds for blocking operations
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- try downloading from ftp
    try:
        #-- try to connect to ftp host
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    else:
        ftp.login(username,password)
        #-- remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary('RETR {0}'.format(ftp_remote_path),
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd('MDTM {0}'.format(ftp_remote_path))
        remote_mtime = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            if verbose:
                args = (posixpath.join(*HOST),local)
                print('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
            #-- keep remote modification time of file and local access time
            os.utime(local, (os.stat(local).st_atime, remote_mtime))
        #-- close the ftp connection
        ftp.close()
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: check internet connection
def check_connection(HOST):
    """
    Check internet connection with http host

    Arguments
    ---------
    HOST: remote http host
    """
    #-- attempt to connect to http host
    try:
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: download a file from a http host
def from_http(HOST,timeout=None,context=ssl.SSLContext(),local=None,hash='',
    chunk=16384,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a http host

    Arguments
    ---------
    HOST: remote http host path split as list

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    context: SSL context for url opener object
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- try downloading from http
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout,context=context)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            if verbose:
                args = (posixpath.join(*HOST),local)
                print('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: "login" to JPL PO.DAAC Drive with supplied credentials
def build_opener(username, password, context=ssl.SSLContext(),
    password_manager=False, get_ca_certs=False, redirect=False,
    authorization_header=True, urs='https://urs.earthdata.nasa.gov'):
    """
    build urllib opener for NASA Earthdata or JPL PO.DAAC Drive with
    supplied credentials

    Arguments
    ---------
    username: NASA Earthdata username
    password: NASA Earthdata or JPL PO.DAAC WebDAV password

    Keyword arguments
    -----------------
    context: SSL context for opener object
    password_manager: create password manager context using default realm
    get_ca_certs: get list of loaded “certification authority” certificates
    redirect: create redirect handler object
    authorization_header: add base64 encoded authorization header to opener
    urs: Earthdata login URS 3 host
    """
    #-- https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    #-- create a password manager
    if password_manager:
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        #-- Add the username and password for NASA Earthdata Login system
        password_mgr.add_password(None,urs,username,password)
        handler.append(urllib2.HTTPBasicAuthHandler(password_mgr))
    #-- Create cookie jar for storing cookies. This is used to store and return
    #-- the session cookie given to use by the data server (otherwise will just
    #-- keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    #-- SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    handler.append(urllib2.HTTPSHandler(context=context))
    #-- redirect handler
    if redirect:
        handler.append(urllib2.HTTPRedirectHandler())
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    #-- Encode username/password for request authorization headers
    #-- add Authorization header to opener
    if authorization_header:
        b64 = base64.b64encode('{0}:{1}'.format(username,password).encode())
        opener.addheaders = [("Authorization","Basic {0}".format(b64.decode()))]
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.

#-- PURPOSE: check that entered JPL PO.DAAC/ECCO Drive credentials are valid
def check_credentials(HOST='https://podaac-tools.jpl.nasa.gov'):
    """
    Check that entered JPL PO.DAAC or ECCO Drive credentials are valid

    Keyword arguments
    -----------------
    HOST: PO.DAAC or ECCO Drive host
    """
    try:
        request = urllib2.Request(url=posixpath.join(HOST,'drive','files'))
        response = urllib2.urlopen(request, timeout=20)
    except urllib2.HTTPError:
        raise RuntimeError('Check your JPL PO.DAAC Drive credentials')
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: list a directory on JPL PO.DAAC/ECCO Drive https server
def drive_list(HOST,username=None,password=None,build=True,timeout=None,
    urs='podaac-tools.jpl.nasa.gov',parser=lxml.etree.HTMLParser(),
    pattern='',sort=False):
    """
    List a directory on JPL PO.DAAC or ECCO Drive

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    username: NASA Earthdata username
    password: JPL PO.DAAC Drive WebDAV password
    build: Build opener and check WebDAV credentials
    timeout: timeout in seconds for blocking operations
    urs: JPL PO.DAAC or ECCO login URS 3 host
    parser: HTML parser for lxml
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    colnames: list of column names in a directory
    collastmod: list of last modification times for items in the directory
    """
    #-- use netrc credentials
    if build and not (username or password):
        username,login,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- try listing from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        tree = lxml.etree.parse(urllib2.urlopen(request,timeout=timeout),parser)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- read and parse request for files (column names and modified times)
        colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
        #-- get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i) for i in tree.xpath('//tr/td[3]/text()')]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            #-- reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            #-- sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- return the list of column names and last modified times
        return (colnames,collastmod)

#-- PURPOSE: download a file from a PO.DAAC/Ecco Drive https server
def from_drive(HOST,username=None,password=None,build=True,timeout=None,
    urs='podaac-tools.jpl.nasa.gov',local=None,hash='',chunk=16384,
    verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a JPL PO.DAAC or ECCO Drive https server

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    username: NASA Earthdata username
    password: JPL PO.DAAC Drive WebDAV password
    build: Build opener and check WebDAV credentials
    timeout: timeout in seconds for blocking operations
    urs: JPL PO.DAAC or ECCO login URS 3 host
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- use netrc credentials
    if build and not (username or password):
        username,login,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- try downloading from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            if verbose:
                args = (posixpath.join(*HOST),local)
                print('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: download geocenter files from Sutterley and Velicogna (2019)
#-- https://doi.org/10.3390/rs11182108
#-- https://doi.org/10.6084/m9.figshare.7388540
def from_figshare(directory,article='7388540',timeout=None,
    context=ssl.SSLContext(),chunk=16384,verbose=False,fid=sys.stdout,
    pattern=r'(CSR|GFZ|JPL)_(RL\d+)_(.*?)_SLF_iter.txt$',mode=0o775):
    """
    Download Sutterley and Velicogna (2019) geocenter files from figshare

    Arguments
    ---------
    directory: download directory

    Keyword arguments
    -----------------
    article: figshare article number
    timeout: timeout in seconds for blocking operations
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    pattern: regular expression pattern for reducing list
    mode: permissions mode of output local file
    """
    #-- figshare host
    HOST=['https://api.figshare.com','v2','articles',article]
    #-- recursively create directory if non-existent
    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.access(os.path.join(directory,'geocenter'), os.F_OK):
        os.makedirs(os.path.join(directory,'geocenter'), mode)
    #-- Create and submit request.
    request = urllib2.Request(posixpath.join(*HOST))
    response = urllib2.urlopen(request,timeout=timeout,context=context)
    resp = json.loads(response.read())
    #-- reduce list of geocenter files
    geocenter_files = [f for f in resp['files'] if re.match(pattern,f['name'])]
    for f in geocenter_files:
        #-- download geocenter file
        original_md5 = get_hash(os.path.join(directory,'geocenter',f['name']))
        from_http(url_split(f['download_url']),timeout=timeout,context=context,
            local=os.path.join(directory,'geocenter',f['name']),
            hash=original_md5,chunk=chunk,verbose=verbose,fid=fid,mode=mode)
        #-- verify MD5 checksums
        computed_md5 = get_hash(os.path.join(directory,'geocenter',f['name']))
        if (computed_md5 != f['supplied_md5']):
            raise Exception('Checksum mismatch: {0}'.format(f['download_url']))

#-- PURPOSE: download satellite laser ranging files from CSR
#-- http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
#-- http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s
def from_csr(directory,timeout=None,context=ssl.SSLContext(),
    chunk=16384,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download satellite laser ranging (SLR) files from the
        University of Texas Center for Space Research (UTCSR)

    Arguments
    ---------
    directory: download directory

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    context: SSL context for url opener object
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file
    """
    #-- create directory if non-existent
    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.access(os.path.join(directory,'geocenter'), os.F_OK):
        os.makedirs(os.path.join(directory,'geocenter'), mode)
    #-- download SLR 5x5 file
    HOST = ['http://download.csr.utexas.edu','pub','slr','degree_5',
        'CSR_Monthly_5x5_Gravity_Harmonics.txt']
    original_md5 = get_hash(os.path.join(directory,HOST[-1]))
    from_http(HOST,timeout=timeout,context=context,
        local=os.path.join(directory,HOST[-1]),hash=original_md5,
        chunk=chunk,verbose=verbose,fid=fid,mode=mode)
    #-- download CF-CM SLR geocenter file
    HOST = ['http://download.csr.utexas.edu','pub','slr','geocenter',
        'GCN_L1_L2_30d_CF-CM.txt']
    original_md5 = get_hash(os.path.join(directory,'geocenter',HOST[-1]))
    from_http(HOST,timeout=timeout,context=context,
        local=os.path.join(directory,'geocenter',HOST[-1]),hash=original_md5,
        chunk=chunk,verbose=verbose,fid=fid,mode=mode)
    #-- download updated SLR geocenter file from Minkang Cheng
    HOST = ['http://download.csr.utexas.edu','outgoing','cheng',
        'gct2est.220_5s']
    original_md5 = get_hash(os.path.join(directory,'geocenter',HOST[-1]))
    from_http(HOST,timeout=timeout,context=context,
        local=os.path.join(directory,'geocenter',HOST[-1]),hash=original_md5,
        chunk=chunk,verbose=verbose,fid=fid,mode=mode)

#-- PURPOSE: list a directory on the GFZ ICGEM https server
#-- http://icgem.gfz-potsdam.de
def icgem_list(host='http://icgem.gfz-potsdam.de/tom_longtime',timeout=None,
    parser=lxml.etree.HTMLParser()):
    """
    Parse the table of static gravity field models on the GFZ
    International Centre for Global Earth Models (ICGEM) server

    Keyword arguments
    -----------------
    host: url for the GFZ ICGEM gravity field table
    timeout: timeout in seconds for blocking operations
    parser: HTML parser for lxml

    Returns
    -------
    colfiles: dictionary of static file urls mapped by field name
    """
    #-- try listing from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(host)
        tree = lxml.etree.parse(urllib2.urlopen(request,timeout=timeout),parser)
    except:
        raise Exception('List error from {0}'.format(host))
    else:
        #-- read and parse request for files
        colfiles = tree.xpath('//td[@class="tom-cell-modelfile"]//a/@href')
        #-- reduce list of files to find gfc files
        #-- return the dict of model files mapped by name
        return {re.findall('(.*?).gfc',posixpath.basename(f)).pop():url_split(f)
            for i,f in enumerate(colfiles) if re.search('gfc$',f)}
