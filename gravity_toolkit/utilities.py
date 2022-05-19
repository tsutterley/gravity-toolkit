#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (05/2022)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

UPDATE HISTORY:
    Updated 05/2022: function for extracting bucket name from presigned url
    Updated 04/2022: updated docstrings to numpy documentation format
        update CMR queries to prepare for version 1 of RL06
    Updated 03/2022: add NASA Common Metadata Repository (CMR) queries
        added attempt login function to recursively check credentials
    Updated 11/2021: add CSR satellite laser ranging oblateness file
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: added generic list from Apache http server
    Updated 07/2021: added unique filename opener for log files
    Updated 06/2021: add parser for converting file lines to arguments
    Updated 05/2021: download GFZ satellite laser ranging and GravIS files
    Updated 04/2021: download CSR SLR figure axis and azimuthal dependence files
    Updated 03/2021: added sha1 option for retrieving file hashes
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
import boto3
import netrc
import ftplib
import shutil
import base64
import socket
import getpass
import inspect
import hashlib
import logging
import builtins
import dateutil
import warnings
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

#-- PURPOSE: get absolute path within a package from a relative path
def get_data_path(relpath):
    """
    Get the absolute path within a package from a relative path

    Parameters
    ----------
    relpath: str,
        relative path
    """
    #-- current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    if isinstance(relpath,list):
        #-- use *splat operator to extract from list
        return os.path.join(filepath,*relpath)
    elif isinstance(relpath,str):
        return os.path.join(filepath,relpath)

#-- PURPOSE: get the hash value of a file
def get_hash(local, algorithm='MD5'):
    """
    Get the hash value from a local file or BytesIO object

    Parameters
    ----------
    local: obj or str
        BytesIO object or path to file
    algorithm: str, default 'MD5'
        hashing algorithm for checksum validation

            - ``'MD5'``: Message Digest
            - ``'sha1'``: Secure Hash Algorithm
    """
    #-- check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        if (algorithm == 'MD5'):
            return hashlib.md5(local.getvalue()).hexdigest()
        elif (algorithm == 'sha1'):
            return hashlib.sha1(local.getvalue()).hexdigest()
    elif os.access(os.path.expanduser(local),os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(os.path.expanduser(local), 'rb') as local_buffer:
            #-- generate checksum hash for a given type
            if (algorithm == 'MD5'):
                return hashlib.md5(local_buffer.read()).hexdigest()
            elif (algorithm == 'sha1'):
                return hashlib.sha1(local_buffer.read()).hexdigest()
    else:
        return ''

#-- PURPOSE: recursively split a url path
def url_split(s):
    """
    Recursively split a url path into a list

    Parameters
    ----------
    s: str
        url string
    """
    head, tail = posixpath.split(s)
    if head in ('http:','https:','ftp:','s3:'):
        return s,
    elif head in ('', posixpath.sep):
        return tail,
    return url_split(head) + (tail,)

#-- PURPOSE: convert file lines to arguments
def convert_arg_line_to_args(arg_line):
    """
    Convert file lines to arguments

    Parameters
    ----------
    arg_line: str
        line string containing a single argument and/or comments
    """
    #-- remove commented lines and after argument comments
    for arg in re.sub(r'\#(.*?)$',r'',arg_line).split():
        if not arg.strip():
            continue
        yield arg

#-- PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S'):
    """
    Get the Unix timestamp value for a formatted date string

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    format: str, default '%Y-%m-%d %H:%M:%S'
        format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        pass
    else:
        return calendar.timegm(parsed_time)
    #-- try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.timestamp()

#-- PURPOSE: output a time string in isoformat
def isoformat(time_string):
    """
    Reformat a date string to ISO formatting

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    """
    #-- try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.isoformat()

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(value):
    """
    Rounds a number to an even number less than or equal to original

    Parameters
    ----------
    value: float
        number to be rounded
    """
    return 2*int(value//2)

#-- PURPOSE: rounds a number upward to its nearest integer
def ceil(value):
    """
    Rounds a number upward to its nearest integer

    Parameters
    ----------
    value: float
        number to be rounded upward
    """
    return -int(-value//1)

#-- PURPOSE: make a copy of a file with all system information
def copy(source, destination, move=False, **kwargs):
    """
    Copy or move a file with all system information

    Parameters
    ----------
    source: str
        source file
    destination: str
        copied destination file
    move: bool, default False
        remove the source file
    """
    source = os.path.abspath(os.path.expanduser(source))
    destination = os.path.abspath(os.path.expanduser(destination))
    #-- log source and destination
    logging.info('{0} -->\n\t{1}'.format(source,destination))
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    if move:
        os.remove(source)

#-- PURPOSE: open a unique file adding a numerical instance if existing
def create_unique_file(filename):
    """
    Open a unique file adding a numerical instance if existing

    Parameters
    ----------
    filename: str
        full path to output file
    """
    #-- split filename into fileBasename and fileExtension
    fileBasename, fileExtension = os.path.splitext(filename)
    #-- create counter to add to the end of the filename if existing
    counter = 1
    while counter:
        try:
            #-- open file descriptor only if the file doesn't exist
            fd = os.open(filename, os.O_CREAT | os.O_EXCL | os.O_RDWR)
        except OSError:
            pass
        else:
            return os.fdopen(fd, 'w+')
        #-- new filename adds counter the between fileBasename and fileExtension
        filename = '{0}_{1:d}{2}'.format(fileBasename, counter, fileExtension)
        counter += 1

#-- PURPOSE: check ftp connection
def check_ftp_connection(HOST, username=None, password=None):
    """
    Check internet connection with ftp host

    Parameters
    ----------
    HOST: str
        remote ftp host
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
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
def ftp_list(HOST, username=None, password=None, timeout=None,
    basename=False, pattern=None, sort=False):
    """
    List a directory on a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path split as list
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    basename: bool, default False
        return the file or directory basename instead of the full path
    pattern: str or NoneType, default None
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    output: list
        items in a directory
    mtimes: list
        last modification times for items in the directory
    """
    #-- verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
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
def from_ftp(HOST, username=None, password=None, timeout=None,
    local=None, hash='', chunk=8192, verbose=False, fid=sys.stdout,
    mode=0o775):
    """
    Download a file from a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    local: str or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 8192
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
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
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args))
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

    Parameters
    ----------
    HOST: str
        remote http host
    """
    #-- attempt to connect to http host
    try:
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: list a directory on an Apache http Server
def http_list(HOST, timeout=None, context=ssl.SSLContext(),
    parser=lxml.etree.HTMLParser(), format='%Y-%m-%d %H:%M',
    pattern='', sort=False):
    """
    List a directory on an Apache http Server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for lxml
    format: str, default '%Y-%m-%d %H:%M'
        format for input time string
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    #-- verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try listing from http
    try:
        #-- Create and submit request.
        request=urllib2.Request(posixpath.join(*HOST))
        response=urllib2.urlopen(request,timeout=timeout,context=context)
    except (urllib2.HTTPError, urllib2.URLError) as e:
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response,parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        #-- get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
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

#-- PURPOSE: download a file from a http host
def from_http(HOST, timeout=None, context=ssl.SSLContext(),
    local=None, hash='', chunk=16384, verbose=False, fid=sys.stdout,
    mode=0o775):
    """
    Download a file from a http host

    Parameters
    ----------
    HOST: str or list
        remote http host path split as list
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    local: str or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
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
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args))
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: attempt to build an opener with netrc
def attempt_login(urs, context=ssl.SSLContext(),
    password_manager=True, get_ca_certs=False, redirect=False,
    authorization_header=True, **kwargs):
    """
    attempt to build a urllib opener for NASA Earthdata

    Parameters
    ----------
    urs: str
        Earthdata login URS 3 host
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    password_manager: bool, default True
        Create password manager context using default realm
    get_ca_certs: bool, default False
        Get list of loaded “certification authority” certificates
    redirect: bool, default False
        Create redirect handler object
    authorization_header: bool, default True
        Add base64 encoded authorization header to opener
    username: str, default from environmental variable
        NASA Earthdata username
    password: str, default from environmental variable
        NASA Earthdata password
    retries: int, default 5
        number of retry attempts
    netrc: str, default ~/.netrc
        path to .netrc file for authentication

    Returns
    -------
    opener: obj
        OpenerDirector instance
    """
    # set default keyword arguments
    kwargs.setdefault('username', os.environ.get('EARTHDATA_USERNAME'))
    kwargs.setdefault('password', os.environ.get('EARTHDATA_PASSWORD'))
    kwargs.setdefault('retries', 5)
    kwargs.setdefault('netrc', os.path.expanduser('~/.netrc'))
    try:
        # only necessary on jupyterhub
        os.chmod(kwargs['netrc'], 0o600)
        # try retrieving credentials from netrc
        username, _, password = netrc.netrc(kwargs['netrc']).authenticators(urs)
    except Exception as e:
        # try retrieving credentials from environmental variables
        username, password = (kwargs['username'], kwargs['password'])
        pass
    # if username or password are not available
    if not username:
        username = builtins.input('Username for {0}: '.format(urs))
    if not password:
        prompt = 'Password for {0}@{1}: '.format(username, urs)
        password = getpass.getpass(prompt=prompt)
    # for each retry
    for retry in range(kwargs['retries']):
        # build an opener for urs with credentials
        opener = build_opener(username, password,
            context=context,
            password_manager=password_manager,
            get_ca_certs=get_ca_certs,
            redirect=redirect,
            authorization_header=authorization_header,
            urs=urs)
        # try logging in by check credentials
        HOST = 'https://archive.podaac.earthdata.nasa.gov/s3credentials'
        try:
            check_credentials(HOST)
        except Exception as e:
            pass
        else:
            return opener
        # reattempt login
        username = builtins.input('Username for {0}: '.format(urs))
        password = getpass.getpass(prompt=prompt)
    # reached end of available retries
    raise RuntimeError('End of Retries: Check NASA Earthdata credentials')

#-- PURPOSE: "login" to NASA Earthdata with supplied credentials
def build_opener(username, password, context=ssl.SSLContext(),
    password_manager=False, get_ca_certs=False, redirect=False,
    authorization_header=True, urs='https://urs.earthdata.nasa.gov'):
    """
    build urllib opener for NASA Earthdata with supplied credentials

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    password_manager: bool, default False
        Create password manager context using default realm
    get_ca_certs: bool, default False
        Get list of loaded “certification authority” certificates
    redirect: bool, default False
        Create redirect handler object
    authorization_header: bool, default True
        Add base64 encoded authorization header to opener
    urs: str, default 'https://urs.earthdata.nasa.gov'
        Earthdata login URS 3 host

    Returns
    -------
    opener: obj
        OpenerDirector instance
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
    return opener

#-- PURPOSE: get AWS s3 client for PO.DAAC Cumulus
def s3_client(HOST=None, timeout=None, region_name='us-west-2'):
    """
    Get AWS s3 client for PO.DAAC Cumulus

    Parameters
    ----------
    HOST: str
        PO.DAAC or ECCO AWS S3 credential host
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    region_name: str, default 'us-west-2'
        AWS region name

    Returns
    -------
    client: obj
        AWS s3 client for PO.DAAC Cumulus
    """
    request = urllib2.Request(HOST)
    response = urllib2.urlopen(request, timeout=timeout)
    cumulus = json.loads(response.read())
    #-- get AWS client object
    client = boto3.client('s3',
        aws_access_key_id=cumulus['accessKeyId'],
        aws_secret_access_key=cumulus['secretAccessKey'],
        aws_session_token=cumulus['sessionToken'],
        region_name=region_name)
    #-- return the AWS client for region
    return client

#-- PURPOSE: get a s3 bucket name from a presigned url
def s3_bucket(presigned_url):
    """
    Get a s3 bucket name from a presigned url

    Parameters
    ----------
    presigned_url: str
        s3 presigned url

    Returns
    -------
    bucket: str
        s3 bucket name
    """
    host = url_split(presigned_url)
    bucket = re.sub(r's3:\/\/', r'', host[0], re.IGNORECASE)
    return bucket

#-- PURPOSE: get a s3 bucket key from a presigned url
def s3_key(presigned_url):
    """
    Get a s3 bucket key from a presigned url

    Parameters
    ----------
    presigned_url: str
        s3 presigned url

    Returns
    -------
    key: str
        s3 bucket key for object
    """
    host = url_split(presigned_url)
    key = posixpath.join(*host[1:])
    return key

#-- PURPOSE: check that entered NASA Earthdata credentials are valid
def check_credentials(HOST='https://podaac-tools.jpl.nasa.gov/drive/files'):
    """
    Check that entered NASA Earthdata credentials are valid

    HOST: str
        full url to protected credential website
    """
    try:
        request = urllib2.Request(HOST)
        response = urllib2.urlopen(request, timeout=20)
    except urllib2.HTTPError:
        raise RuntimeError('Check your NASA Earthdata credentials')
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: list a directory on JPL PO.DAAC/ECCO Drive https server
def drive_list(HOST, username=None, password=None, build=True,
    timeout=None, urs='podaac-tools.jpl.nasa.gov',
    parser=lxml.etree.HTMLParser(), pattern='', sort=False):
    """
    List a directory on
    `JPL PO.DAAC <https://podaac-tools.jpl.nasa.gov/drive>`_ or
    `ECCO Drive <https://ecco.jpl.nasa.gov/drive/>`_

    Parameters
    ----------
    HOST: str or list
        remote https host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        JPL PO.DAAC Drive WebDAV password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'podaac-tools.jpl.nasa.gov'
        JPL PO.DAAC or ECCO login URS 3 host
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for lxml
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    #-- use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try listing from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        tree = lxml.etree.parse(urllib2.urlopen(request,timeout=timeout),parser)
    except (urllib2.HTTPError, urllib2.URLError) as e:
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

#-- PURPOSE: download a file from a PO.DAAC/ECCO Drive https server
def from_drive(HOST, username=None, password=None, build=True,
    timeout=None, urs='podaac-tools.jpl.nasa.gov', local=None,
    hash='', chunk=16384, verbose=False, fid=sys.stdout, mode=0o775):
    """
    Download a file from a
    `JPL PO.DAAC <https://podaac-tools.jpl.nasa.gov/drive>`_ or
    `ECCO Drive <https://ecco.jpl.nasa.gov/drive/>`_ https server

    Parameters
    ----------
    HOST: str or list
        remote https host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        JPL PO.DAAC Drive WebDAV password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'podaac-tools.jpl.nasa.gov'
        JPL PO.DAAC or ECCO login URS 3 host
    local: str or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try downloading from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout)
    except (urllib2.HTTPError, urllib2.URLError) as e:
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
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args))
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: retrieve shortnames for GRACE/GRACE-FO products
def cmr_product_shortname(mission, center, release, level='L2', version='0'):
    """
    Create a list of product shortnames for CMR queries

    Parameters
    ----------
    mission: str
        GRACE (grace) or GRACE Follow-On (grace-fo)
    center: str
        GRACE/GRACE-FO processing center
    release: str
        GRACE/GRACE-FO data release
    level: str, default 'L2'
        GRACE/GRACE-FO product level

            - ``'L1A'``
            - ``'L1B'``
            - ``'L2'``
    version: str, default '0'
        GRACE/GRACE-FO Level-2 data version

    Returns
    -------
    cmr_shortnames: list
        shortnames for CMR queries
    """
    #-- build dictionary for GRACE/GRACE-FO shortnames
    cmr_shortname = {}
    cmr_shortname['grace'] = {}
    cmr_shortname['grace-fo'] = {}
    #-- format of GRACE/GRACE-FO shortnames
    grace_l1_format = 'GRACE_{0}_GRAV_{1}_{2}'
    grace_l2_format = 'GRACE_{0}_{1}_GRAV_{2}_{3}'
    gracefo_l1_format = 'GRACEFO_{0}_{1}_GRAV_{2}_{3}'
    gracefo_l2_format = 'GRACEFO_{0}_{1}_MONTHLY_{2}{3}'
    #-- dictionary entries for each product level
    cmr_shortname['grace']['L1B'] = dict(GFZ={},JPL={})
    cmr_shortname['grace']['L2'] = dict(CSR={},GFZ={},JPL={})
    cmr_shortname['grace-fo']['L1A'] = dict(JPL={})
    cmr_shortname['grace-fo']['L1B'] = dict(JPL={})
    cmr_shortname['grace-fo']['L2'] = dict(CSR={},GFZ={},JPL={})

    #-- dictionary entry for GRACE Level-1B deliasing products
    #-- for each data release
    for rl in ['RL06']:
        shortname = grace_l1_format.format('AOD1B','GFZ',rl)
        cmr_shortname['grace']['L1B']['GFZ'][rl] = [shortname]

    #-- dictionary entries for GRACE Level-1B ranging data products
    #-- for each data release
    for rl in ['RL02','RL03']:
        shortname = grace_l1_format.format('L1B','JPL',rl)
        cmr_shortname['grace']['L1B']['JPL'][rl] = [shortname]

    #-- dictionary entries for GRACE Level-2 products
    #-- for each data release
    for rl in ['RL06']:
        #-- Center for Space Research (CSR)
        cmr_shortname['grace']['L2']['CSR'][rl] = []
        #-- German Research Centre for Geosciences (GFZ)
        cmr_shortname['grace']['L2']['GFZ'][rl] = []
        #-- NASA Jet Propulsion Laboratory (JPL)
        cmr_shortname['grace']['L2']['JPL'][rl] = []
        #-- create list of product shortnames for GRACE level-2 products
        #-- for each L2 data processing center
        for c in ['CSR','GFZ','JPL']:
            #-- for each level-2 product
            for p in ['GAA', 'GAB', 'GAC', 'GAD', 'GSM']:
                #-- skip atmospheric and oceanic dealiasing products for CSR
                if (c == 'CSR') and p in ('GAA', 'GAB'):
                    continue
                #-- shortname for center and product
                shortname = grace_l2_format.format(p,'L2',c,rl)
                cmr_shortname['grace']['L2'][c][rl].append(shortname)

    #-- dictionary entries for GRACE-FO Level-1 ranging data products
    #-- for each data release
    for rl in ['RL04']:
        for l in ['L1A','L1B']:
            shortname = gracefo_l1_format.format(l,'ASCII','JPL',rl)
            cmr_shortname['grace-fo'][l]['JPL'][rl] = [shortname]

    #-- dictionary entries for GRACE-FO Level-2 products
    #-- for each data release
    for rl in ['RL06']:
        rs = re.findall('\d+',rl).pop().zfill(3)
        for c in ['CSR','GFZ','JPL']:
            shortname = gracefo_l2_format.format('L2',c,rs,version)
            cmr_shortname['grace-fo']['L2'][c][rl] = [shortname]

    #-- try to retrieve the shortname for a given mission
    try:
        cmr_shortnames = cmr_shortname[mission][level][center][release]
    except Exception as e:
        raise Exception('NASA CMR shortname not found')
    else:
        return cmr_shortnames

def cmr_readable_granules(product, level='L2', solution='BA01', version='0'):
    """
    Create readable granule names pattern for CMR queries

    Parameters
    ----------
    product: str
        GRACE/GRACE-FO data product
    level: str, default 'L2'
        GRACE/GRACE-FO product level

            - ``'L1A'``
            - ``'L1B'``
            - ``'L2'``
    solution: str, default 'BA01'
        monthly gravity field solution for Release-06

            - ``'BA01'``: unconstrained monthly gravity field solution to d/o 60
            - ``'BB01'``: unconstrained monthly gravity field solution to d/o 96
            - ``'BC01'``: computed monthly dealiasing solution to d/o 180
    version: str, default '0'
        GRACE/GRACE-FO Level-2 data version

    Returns
    -------
    pattern: str
        readable granule names pattern for CMR queries
    """
    if (level == 'L1B') and (product == 'AOD1B'):
        pattern = 'AOD1B_*'
    elif (level == 'L1A') or (level == 'L1B'):
        pattern = 'grace*'
    elif (level == 'L2') and (product == 'GSM'):
        args = (product, solution, version)
        pattern = '{0}-2_???????-???????_????_?????_{1}_???{2}*'.format(*args)
    elif (level == 'L2'):
        args = (product, 'BC01', version)
        pattern = '{0}-2_???????-???????_????_?????_{1}_???{2}*'.format(*args)
    else:
        pattern = '*'
    #-- return readable granules pattern
    return pattern

#-- PURPOSE: filter the CMR json response for desired data files
def cmr_filter_json(search_results, endpoint="data"):
    """
    Filter the CMR json response for desired data files

    Parameters
    ----------
    search_results: dict
        json response from CMR query
    endpoint: str, default 'data'
        url endpoint type

            - ``'data'``: PO.DAAC https archive
            - ``'s3'``: PO.DAAC Cumulus AWS S3 bucket

    Returns
    -------
    granule_names: list
        GRACE/GRACE-FO granule names
    granule_urls: list
        GRACE/GRACE-FO granule urls
    granule_mtimes: list
        GRACE/GRACE-FO granule modification times
    """
    #-- output list of granule ids, urls and modified times
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    #-- check that there are urls for request
    if ('feed' not in search_results) or ('entry' not in search_results['feed']):
        return (granule_names,granule_urls)
    # descriptor links for each endpoint
    rel = {}
    rel['data'] = "http://esipfed.org/ns/fedsearch/1.1/data#"
    rel['s3'] = "http://esipfed.org/ns/fedsearch/1.1/s3#"
    #-- iterate over references and get cmr location
    for entry in search_results['feed']['entry']:
        granule_names.append(entry['title'])
        granule_mtimes.append(get_unix_time(entry['updated'],
            format='%Y-%m-%dT%H:%M:%S.%f%z'))
        for link in entry['links']:
            if (link['rel'] == rel[endpoint]):
                granule_urls.append(link['href'])
                break
    #-- return the list of urls, granule ids and modified times
    return (granule_names,granule_urls,granule_mtimes)

#-- PURPOSE: cmr queries for GRACE/GRACE-FO products
def cmr(mission=None, center=None, release=None, level='L2', product=None,
    solution='BA01', version='0', start_date=None, end_date=None,
    provider='POCLOUD', endpoint='data', verbose=False, fid=sys.stdout):
    """
    Query the NASA Common Metadata Repository (CMR) for GRACE/GRACE-FO data

    Parameters
    ----------
    mission: str or NoneType, default None
        GRACE (``'grace'``) or GRACE Follow-On (``'grace-fo'``)
    center: str or NoneType, default None
        GRACE/GRACE-FO processing center
    release: str or NoneType, default None
        GRACE/GRACE-FO data release
    level: str, default 'L2'
        GRACE/GRACE-FO product level
    product: str or NoneType, default None
        GRACE/GRACE-FO data product
    solution: str, default 'BA01'
        monthly gravity field solution for Release-06
    version: str, default '0'
        GRACE/GRACE-FO Level-2 data version
    start_date: str or NoneType, default None
        starting date for CMR product query
    end_date: str or NoneType, default None
        ending date for CMR product query
    provider: str, default 'POCLOUD'
        CMR data provider

            - ``'PODAAC'``: PO.DAAC Drive
            - ``'POCLOUD'``: PO.DAAC Cumulus
    endpoint: str, default 'data'
        url endpoint type

            - ``'data'``: PO.DAAC https archive
            - ``'s3'``: PO.DAAC Cumulus AWS S3 bucket
    verbose: bool, default False
        print CMR query information
    fid: obj, default sys.stdout
        open file object to print if verbose

    Returns
    -------
    granule_names: list
        GRACE/GRACE-FO granule names
    granule_urls: list
        GRACE/GRACE-FO granule urls
    granule_mtimes: list
        GRACE/GRACE-FO granule modification times
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- build urllib2 opener with SSL context
    #-- https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    #-- Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    handler.append(urllib2.HTTPSHandler(context=ssl.SSLContext()))
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    #-- build CMR query
    cmr_format = 'json'
    cmr_page_size = 2000
    CMR_HOST = ['https://cmr.earthdata.nasa.gov','search',
        'granules.{0}'.format(cmr_format)]
    #-- build list of CMR query parameters
    CMR_KEYS = []
    CMR_KEYS.append('?provider={0}'.format(provider))
    CMR_KEYS.append('&sort_key[]=start_date')
    CMR_KEYS.append('&sort_key[]=producer_granule_id')
    CMR_KEYS.append('&scroll=true')
    CMR_KEYS.append('&page_size={0}'.format(cmr_page_size))
    #-- dictionary of product shortnames
    short_names = cmr_product_shortname(mission, center, release,
        level=level, version=version)
    for short_name in short_names:
        CMR_KEYS.append('&short_name={0}'.format(short_name))
    #-- append keys for start and end time
    #-- verify that start and end times are in ISO format
    start_date = isoformat(start_date) if start_date else ''
    end_date = isoformat(end_date) if end_date else ''
    CMR_KEYS.append('&temporal={0},{1}'.format(start_date, end_date))
    #-- append keys for querying specific products
    CMR_KEYS.append("&options[readable_granule_name][pattern]=true")
    CMR_KEYS.append("&options[spatial][or]=true")
    readable_granule = cmr_readable_granules(product,
        level=level, solution=solution, version=version)
    CMR_KEYS.append("&readable_granule_name[]={0}".format(readable_granule))
    #-- full CMR query url
    cmr_query_url = "".join([posixpath.join(*CMR_HOST),*CMR_KEYS])
    logging.info('CMR request={0}'.format(cmr_query_url))
    #-- output list of granule names and urls
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    cmr_scroll_id = None
    while True:
        req = urllib2.Request(cmr_query_url)
        if cmr_scroll_id:
            req.add_header('cmr-scroll-id', cmr_scroll_id)
        response = opener.open(req)
        #-- get scroll id for next iteration
        if not cmr_scroll_id:
            headers = {k.lower():v for k,v in dict(response.info()).items()}
            cmr_scroll_id = headers['cmr-scroll-id']
        #-- read the CMR search as JSON
        search_page = json.loads(response.read().decode('utf8'))
        ids,urls,mtimes = cmr_filter_json(search_page, endpoint=endpoint)
        if not urls:
            break
        #-- extend lists
        granule_names.extend(ids)
        granule_urls.extend(urls)
        granule_mtimes.extend(mtimes)
    #-- return the list of granule ids, urls and modification times
    return (granule_names, granule_urls, granule_mtimes)

#-- PURPOSE: download geocenter files from Sutterley and Velicogna (2019)
#-- https://doi.org/10.3390/rs11182108
#-- https://doi.org/10.6084/m9.figshare.7388540
def from_figshare(directory,article='7388540',timeout=None,
    context=ssl.SSLContext(),chunk=16384,verbose=False,fid=sys.stdout,
    pattern=r'(CSR|GFZ|JPL)_(RL\d+)_(.*?)_SLF_iter.txt$',mode=0o775):
    """
    Download [Sutterley2019]_ geocenter files from
    `figshare <https://doi.org/10.6084/m9.figshare.7388540>`_

    Parameters
    ----------
    directory: str
        download directory
    article: str
        figshare article number
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    pattern: str, default '(CSR|GFZ|JPL)_(RL\d+)_(.*?)_SLF_iter.txt$'
        regular expression pattern for reducing list
    mode: oct, default 0o775
        permissions mode of output local file

    References
    ----------
    .. [Sutterley2019] T. C. Sutterley and I. Velicogna,
        "Improved Estimates of Geocenter Variability from Time-Variable Gravity
        and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019).
        `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_
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

#-- PURPOSE: send files to figshare using secure FTP uploader
def to_figshare(files,username=None,password=None,directory=None,
    timeout=None,context=ssl.SSLContext(ssl.PROTOCOL_TLS),
    get_ca_certs=False,verbose=False,chunk=8192):
    """
    Send files to figshare using secure `FTP uploader
    <https://help.figshare.com/article/upload-large-datasets-and-bulk-upload-using-the-ftp-uploader-desktop-uploader-or-api>`_

    Parameters
    ----------
    files: list
        files to upload
    username: str or NoneType, default None
        ftp username
    password: str or NoneType, default None
        ftp password
    directory: str or NoneType, default None
        figshare subdirectory for sending data
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for url opener object
    get_ca_certs: bool, default False
        get list of loaded “certification authority” certificates
    verbose: bool, default False
        print ftp transfer information
    chunk: int, default 8192
        chunk size for transfer encoding
    """
    #-- SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    #-- connect to figshare secure ftp host
    ftps = ftplib.FTP_TLS(host='ftps.figshare.com',
        user=username,
        passwd=password,
        context=context,
        timeout=timeout)
    #-- set the verbosity level
    ftps.set_debuglevel(1) if verbose else None
    #-- encrypt data connections
    ftps.prot_p()
    #-- try to create project directory
    try:
        #-- will only create the directory if non-existent
        ftps.mkd(posixpath.join('data',directory))
    except:
        pass
    #-- upload each file
    for local_file in files:
        #-- remote ftp file
        ftp_remote_path = posixpath.join('data',directory,
            os.path.basename(local_file))
        #-- open local file and send bytes
        with open(os.path.expanduser(local_file),'rb') as fp:
            ftps.storbinary('STOR {0}'.format(ftp_remote_path), fp,
                blocksize=chunk, callback=None, rest=None)

#-- PURPOSE: download satellite laser ranging files from CSR
#-- http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
#-- http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s
def from_csr(directory,timeout=None,context=ssl.SSLContext(),
    chunk=16384,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download `satellite laser ranging (SLR) <http://download.csr.utexas.edu/pub/slr/>`_
    files from the University of Texas Center for Space Research (UTCSR)

    Parameters
    ----------
    directory: str
        download directory
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext()
        SSL context for url opener object
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default fid.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file
    """
    #-- CSR download http server
    HOST = 'http://download.csr.utexas.edu'
    #-- recursively create directory if non-existent
    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.access(os.path.join(directory,'geocenter'), os.F_OK):
        os.makedirs(os.path.join(directory,'geocenter'), mode)
    #-- download SLR 5x5, figure axis and azimuthal dependence files
    FILES = []
    FILES.append([HOST,'pub','slr','degree_5',
        'CSR_Monthly_5x5_Gravity_Harmonics.txt'])
    FILES.append([HOST,'pub','slr','degree_2','C20_RL06.txt'])
    FILES.append([HOST,'pub','slr','degree_2','C21_S21_RL06.txt'])
    FILES.append([HOST,'pub','slr','degree_2','C22_S22_RL06.txt'])
    #-- for each SLR file
    for FILE in FILES:
        original_md5 = get_hash(os.path.join(directory,FILE[-1]))
        from_http(FILE,timeout=timeout,context=context,
            local=os.path.join(directory,FILE[-1]),
            hash=original_md5,chunk=chunk,verbose=verbose,
            fid=fid,mode=mode)
    #-- download CF-CM SLR and updated SLR geocenter files from Minkang Cheng
    FILES = []
    FILES.append([HOST,'pub','slr','geocenter','GCN_L1_L2_30d_CF-CM.txt'])
    FILES.append([HOST,'outgoing','cheng','gct2est.220_5s'])
    #-- for each SLR geocenter file
    for FILE in FILES:
        original_md5 = get_hash(os.path.join(directory,'geocenter',FILE[-1]))
        from_http(FILE,timeout=timeout,context=context,
            local=os.path.join(directory,'geocenter',FILE[-1]),
            hash=original_md5,chunk=chunk,verbose=verbose,
            fid=fid,mode=mode)

#-- PURPOSE: download GravIS and satellite laser ranging files from GFZ
#-- ftp://isdcftp.gfz-potsdam.de/grace/Level-2/GFZ/RL06_SLR_C20/
#-- ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
def from_gfz(directory,timeout=None,chunk=8192,verbose=False,fid=sys.stdout,
    mode=0o775):
    """
    Download GravIS and satellite laser ranging (SLR) files from the
        German Research Centre for Geosciences (GeoForschungsZentrum, GFZ)

    Parameters
    ----------
    directory: str
        download directory
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    chunk: int, default 8192
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file
    """
    #-- recursively create directories if non-existent
    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.access(os.path.join(directory,'geocenter'), os.F_OK):
        os.makedirs(os.path.join(directory,'geocenter'), mode)
    #-- SLR oblateness and combined low-degree harmonic files
    FILES = []
    FILES.append(['isdcftp.gfz-potsdam.de','grace','Level-2','GFZ',
        'RL06_SLR_C20','GFZ_RL06_C20_SLR.dat'])
    FILES.append(['isdcftp.gfz-potsdam.de','grace','GravIS','GFZ',
        'Level-2B','aux_data','GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat'])
    #-- get each file
    for FILE in FILES:
        local = os.path.join(directory,FILE[-1])
        from_ftp(FILE,timeout=timeout,local=local,hash=get_hash(local),
            chunk=chunk,verbose=verbose,fid=fid,mode=mode)
    #-- GravIS geocenter file
    FILE = ['isdcftp.gfz-potsdam.de','grace','GravIS','GFZ','Level-2B',
        'aux_data','GRAVIS-2B_GFZOP_GEOCENTER_0002.dat']
    local = os.path.join(directory,'geocenter',FILE[-1])
    from_ftp(FILE,timeout=timeout,local=local,hash=get_hash(local),
        chunk=chunk,verbose=verbose,fid=fid,mode=mode)

#-- PURPOSE: list a directory on the GFZ ICGEM https server
#-- http://icgem.gfz-potsdam.de
def icgem_list(host='http://icgem.gfz-potsdam.de/tom_longtime',timeout=None,
    parser=lxml.etree.HTMLParser()):
    """
    Parse the table of static gravity field models on the GFZ
    `International Centre for Global Earth Models (ICGEM) <http://icgem.gfz-potsdam.de/>`_
    server

    Parameters
    ----------
    host: str
        url for the GFZ ICGEM gravity field table
    timeout: int or NoneType
        timeout in seconds for blocking operations
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for lxml

    Returns
    -------
    colfiles: dict
        Static gravity field file urls mapped by field name
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
