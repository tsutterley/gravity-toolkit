#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (11/2024)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

UPDATE HISTORY:
    Updated 11/2024: simplify unique file name function
        add function to scrape GSFC website for GRACE mascon urls
    Updated 10/2024: update CMR search utility to replace deprecated scrolling
        https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html
    Updated 08/2024: generalize hash function to use any available algorithm
    Updated 06/2024: added wrapper to importlib for optional dependencies
        make default case for an import exception be a class
    Updated 04/2024: added argument for products in CMR shortname query
    Updated 11/2023: updated ssl context to fix deprecation error
    Updated 10/2023: add capability to download CSR LRI solutions
    Updated 06/2023: add functions to retrieve and revoke Earthdata tokens
        add TN11e.txt file to list of CSR SLR downloads
    Updated 05/2023: add reify decorator for evaluation of properties
        use pathlib to define and operate on paths
    Updated 04/2023: use release-03 GFZ GravIS SLR and geocenter files
    Updated 03/2023: place boto3 import within try/except statement
    Updated 01/2023: add default ssl context attribute with protocol
    Updated 12/2022: add variables for NASA DAAC and s3 providers
        add functions for managing and maintaining git repositories
    Updated 11/2022: add CMR queries for collection metadata
        exposed GSFC SLR url for weekly 5x5 harmonics as an option
    Updated 08/2022: add regular expression function for finding files
    Updated 07/2022: add s3 endpoints and buckets for Earthdata Cumulus
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
from __future__ import print_function, division, annotations

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
import getpass
import inspect
import hashlib
import logging
import pathlib
import builtins
import dateutil
import warnings
import importlib
import posixpath
import lxml.etree
import subprocess
import calendar,time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
    import urllib.request as urllib2

# PURPOSE: get absolute path within a package from a relative path
def get_data_path(relpath: list | str | pathlib.Path):
    """
    Get the absolute path within a package from a relative path

    Parameters
    ----------
    relpath: list, str or pathlib.Path
        relative path
    """
    # current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = pathlib.Path(filename).absolute().parent
    if isinstance(relpath, list):
        # use *splat operator to extract from list
        return filepath.joinpath(*relpath)
    elif isinstance(relpath, str):
        return filepath.joinpath(relpath)


def import_dependency(
        name: str,
        extra: str = "",
        raise_exception: bool = False
    ):
    """
    Import an optional dependency

    Adapted from ``pandas.compat._optional::import_optional_dependency``

    Parameters
    ----------
    name: str
        Module name
    extra: str, default ""
        Additional text to include in the ``ImportError`` message
    raise_exception: bool, default False
        Raise an ``ImportError`` if the module is not found

    Returns
    -------
    module: obj
        Imported module
    """
    # check if the module name is a string
    msg = f"Invalid module name: '{name}'; must be a string"
    assert isinstance(name, str), msg
    # default error if module cannot be imported
    err = f"Missing optional dependency '{name}'. {extra}"
    module = type('module', (), {})
    # try to import the module
    try:
        module = importlib.import_module(name)
    except (ImportError, ModuleNotFoundError) as exc:
        if raise_exception:
            raise ImportError(err) from exc
        else:
            logging.debug(err)
    # return the module
    return module

class reify(object):
    """Class decorator that puts the result of the method it
    decorates into the instance"""
    def __init__(self, wrapped):
        self.wrapped = wrapped
        self.__name__ = wrapped.__name__
        self.__doc__ = wrapped.__doc__

    def __get__(self, inst, objtype=None):
        if inst is None:
            return self
        val = self.wrapped(inst)
        setattr(inst, self.wrapped.__name__, val)
        return val

# PURPOSE: get the hash value of a file
def get_hash(
        local: str | io.IOBase | pathlib.Path,
        algorithm: str = 'md5'
    ):
    """
    Get the hash value from a local file or ``BytesIO`` object

    Parameters
    ----------
    local: obj, str or pathlib.Path
        BytesIO object or path to file
    algorithm: str, default 'md5'
        hashing algorithm for checksum validation
    """
    # check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        # generate checksum hash for a given type
        if algorithm in hashlib.algorithms_available:
            return hashlib.new(algorithm, local.getvalue()).hexdigest()
        else:
            raise ValueError(f'Invalid hashing algorithm: {algorithm}')
    elif isinstance(local, (str, pathlib.Path)):
        # generate checksum hash for local file
        local = pathlib.Path(local).expanduser()
        # if file currently doesn't exist, return empty string
        if not local.exists():
            return ''
        # open the local_file in binary read mode
        with local.open(mode='rb') as local_buffer:
            # generate checksum hash for a given type
            if algorithm in hashlib.algorithms_available:
                return hashlib.new(algorithm, local_buffer.read()).hexdigest()
            else:
                raise ValueError(f'Invalid hashing algorithm: {algorithm}')
    else:
        return ''

# PURPOSE: get the git hash value
def get_git_revision_hash(
        refname: str = 'HEAD',
        short: bool = False
    ):
    """
    Get the ``git`` hash value for a particular reference

    Parameters
    ----------
    refname: str, default HEAD
        Symbolic reference name
    short: bool, default False
        Return the shorted hash value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'rev-parse']
    cmd.append('--short') if short else None
    cmd.append(refname)
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

# PURPOSE: get the current git status
def get_git_status():
    """Get the status of a ``git`` repository as a boolean value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'status', '--porcelain']
    with warnings.catch_warnings():
        return bool(subprocess.check_output(cmd))

# PURPOSE: recursively split a url path
def url_split(s: str):
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

# PURPOSE: convert file lines to arguments
def convert_arg_line_to_args(arg_line):
    """
    Convert file lines to arguments

    Parameters
    ----------
    arg_line: str
        line string containing a single argument and/or comments
    """
    # remove commented lines and after argument comments
    for arg in re.sub(r'\#(.*?)$',r'',arg_line).split():
        if not arg.strip():
            continue
        yield arg

# PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(
        time_string: str,
        format: str = '%Y-%m-%d %H:%M:%S'
    ):
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
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.timestamp()

# PURPOSE: output a time string in isoformat
def isoformat(time_string: str):
    """
    Reformat a date string to ISO formatting

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    """
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.isoformat()

# PURPOSE: rounds a number to an even number less than or equal to original
def even(value: float):
    """
    Rounds a number to an even number less than or equal to original

    Parameters
    ----------
    value: float
        number to be rounded
    """
    return 2*int(value//2)

# PURPOSE: rounds a number upward to its nearest integer
def ceil(value: float):
    """
    Rounds a number upward to its nearest integer

    Parameters
    ----------
    value: float
        number to be rounded upward
    """
    return -int(-value//1)

# PURPOSE: make a copy of a file with all system information
def copy(
        source: str | pathlib.Path,
        destination: str | pathlib.Path,
        move: bool = False,
        **kwargs
    ):
    """
    Copy or move a file with all system information

    Parameters
    ----------
    source: str or pathlib.Path
        source file
    destination: str or pathlib.Path
        copied destination file
    move: bool, default False
        remove the source file
    """
    source = pathlib.Path(source).expanduser().absolute()
    destination = pathlib.Path(destination).expanduser().absolute()
    # log source and destination
    logging.info(f'{str(source)} -->\n\t{str(destination)}')
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    # remove the original file if moving
    if move:
        source.unlink()

# PURPOSE: open a unique file adding a numerical instance if existing
def create_unique_file(filename: str | pathlib.Path):
    """
    Open a unique file adding a numerical instance if existing

    Parameters
    ----------
    filename: str or pathlib.Path
        full path to output file
    """
    # validate input filename
    filename = pathlib.Path(filename).expanduser().absolute()
    stem, suffix = filename.stem, filename.suffix
    # create counter to add to the end of the filename if existing
    counter = 1
    while counter:
        try:
            # open file descriptor only if the file doesn't exist
            fd = filename.open(mode='xb')
        except OSError:
            pass
        else:
            # return the file descriptor
            return fd
        # new filename adds a counter before the file extension
        filename = filename.with_name(f'{stem}_{counter:d}{suffix}')
        counter += 1

# PURPOSE: check ftp connection
def check_ftp_connection(
        HOST: str,
        username: str | None = None,
        password: str | None = None
    ):
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
    # attempt to connect to ftp host
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

# PURPOSE: list a directory on a ftp host
def ftp_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        basename: bool = False,
        pattern: str | None = None,
        sort: bool = False
    ):
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
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        # get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        # iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                # try sending modification time command
                mdtm = ftp.sendcmd(f'MDTM {f}')
            except ftplib.error_perm:
                # directories will return with an error
                pass
            else:
                # convert the modification time into unix time
                mtimes[i] = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            # reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            # sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # close the ftp connection
        ftp.close()
        # return the list of items and last modified times
        return (output, mtimes)

# PURPOSE: download a file from a ftp host
def from_ftp(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 8192,
        verbose: bool = False,
        fid=sys.stdout,
        mode: oct = 0o775
    ):
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
    local: str, pathlib.Path or NoneType, default None
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
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from ftp
    try:
        # try to connect to ftp host
        ftp = ftplib.FTP(HOST[0], timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary(f'RETR {ftp_remote_path}',
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd(f'MDTM {ftp_remote_path}')
        remote_mtime = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
            # keep remote modification time of file and local access time
            os.utime(local, (local.stat().st_atime, remote_mtime))
        # close the ftp connection
        ftp.close()
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

def _create_default_ssl_context() -> ssl.SSLContext:
    """Creates the default SSL context
    """
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    _set_ssl_context_options(context)
    context.options |= ssl.OP_NO_COMPRESSION
    return context

def _create_ssl_context_no_verify() -> ssl.SSLContext:
    """Creates an SSL context for unverified connections
    """
    context = _create_default_ssl_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    return context

def _set_ssl_context_options(context: ssl.SSLContext) -> None:
    """Sets the default options for the SSL context
    """
    if sys.version_info >= (3, 10) or ssl.OPENSSL_VERSION_INFO >= (1, 1, 0, 7):
        context.minimum_version = ssl.TLSVersion.TLSv1_2
    else:
        context.options |= ssl.OP_NO_SSLv2
        context.options |= ssl.OP_NO_SSLv3
        context.options |= ssl.OP_NO_TLSv1
        context.options |= ssl.OP_NO_TLSv1_1

# default ssl context
_default_ssl_context = _create_ssl_context_no_verify()

# PURPOSE: check internet connection
def check_connection(
        HOST: str,
        context: ssl.SSLContext = _default_ssl_context,
    ):
    """
    Check internet connection with http host

    Parameters
    ----------
    HOST: str
        remote http host
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    """
    # attempt to connect to http host
    try:
        urllib2.urlopen(HOST, timeout=20, context=context)
    except urllib2.URLError as exc:
        raise RuntimeError('Check internet connection') from exc
    else:
        return True

# PURPOSE: list a directory on an Apache http Server
def http_list(
        HOST: str | list,
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        parser = lxml.etree.HTMLParser(),
        format: str = '%Y-%m-%d %H:%M',
        pattern: str = '',
        sort: bool = False
    ):
    """
    List a directory on an Apache http Server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
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
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response, parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        # get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern, f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames, collastmod)

# PURPOSE: download a file from a http host
def from_http(
        HOST: str | list,
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from a http host

    Parameters
    ----------
    HOST: str or list
        remote http host path split as list
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    local: str, pathlib.Path or NoneType, default None
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
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except:
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: load a JSON response from a http host
def from_json(
        HOST: str | list,
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context
    ) -> dict:
    """
    Load a JSON response from a http host

    Parameters
    ----------
    HOST: str or list
        remote http host path split as list
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default pyTMD.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    """
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try loading JSON from http
    try:
        # Create and submit request for JSON response
        request = urllib2.Request(posixpath.join(*HOST))
        request.add_header('Accept', 'application/json')
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        msg = 'Load error from {0}'.format(posixpath.join(*HOST))
        raise Exception(msg) from exc
    else:
        # load JSON response
        return json.loads(response.read())

# PURPOSE: attempt to build an opener with netrc
def attempt_login(
        urs: str,
        context: ssl.SSLContext = _default_ssl_context,
        password_manager: bool = True,
        get_ca_certs: bool = False,
        redirect: bool = False,
        authorization_header: bool = True,
        **kwargs
    ):
    """
    Attempt to build a ``urllib`` opener for NASA Earthdata

    Parameters
    ----------
    urs: str
        Earthdata login URS 3 host
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
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
    kwargs.setdefault('netrc', pathlib.Path.home().joinpath('.netrc'))
    try:
        # verify permissions level of netrc file
        # only necessary on jupyterhub
        nc = pathlib.Path(kwargs['netrc']).expanduser().absolute()
        nc.chmod(mode=0o600)
        # try retrieving credentials from netrc
        username, _, password = netrc.netrc(nc).authenticators(urs)
    except Exception as exc:
        # try retrieving credentials from environmental variables
        username, password = (kwargs['username'], kwargs['password'])
        pass
    # if username or password are not available
    if not username:
        username = builtins.input(f'Username for {urs}: ')
    if not password:
        prompt = f'Password for {username}@{urs}: '
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
        except Exception as exc:
            pass
        else:
            return opener
        # reattempt login
        username = builtins.input(f'Username for {urs}: ')
        password = getpass.getpass(prompt=prompt)
    # reached end of available retries
    raise RuntimeError('End of Retries: Check NASA Earthdata credentials')

# PURPOSE: "login" to NASA Earthdata with supplied credentials
def build_opener(
        username: str,
        password: str,
        context: ssl.SSLContext = _default_ssl_context,
        password_manager: bool = False,
        get_ca_certs: bool = False,
        redirect: bool = False,
        authorization_header: bool = True,
        urs: str = 'https://urs.earthdata.nasa.gov'
    ):
    """
    Build ``urllib`` opener for NASA Earthdata with supplied credentials

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
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
    # https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    # create a password manager
    if password_manager:
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        # Add the username and password for NASA Earthdata Login system
        password_mgr.add_password(None, urs, username, password)
        handler.append(urllib2.HTTPBasicAuthHandler(password_mgr))
    # Create cookie jar for storing cookies. This is used to store and return
    # the session cookie given to use by the data server (otherwise will just
    # keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    # SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    handler.append(urllib2.HTTPSHandler(context=context))
    # redirect handler
    if redirect:
        handler.append(urllib2.HTTPRedirectHandler())
    # create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    # Encode username/password for request authorization headers
    # add Authorization header to opener
    if authorization_header:
        b64 = base64.b64encode(f'{username}:{password}'.encode())
        opener.addheaders = [("Authorization", f"Basic {b64.decode()}")]
    # Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    # All calls to urllib2.urlopen will now use handler
    # Make sure not to include the protocol in with the URL, or
    # HTTPPasswordMgrWithDefaultRealm will be confused.
    return opener

# PURPOSE: generate a NASA Earthdata user token
def get_token(
        HOST: str = 'https://urs.earthdata.nasa.gov/api/users/token',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        context: ssl.SSLContext = _default_ssl_context,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    Generate a NASA Earthdata User Token

    Parameters
    ----------
    HOST: str or list
        NASA Earthdata token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host

    Returns
    -------
    token: dict
        JSON response with NASA Earthdata User Token
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            context=context,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # create post response with Earthdata token API
    try:
        request = urllib2.Request(HOST, method='POST')
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # read and return JSON response
    return json.loads(response.read())

# PURPOSE: generate a NASA Earthdata user token
def list_tokens(
        HOST: str = 'https://urs.earthdata.nasa.gov/api/users/tokens',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        context: ssl.SSLContext = _default_ssl_context,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    List the current associated NASA Earthdata User Tokens

    Parameters
    ----------
    HOST: str
        NASA Earthdata list token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host

    Returns
    -------
    tokens: list
        JSON response with NASA Earthdata User Tokens
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            context=context,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # create get response with Earthdata list tokens API
    try:
        request = urllib2.Request(HOST)
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # read and return JSON response
    return json.loads(response.read())

# PURPOSE: revoke a NASA Earthdata user token
def revoke_token(
        token: str,
        HOST: str = f'https://urs.earthdata.nasa.gov/api/users/revoke_token',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        context: ssl.SSLContext = _default_ssl_context,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    Generate a NASA Earthdata User Token

    Parameters
    ----------
    token: str
        NASA Earthdata token to be revoked
    HOST: str
        NASA Earthdata revoke token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            context=context,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # full path for NASA Earthdata revoke token API
    url = f'{HOST}?token={token}'
    # create post response with Earthdata revoke tokens API
    try:
        request = urllib2.Request(url, method='POST')
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # verbose response
    logging.debug(f'Token Revoked: {token}')

# NASA on-prem DAAC providers
_daac_providers = {
    'gesdisc': 'GES_DISC',
    'ghrcdaac': 'GHRC_DAAC',
    'lpdaac': 'LPDAAC_ECS',
    'nsidc': 'NSIDC_ECS',
    'ornldaac': 'ORNL_DAAC',
    'podaac': 'PODAAC',
}

# NASA Cumulus AWS providers
_s3_providers = {
    'gesdisc': 'GES_DISC',
    'ghrcdaac': 'GHRC_DAAC',
    'lpdaac': 'LPCLOUD',
    'nsidc': 'NSIDC_CPRD',
    'ornldaac': 'ORNL_CLOUD',
    'podaac': 'POCLOUD',
}

# NASA Cumulus AWS S3 credential endpoints
_s3_endpoints = {
    'gesdisc': 'https://data.gesdisc.earthdata.nasa.gov/s3credentials',
    'ghrcdaac': 'https://data.ghrc.earthdata.nasa.gov/s3credentials',
    'lpdaac': 'https://data.lpdaac.earthdatacloud.nasa.gov/s3credentials',
    'nsidc': 'https://data.nsidc.earthdatacloud.nasa.gov/s3credentials',
    'ornldaac': 'https://data.ornldaac.earthdata.nasa.gov/s3credentials',
    'podaac': 'https://archive.podaac.earthdata.nasa.gov/s3credentials'
}

# NASA Cumulus AWS S3 buckets
_s3_buckets = {
    'gesdisc': 'gesdisc-cumulus-prod-protected',
    'ghrcdaac': 'ghrc-cumulus-dev',
    'lpdaac': 'lp-prod-protected',
    'nsidc': 'nsidc-cumulus-prod-protected',
    'ornldaac': 'ornl-cumulus-prod-protected',
    'podaac': 'podaac-ops-cumulus-protected',
    'podaac-doc': 'podaac-ops-cumulus-docs'
}

def s3_region():
    """
    Get AWS s3 region for EC2 instance

    Returns
    -------
    region_name: str
        AWS region name
    """
    boto3 = import_dependency('boto3')
    region_name = boto3.session.Session().region_name
    return region_name

# PURPOSE: get AWS s3 client for PO.DAAC Cumulus
def s3_client(
        HOST: str = _s3_endpoints['podaac'],
        timeout: int | None = None,
        region_name: str = 'us-west-2'
    ):
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
    # get AWS client object
    boto3 = import_dependency('boto3')
    client = boto3.client('s3',
        aws_access_key_id=cumulus['accessKeyId'],
        aws_secret_access_key=cumulus['secretAccessKey'],
        aws_session_token=cumulus['sessionToken'],
        region_name=region_name)
    # return the AWS client for region
    return client

# PURPOSE: get a s3 bucket name from a presigned url
def s3_bucket(presigned_url: str) -> str:
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

# PURPOSE: get a s3 bucket key from a presigned url
def s3_key(presigned_url: str) -> str:
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

# PURPOSE: check that entered NASA Earthdata credentials are valid
def check_credentials(HOST: str = _s3_endpoints['podaac']):
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

# PURPOSE: list a directory on JPL PO.DAAC/ECCO Drive https server
def drive_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        timeout: int | None = None,
        urs: str = 'podaac-tools.jpl.nasa.gov',
        parser = lxml.etree.HTMLParser(),
        pattern: str = '',
        sort: bool = False
    ):
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
        HTML parser for ``lxml``
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
    # use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    # build urllib2 opener and check credentials
    if build:
        # build urllib2 opener with credentials
        build_opener(username, password)
        # check credentials
        check_credentials()
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from https
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        tree = lxml.etree.parse(urllib2.urlopen(request, timeout=timeout),parser)
    except (urllib2.HTTPError, urllib2.URLError) as exc:
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for files (column names and modified times)
        colnames = tree.xpath('//tr/td//a[@class="text-left"]/text()')
        # get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i) for i in tree.xpath('//tr/td[3]/text()')]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames,collastmod)

# PURPOSE: download a file from a PO.DAAC/ECCO Drive https server
def from_drive(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        timeout: int | None = None,
        urs: str = 'podaac-tools.jpl.nasa.gov',
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
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
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    # build urllib2 opener and check credentials
    if build:
        # build urllib2 opener with credentials
        build_opener(username, password)
        # check credentials
        check_credentials()
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from https
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout)
    except (urllib2.HTTPError, urllib2.URLError) as exc:
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode=mode)
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: retrieve shortnames for GRACE/GRACE-FO products
def cmr_product_shortname(
        mission: str,
        center: str,
        release: str,
        level: str = 'L2',
        version: str = '0',
        product: list = ['GAA','GAB','GAC','GAD','GSM']
    ):
    """
    Create a list of product shortnames for NASA Common Metadata
    Repository (CMR) queries

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
    product: list, default ['GAA','GAB','GAC','GAD','GSM']
        GRACE/GRACE-FO Level-2 data products

    Returns
    -------
    cmr_shortnames: list
        shortnames for CMR queries
    """
    # build dictionary for GRACE/GRACE-FO shortnames
    cmr_shortname = {}
    cmr_shortname['grace'] = {}
    cmr_shortname['grace-fo'] = {}
    # format of GRACE/GRACE-FO shortnames
    grace_l1_format = 'GRACE_{0}_GRAV_{1}_{2}'
    grace_l2_format = 'GRACE_{0}_{1}_GRAV_{2}_{3}'
    gracefo_l1_format = 'GRACEFO_{0}_{1}_GRAV_{2}_{3}'
    gracefo_l2_format = 'GRACEFO_{0}_{1}_MONTHLY_{2}{3}'
    # dictionary entries for each product level
    cmr_shortname['grace']['L1B'] = dict(GFZ={},JPL={})
    cmr_shortname['grace']['L2'] = dict(CSR={},GFZ={},JPL={})
    cmr_shortname['grace-fo']['L1A'] = dict(JPL={})
    cmr_shortname['grace-fo']['L1B'] = dict(JPL={})
    cmr_shortname['grace-fo']['L2'] = dict(CSR={},GFZ={},JPL={})

    # dictionary entry for GRACE Level-1B dealiasing products
    # for each data release
    for rl in ['RL06']:
        shortname = grace_l1_format.format('AOD1B','GFZ',rl)
        cmr_shortname['grace']['L1B']['GFZ'][rl] = [shortname]

    # dictionary entries for GRACE Level-1B ranging data products
    # for each data release
    for rl in ['RL02','RL03']:
        shortname = grace_l1_format.format('L1B','JPL',rl)
        cmr_shortname['grace']['L1B']['JPL'][rl] = [shortname]

    # dictionary entries for GRACE Level-2 products
    # for each data release
    for rl in ['RL06']:
        # Center for Space Research (CSR)
        cmr_shortname['grace']['L2']['CSR'][rl] = []
        # German Research Centre for Geosciences (GFZ)
        cmr_shortname['grace']['L2']['GFZ'][rl] = []
        # NASA Jet Propulsion Laboratory (JPL)
        cmr_shortname['grace']['L2']['JPL'][rl] = []
        # check that product is iterable
        if isinstance(product, str):
            product = [product]
        # create list of product shortnames for GRACE level-2 products
        # for each L2 data processing center
        for c in ['CSR','GFZ','JPL']:
            # for each level-2 product
            for p in product:
                # skip atmospheric and oceanic dealiasing products for CSR
                if (c == 'CSR') and p in ('GAA', 'GAB'):
                    continue
                # shortname for center and product
                shortname = grace_l2_format.format(p,'L2',c,rl)
                cmr_shortname['grace']['L2'][c][rl].append(shortname)

    # dictionary entries for GRACE-FO Level-1 ranging data products
    # for each data release
    for rl in ['RL04']:
        for l in ['L1A','L1B']:
            shortname = gracefo_l1_format.format(l,'ASCII','JPL',rl)
            cmr_shortname['grace-fo'][l]['JPL'][rl] = [shortname]

    # dictionary entries for GRACE-FO Level-2 products
    # for each data release
    for rl in ['RL06']:
        rs = re.findall(r'\d+',rl).pop().zfill(3)
        for c in ['CSR','GFZ','JPL']:
            shortname = gracefo_l2_format.format('L2',c,rs,version)
            cmr_shortname['grace-fo']['L2'][c][rl] = [shortname]

    # try to retrieve the shortname for a given mission
    try:
        cmr_shortnames = cmr_shortname[mission][level][center][release]
    except Exception as exc:
        raise Exception('NASA CMR shortname not found')
    else:
        return cmr_shortnames

def cmr_readable_granules(
        product: str,
        level: str = 'L2',
        solution: str = 'BA01',
        version: str = '0'
    ):
    """
    Create readable granule names pattern for NASA Common Metadata
    Repository (CMR) queries

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
    # return readable granules pattern
    return pattern

# PURPOSE: filter the CMR json response for desired data files
def cmr_filter_json(
        search_results: dict,
        endpoint: str = 'data'
    ):
    """
    Filter the NASA Common Metadata Repository (CMR) json
    response for desired data files

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
    # output list of granule ids, urls and modified times
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    # check that there are urls for request
    if ('feed' not in search_results) or ('entry' not in search_results['feed']):
        return (granule_names,granule_urls)
    # descriptor links for each endpoint
    rel = {}
    rel['data'] = "http://esipfed.org/ns/fedsearch/1.1/data#"
    rel['s3'] = "http://esipfed.org/ns/fedsearch/1.1/s3#"
    # iterate over references and get cmr location
    for entry in search_results['feed']['entry']:
        granule_names.append(entry['title'])
        granule_mtimes.append(get_unix_time(entry['updated'],
            format='%Y-%m-%dT%H:%M:%S.%f%z'))
        for link in entry['links']:
            if (link['rel'] == rel[endpoint]):
                granule_urls.append(link['href'])
                break
    # return the list of urls, granule ids and modified times
    return (granule_names,granule_urls,granule_mtimes)

# PURPOSE: filter the CMR json response for desired metadata files
def cmr_metadata_json(
        search_results: dict,
        endpoint: str = 'data'
    ):
    """
    Filter the NASA Common Metadata Repository (CMR) json response
    for desired metadata files

    Parameters
    ----------
    search_results: dict
        json response from CMR query
    endpoint: str, default 'data'
        url endpoint type

            - ``'documentation'``: PO.DAAC documentation archive
            - ``'data'``: PO.DAAC https archive
            - ``'s3'``: PO.DAAC Cumulus AWS S3 bucket

    Returns
    -------
    collection_urls: list
        urls from collection of endpoint type
    """
    # output list of collection urls
    collection_urls = []
    # check that there are urls for request
    if ('feed' not in search_results) or ('entry' not in search_results['feed']):
        return collection_urls
    # descriptor links for each endpoint
    rel = {}
    rel['documentation'] = "http://esipfed.org/ns/fedsearch/1.1/documentation#"
    rel['data'] = "http://esipfed.org/ns/fedsearch/1.1/data#"
    rel['s3'] = "http://esipfed.org/ns/fedsearch/1.1/s3#"
    # iterate over references and get cmr location
    for entry in search_results['feed']['entry']:
        for link in entry['links']:
            if (link['rel'] == rel[endpoint]):
                collection_urls.append(link['href'])
    # return the list of urls
    return collection_urls

# PURPOSE: cmr queries for GRACE/GRACE-FO products
def cmr(
        mission: str | None = None,
        center: str | None = None,
        release: str | None = None,
        level: str | None = 'L2',
        product: str | None = None,
        solution: str | None = 'BA01',
        version: str | None = '0',
        start_date: str | None = None,
        end_date: str | None = None,
        provider: str | None = 'POCLOUD',
        endpoint: str | None = 'data',
        context: ssl.SSLContext = _default_ssl_context,
        verbose: bool = False,
        fid = sys.stdout
    ):
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
    level: str or NoneType, default 'L2'
        GRACE/GRACE-FO product level
    product: str or NoneType, default None
        GRACE/GRACE-FO data product
    solution: str or NoneType, default 'BA01'
        monthly gravity field solution for Release-06
    version: str or NoneType, default '0'
        GRACE/GRACE-FO Level-2 data version
    start_date: str or NoneType, default None
        starting date for CMR product query
    end_date: str or NoneType, default None
        ending date for CMR product query
    provider: str or NoneType, default 'POCLOUD'
        CMR data provider

            - ``'PODAAC'``: PO.DAAC Drive
            - ``'POCLOUD'``: PO.DAAC Cumulus
    endpoint: str or NoneType, default 'data'
        url endpoint type

            - ``'data'``: PO.DAAC https archive
            - ``'s3'``: PO.DAAC Cumulus AWS S3 bucket
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
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
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # build urllib2 opener with SSL context
    # https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    # Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    handler.append(urllib2.HTTPSHandler(context=context))
    # create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    # build CMR query
    cmr_query_type = 'granules'
    cmr_format = 'json'
    cmr_page_size = 2000
    CMR_HOST = ['https://cmr.earthdata.nasa.gov','search',
        f'{cmr_query_type}.{cmr_format}']
    # build list of CMR query parameters
    CMR_KEYS = []
    CMR_KEYS.append(f'?provider={provider}')
    CMR_KEYS.append('&sort_key[]=start_date')
    CMR_KEYS.append('&sort_key[]=producer_granule_id')
    CMR_KEYS.append(f'&page_size={cmr_page_size}')
    # dictionary of product shortnames
    short_names = cmr_product_shortname(mission, center, release,
        level=level, version=version)
    for short_name in short_names:
        CMR_KEYS.append(f'&short_name={short_name}')
    # append keys for start and end time
    # verify that start and end times are in ISO format
    start_date = isoformat(start_date) if start_date else ''
    end_date = isoformat(end_date) if end_date else ''
    CMR_KEYS.append(f'&temporal={start_date},{end_date}')
    # append keys for querying specific products
    CMR_KEYS.append("&options[readable_granule_name][pattern]=true")
    CMR_KEYS.append("&options[spatial][or]=true")
    readable_granule = cmr_readable_granules(product,
        level=level, solution=solution, version=version)
    CMR_KEYS.append(f"&readable_granule_name[]={readable_granule}")
    # full CMR query url
    cmr_query_url = "".join([posixpath.join(*CMR_HOST),*CMR_KEYS])
    logging.info(f'CMR request={cmr_query_url}')
    # output list of granule names and urls
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    cmr_search_after = None
    while True:
        req = urllib2.Request(cmr_query_url)
        # add CMR search after header
        if cmr_search_after:
            req.add_header('CMR-Search-After', cmr_search_after)
            logging.debug(f'CMR-Search-After: {cmr_search_after}')
        response = opener.open(req)
        # get search after index for next iteration
        headers = {k.lower():v for k,v in dict(response.info()).items()}
        cmr_search_after = headers.get('cmr-search-after')
        # read the CMR search as JSON
        search_page = json.loads(response.read().decode('utf8'))
        ids,urls,mtimes = cmr_filter_json(search_page, endpoint=endpoint)
        if not urls or cmr_search_after is None:
            break
        # extend lists
        granule_names.extend(ids)
        granule_urls.extend(urls)
        granule_mtimes.extend(mtimes)
    # return the list of granule ids, urls and modification times
    return (granule_names, granule_urls, granule_mtimes)

# PURPOSE: cmr queries for GRACE/GRACE-FO auxiliary data and documentation
def cmr_metadata(
        mission: str | None = None,
        center: str | None = None,
        release: str | None = None,
        level: str | None = 'L2',
        version: str | None = '0',
        provider: str | None = 'POCLOUD',
        endpoint: str | None = 'data',
        pattern: str | None = '',
        context: ssl.SSLContext = _default_ssl_context,
        verbose: bool = False,
        fid = sys.stdout
    ):
    """
    Query the NASA Common Metadata Repository (CMR) for GRACE/GRACE-FO
    auxiliary data and documentation

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
    version: str, default '0'
        GRACE/GRACE-FO Level-2 data version
    provider: str, default 'POCLOUD'
        CMR data provider

            - ``'PODAAC'``: PO.DAAC Drive
            - ``'POCLOUD'``: PO.DAAC Cumulus
    endpoint: str, default 'data'
        url endpoint type

            - ``'documentation'``: PO.DAAC documentation archive
            - ``'data'``: PO.DAAC https archive
            - ``'s3'``: PO.DAAC Cumulus AWS S3 bucket
    pattern: str, default ''
        regular expression pattern for reducing list
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    verbose: bool, default False
        print CMR query information
    fid: obj, default sys.stdout
        open file object to print if verbose

    Returns
    -------
    collection_urls: list
        urls from collection of endpoint type
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # build urllib2 opener with SSL context
    # https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    # Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    handler.append(urllib2.HTTPSHandler(context=context))
    # create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    # build CMR query
    cmr_query_type = 'collections'
    cmr_format = 'json'
    CMR_HOST = ['https://cmr.earthdata.nasa.gov','search',
        f'{cmr_query_type}.{cmr_format}']
    # build list of CMR query parameters
    CMR_KEYS = []
    CMR_KEYS.append(f'?provider={provider}')
    # dictionary of product shortnames
    short_names = cmr_product_shortname(mission, center, release,
        level=level, version=version)
    for short_name in short_names:
        CMR_KEYS.append(f'&short_name={short_name}')
    # full CMR query url
    cmr_query_url = "".join([posixpath.join(*CMR_HOST),*CMR_KEYS])
    logging.info(f'CMR request={cmr_query_url}')
    # query CMR for collection metadata
    req = urllib2.Request(cmr_query_url)
    response = opener.open(req)
    # read the CMR search as JSON
    search_page = json.loads(response.read().decode('utf8'))
    # filter the JSON response for desired endpoint links
    collection_urls = cmr_metadata_json(search_page, endpoint=endpoint)
    # reduce using regular expression pattern
    if pattern:
        i = [i for i,f in enumerate(collection_urls) if re.search(pattern,f)]
        # reduce list of collection_urls
        collection_urls = [collection_urls[indice] for indice in i]
    # return the list of collection urls
    return collection_urls

# PURPOSE: create and compile regular expression operator to find GRACE files
def compile_regex_pattern(
        PROC: str,
        DREL: str,
        DSET: str,
        mission: str | None = None,
        solution: str | None = r'BA01',
        version: str | None = r'\d+'
    ):
    """
    Compile regular expressor operators for finding a specified
    subset of GRACE/GRACE-FO Level-2 spherical harmonic files

    Parameters
    ----------
    PROC: str
        GRACE/GRACE-FO data processing center

            - ``'CNES'``: French Centre National D'Etudes Spatiales
            - ``'CSR'``: University of Texas Center for Space Research
            - ``'GFZ'``: German Research Centre for Geosciences (GeoForschungsZentrum)
            - ``'JPL'``: Jet Propulsion Laboratory
    DREL: str
        GRACE/GRACE-FO data release
    DSET: str
        GRACE/GRACE-FO data product

            - ``'GAA'``: non-tidal atmospheric correction
            - ``'GAB'``: non-tidal oceanic correction
            - ``'GAC'``: combined non-tidal atmospheric and oceanic correction
            - ``'GAD'``: ocean bottom pressure product
            - ``'GSM'``: corrected monthly static gravity field product
    mission: str or NoneType, default None
        GRACE/GRACE-FO mission shortname

            - ``'GRAC'``: GRACE
            - ``'GRFO'``: GRACE-FO
    solution: str, default 'BA01'
        monthly gravity field solution for Release-06

            - ``'BA01'``: unconstrained monthly gravity field solution to d/o 60
            - ``'BB01'``: unconstrained monthly gravity field solution to d/o 96
            - ``'BC01'``: computed monthly dealiasing solution to d/o 180
    version: str, default '0'
        GRACE/GRACE-FO Level-2 data version
    """
    # verify inputs
    if mission and mission not in ('GRAC','GRFO'):
        raise ValueError(f'Unknown mission {mission}')
    if PROC not in ('CNES','CSR','GFZ','JPL'):
        raise ValueError(f'Unknown processing center {PROC}')
    if DSET not in ('GAA','GAB','GAC','GAD','GSM'):
        raise ValueError(f'Unknown Level-2 product {DSET}')
    if isinstance(version, int):
        version = str(version).zfill(2)
    # compile regular expression operator for inputs
    if ((DSET == 'GSM') and (PROC == 'CSR') and (DREL in ('RL04','RL05'))):
        # CSR GSM: only monthly degree 60 products
        # not the longterm degree 180, degree 96 dataset or the
        # special order 30 datasets for the high-resonance months
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        pattern = r'{0}-2_\d+-\d+_\d+_UTCSR_0060_000{1:d}(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'CSR') and (DREL == 'RL06')):
        # CSR GSM RL06: monthly products for mission and solution
        release, = re.findall(r'\d+', DREL)
        args = (DSET, mission, solution, release.zfill(2), version.zfill(2))
        pattern = r'{0}-2_\d+-\d+_{1}_UTCSR_{2}_{3}{4}(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'CSR') and (DREL.endswith('LRI'))):
        # CSR GSM LRI solutions: monthly products for mission and solution
        release, version = re.findall(r'(\d+)\.(\d+)', DREL).pop()
        args = (DSET, mission, r'EA01', release.zfill(2), version.zfill(2))
        pattern = r'{0}-2_\d+-\d+_{1}_UTCSR_{2}_{3}{4}(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL04')):
        # GFZ RL04: only unconstrained solutions (not GK2 products)
        args = (DSET,)
        pattern = r'{0}-2_\d+-\d+_\d+_EIGEN_G---_0004(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL05')):
        # GFZ RL05: updated RL05a products which are less constrained to
        # the background model.  Allow regularized fields
        args = (DSET, r'(G---|GK2-)')
        pattern = r'{0}-2_\d+-\d+_\d+_EIGEN_{1}_005a(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'GFZ') and (DREL == 'RL06')):
        # GFZ GSM RL06: monthly products for mission and solution
        release, = re.findall(r'\d+', DREL)
        args = (DSET, mission, solution, release.zfill(2), version.zfill(2))
        pattern = r'{0}-2_\d+-\d+_{1}_GFZOP_{2}_{3}{4}(\.gz)?$'
    elif (PROC == 'JPL') and DREL in ('RL04','RL05'):
        # JPL: RL04a and RL05a products (denoted by 0001)
        release, = re.findall(r'\d+', DREL)
        args = (DSET, int(release))
        pattern = r'{0}-2_\d+-\d+_\d+_JPLEM_0001_000{1:d}(\.gz)?$'
    elif ((DSET == 'GSM') and (PROC == 'JPL') and (DREL == 'RL06')):
        # JPL GSM RL06: monthly products for mission and solution
        release, = re.findall(r'\d+', DREL)
        args = (DSET, mission, solution, release.zfill(2), version.zfill(2))
        pattern = r'{0}-2_\d+-\d+_{1}_JPLEM_{2}_{3}{4}(\.gz)?$'
    elif (PROC == 'CNES'):
        # CNES: use products in standard format
        args = (DSET,)
        pattern = r'{0}-2_\d+-\d+_\d+_GRGS_([a-zA-Z0-9_\-]+)(\.txt)?(\.gz)?$'
    elif mission is not None:
        # dealiasing products with mission listed
        args = (DSET, mission)
        pattern = r'{0}-2_([a-zA-Z0-9_\-]+)_{1}_([a-zA-Z0-9_\-]+)(\.gz)?$'
    else:
        # dealiasing products: use products in standard format
        args = (DSET,)
        pattern = r'{0}-2_([a-zA-Z0-9_\-]+)(\.gz)?$'
    # return the compiled regular expression operator
    return re.compile(pattern.format(*args), re.VERBOSE)

# PURPOSE: download geocenter files from Sutterley and Velicogna (2019)
# https://doi.org/10.3390/rs11182108
# https://doi.org/10.6084/m9.figshare.7388540
def from_figshare(
        directory: str | pathlib.Path,
        article: str = '7388540',
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        chunk: int | None = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        pattern: str = r'(CSR|GFZ|JPL)_(RL\d+)_(.*?)_SLF_iter.txt$',
        mode: oct = 0o775
    ):
    """
    Download :cite:p:`Sutterley:2019bx` geocenter files from
    `figshare <https://doi.org/10.6084/m9.figshare.7388540>`_

    Parameters
    ----------
    directory: str
        download directory
    article: str
        figshare article number
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    pattern: str, default r'(CSR|GFZ|JPL)_(RL\d+)_(.*?)_SLF_iter.txt$'
        regular expression pattern for reducing list
    mode: oct, default 0o775
        permissions mode of output local file
    """
    # figshare host
    HOST=['https://api.figshare.com','v2','articles',article]
    # recursively create directory if non-existent
    directory = pathlib.Path(directory).expanduser().absolute()
    local_dir = directory.joinpath('geocenter')
    local_dir.mkdir(mode=mode, parents=True, exist_ok=True)
    # Create and submit request.
    request = urllib2.Request(posixpath.join(*HOST))
    response = urllib2.urlopen(request, timeout=timeout,context=context)
    resp = json.loads(response.read())
    # reduce list of geocenter files
    geocenter_files = [f for f in resp['files'] if re.match(pattern,f['name'])]
    for f in geocenter_files:
        # download geocenter file
        local_file = local_dir.joinpath(f['name'])
        original_md5 = get_hash(local_file)
        from_http(f['download_url'],
            timeout=timeout,
            context=context,
            local=local_file,
            hash=original_md5,
            chunk=chunk,
            verbose=verbose,
            fid=fid,
            mode=mode)
        # verify MD5 checksums
        computed_md5 = get_hash(local_file)
        if (computed_md5 != f['supplied_md5']):
            raise Exception(f'Checksum mismatch: {f["download_url"]}')

# PURPOSE: send files to figshare using secure FTP uploader
def to_figshare(
        files: list,
        username: str | None = None,
        password: str | None = None,
        directory: str | None | pathlib.Path = None,
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        get_ca_certs: bool = False,
        verbose: bool = False,
        chunk: int = 8192
    ):
    """
    Send files to figshare using secure `FTP uploader
    <https://help.figshare.com/article/upload-large-datasets-and-
    bulk-upload-using-the-ftp-uploader-desktop-uploader-or-api>`_

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
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    get_ca_certs: bool, default False
        get list of loaded “certification authority” certificates
    verbose: bool, default False
        print ftp transfer information
    chunk: int, default 8192
        chunk size for transfer encoding
    """
    # SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    # connect to figshare secure ftp host
    ftps = ftplib.FTP_TLS(host='ftps.figshare.com',
        user=username,
        passwd=password,
        context=context,
        timeout=timeout)
    # set the verbosity level
    ftps.set_debuglevel(1) if verbose else None
    # encrypt data connections
    ftps.prot_p()
    # try to create project directory
    try:
        # will only create the directory if non-existent
        ftps.mkd(posixpath.join('data',directory))
    except:
        pass
    # upload each file
    for local_file in files:
        # local file
        local_file = pathlib.Path(local_file).expanduser().absolute()
        # remote ftp file
        ftp_remote_path = posixpath.join('data',directory,
            local_file.name)
        # open local file and send bytes
        with local_file.open(mode='rb') as fp:
            ftps.storbinary(f'STOR {ftp_remote_path}', fp,
                blocksize=chunk, callback=None, rest=None)

# PURPOSE: download files from CSR
# http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt
# http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s
def from_csr(
        directory: str | pathlib.Path,
        variable: str | list | tuple | None = None,
        version: str = 'RL06.1LRI',
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        chunk: int | None = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download files from the University of Texas Center for
    Space Research (UTCSR)

    Parameters
    ----------
    directory: str
        download directory
    variable: str, list, tuple or NoneType, default None
        CSR variable to download

            - ``'SLR'``: low degree SLR solutions
            - ``'geocenter'``: SLR geocenter solutions
            - ``'LRI'``: level-2 solutions from LRI
    version: str, default 'RL06.1LRI'
        Version of the LRI dataset to download
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default fid.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file
    """
    # CSR download http server
    HOST = 'http://download.csr.utexas.edu'
    # recursively create directory if non-existent
    directory = pathlib.Path(directory).expanduser().absolute()
    directory.mkdir(mode=mode, parents=True, exist_ok=True)
    # verify inputs for variable to be iterable
    if isinstance(variable, str):
        variable = [variable]
    # download SLR files from CSR
    if 'SLR' in variable:
        # download SLR 5x5, figure axis and azimuthal dependence files
        FILES = []
        FILES.append([HOST,'pub','slr','degree_5',
            'CSR_Monthly_5x5_Gravity_Harmonics.txt'])
        FILES.append([HOST,'pub','slr','degree_2','C20_RL06.txt'])
        FILES.append([HOST,'pub','slr','degree_2','C21_S21_RL06.txt'])
        FILES.append([HOST,'pub','slr','degree_2','C22_S22_RL06.txt'])
        FILES.append([HOST,'pub','slr','TN11E','TN11E.txt'])
        # for each SLR file
        for FILE in FILES:
            local_file = directory.joinpath(FILE[-1])
            original_md5 = get_hash(local_file)
            from_http(FILE,
                timeout=timeout,
                context=context,
                local=local_file,
                hash=original_md5,
                chunk=chunk,
                verbose=verbose,
                fid=fid,
                mode=mode)
    # download geocenter files from CSR
    if 'geocenter' in variable:
        # recursively create geocenter directory if non-existent
        local_dir = directory.joinpath('geocenter')
        local_dir.mkdir(mode=mode, parents=True, exist_ok=True)
        # download CF-CM SLR and updated SLR geocenter files from Minkang Cheng
        FILES = []
        FILES.append([HOST,'pub','slr','geocenter','GCN_L1_L2_30d_CF-CM.txt'])
        FILES.append([HOST,'outgoing','cheng','gct2est.220_5s'])
        # for each SLR geocenter file
        for FILE in FILES:
            local_file = local_dir.joinpath(FILE[-1])
            original_md5 = get_hash(local_file)
            from_http(FILE,
                timeout=timeout,
                context=context,
                local=local_file,
                hash=original_md5,
                chunk=chunk,
                verbose=verbose,
                fid=fid,
                mode=mode)
    # download LRI-only solutions
    if 'LRI' in variable:
        remote_path = ['http://download.csr.utexas.edu',
            'outgoing', 'gracefo', version]
        # find years of available LRI data
        years, _ = http_list(remote_path, pattern=r'\d{4}')
        # download each available CSR product
        for PROD in ['GAC','GAD','GSM']:
            # recursively create local directory if non-existent
            local_dir = directory.joinpath('CSR', version, PROD)
            local_dir.mkdir(mode=mode, parents=True, exist_ok=True)
            # for each year
            for year in years:
                # find LRI files
                files, mtimes = http_list([*remote_path, year], pattern=PROD)
                # download each file
                for fi, lmd in zip(files, mtimes):
                    local_file = local_dir.joinpath(fi)
                    original_md5 = get_hash(local_file)
                    from_http([*remote_path, year, fi],
                        timeout=timeout,
                        context=context,
                        local=local_file,
                        hash=original_md5,
                        chunk=chunk,
                        verbose=verbose,
                        fid=fid,
                        mode=mode)

# PURPOSE: download GravIS and satellite laser ranging files from GFZ
# ftp://isdcftp.gfz-potsdam.de/grace/Level-2/GFZ/RL06_SLR_C20/
# ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/
def from_gfz(
        directory: str | pathlib.Path,
        timeout: int | None = None,
        chunk: int | None = 8192,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
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
    # recursively create directories if non-existent
    directory = pathlib.Path(directory).expanduser().absolute()
    local_dir = directory.joinpath('geocenter')
    local_dir.mkdir(mode=mode, parents=True, exist_ok=True)
    # SLR oblateness and combined low-degree harmonic files
    FILES = []
    FILES.append(['isdcftp.gfz-potsdam.de','grace','Level-2','GFZ',
        'RL06_SLR_C20','GFZ_RL06_C20_SLR.dat'])
    FILES.append(['isdcftp.gfz-potsdam.de','grace','GravIS','GFZ',
        'Level-2B','aux_data','GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0003.dat'])
    # get each file
    for FILE in FILES:
        local_file = directory.joinpath(FILE[-1])
        from_ftp(FILE,
            timeout=timeout,
            local=local_file,
            hash=get_hash(local_file),
            chunk=chunk,
            verbose=verbose,
            fid=fid,
            mode=mode)
    # GravIS geocenter file
    FILE = ['isdcftp.gfz-potsdam.de','grace','GravIS','GFZ','Level-2B',
        'aux_data','GRAVIS-2B_GFZOP_GEOCENTER_0003.dat']
    local_file = local_dir.joinpath(FILE[-1])
    from_ftp(FILE,
        timeout=timeout,
        local=local_file,
        hash=get_hash(local_file),
        chunk=chunk,
        verbose=verbose,
        fid=fid,
        mode=mode)

# PURPOSE: lists files by scraping the GSFC grace-mascons website
def gsfc_list(
        HOST: str | list = 'https://earth.gsfc.nasa.gov/geo/data/grace-mascons',
        timeout: int | None = None,
        parser = lxml.etree.HTMLParser(),
        pattern: str = r'',
        sort: bool = False
    ):
    """
    Lists files by scraping the GSFC website for GRACE mascons

    Parameters
    ----------
    HOST: str or list
        remote https host
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    """
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from https
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        tree = lxml.etree.parse(urllib2.urlopen(request, timeout=timeout),parser)
    except (urllib2.HTTPError, urllib2.URLError) as exc:
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for relative links to files
        rellinks = tree.xpath('//tr/td//a/@href')
    # form complete column names
    colnames = [posixpath.join(HOST[0], *url_split(l)) for l in rellinks]
    # reduce using regular expression pattern
    if pattern:
        colnames = [f for i,f in enumerate(colnames) if re.search(pattern,f)]
    # sort list of column names
    if sort:
        colnames = [j for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
    # return the list of column names
    return colnames

# PURPOSE: download satellite laser ranging files from GSFC
# https://earth.gsfc.nasa.gov/geo/data/slr
def from_gsfc(
        directory: str | pathlib.Path,
        host: str = 'https://earth.gsfc.nasa.gov/sites/default/files/geo/slr-weekly',
        timeout: int | None = None,
        context: ssl.SSLContext = _default_ssl_context,
        chunk: int | None = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        copy: bool = False,
        mode: oct = 0o775
    ):
    """
    Download `satellite laser ranging (SLR) <https://earth.gsfc.nasa.gov/geo/data/slr/>`_
    files from NASA Goddard Space Flight Center (GSFC)

    Parameters
    ----------
    directory: str
        download directory
    host: str, default 'https://earth.gsfc.nasa.gov/sites/default/files/geo/slr-weekly'
        url for the GSFC SLR weekly fields
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default gravity_toolkit.utilities._default_ssl_context
        SSL context for ``urllib`` opener object
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default fid.stdout
        open file object to print if verbose
    copy: bool, default False
        create a copy of file for archival purposes
    mode: oct, default 0o775
        permissions mode of output local file
    """
    # recursively create directory if non-existent
    directory = pathlib.Path(directory).expanduser().absolute()
    directory.mkdir(mode=mode, parents=True, exist_ok=True)
    # download GSFC SLR 5x5 file
    FILE = 'gsfc_slr_5x5c61s61.txt'
    local_file = directory.joinpath(FILE)
    original_md5 = get_hash(local_file)
    fileID = from_http(posixpath.join(host,FILE),
        timeout=timeout,
        context=context,
        local=local_file,
        hash=original_md5,
        chunk=chunk,
        verbose=verbose,
        fid=fid,
        mode=mode)
    # create a dated copy for archival purposes
    if copy:
        # create copy of file for archiving
        # read file and extract data date span
        file_contents = fileID.read().decode('utf-8').splitlines()
        data_span, = [l for l in file_contents if l.startswith('Data span:')]
        # extract start and end of data date span
        span_start,span_end = re.findall(r'\d+[\s+]\w{3}[\s+]\d{4}', data_span)
        # create copy of file with date span in filename
        YM1 = time.strftime('%Y%m', time.strptime(span_start, '%d %b %Y'))
        YM2 = time.strftime('%Y%m', time.strptime(span_end, '%d %b %Y'))
        COPY = f'GSFC_SLR_5x5c61s61_{YM1}_{YM2}.txt'
        shutil.copyfile(local_file, directory.joinpath(COPY))
        # copy modification times and permissions for archive file
        shutil.copystat(local_file, directory.joinpath(COPY))

# PURPOSE: list a directory on the GFZ ICGEM https server
# http://icgem.gfz-potsdam.de
def icgem_list(
        host: str = 'http://icgem.gfz-potsdam.de/tom_longtime',
        timeout: int | None = None,
        parser=lxml.etree.HTMLParser()
    ):
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
        HTML parser for ``lxml``

    Returns
    -------
    colfiles: dict
        Static gravity field file urls mapped by field name
    """
    # try listing from https
    try:
        # Create and submit request.
        request = urllib2.Request(host)
        tree = lxml.etree.parse(urllib2.urlopen(request, timeout=timeout),parser)
    except:
        raise Exception(f'List error from {host}')
    else:
        # read and parse request for files
        colfiles = tree.xpath('//td[@class="tom-cell-modelfile"]//a/@href')
        # reduce list of files to find gfc files
        # return the dict of model files mapped by name
        return {re.findall(r'(.*?).gfc',posixpath.basename(f)).pop():url_split(f)
            for i,f in enumerate(colfiles) if re.search(r'gfc$',f)}
