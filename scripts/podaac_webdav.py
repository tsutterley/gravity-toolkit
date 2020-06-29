#!/usr/bin/env python
u"""
podaac_webdav.py
Written by Tyler Sutterley (05/2020)

Retrieves and prints a user's PO.DAAC WebDAV credentials

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add PO.DAAC Drive OPS to NASA Earthdata Applications and get WebDAV Password
https://podaac-tools.jpl.nasa.gov/drive

CALLING SEQUENCE:
    podaac_webdav(<USER=<username>, PASSWORD=<password>)
        or
    python podaac_webdav.py --user=<username>
    where <username> and <password> are your NASA Earthdata credentials

OUTPUTS:
    PODAAC WebDAV credentials

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user=X: username for NASA Earthdata Login
    -N X, --netrc=X: path to .netrc file for authentication
    -A, --append: append .netrc file instead of printing

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

UPDATE HISTORY:
    Written 05/2020 for public release
"""
from __future__ import print_function

import sys
import os
import re
import ssl
import time
import netrc
import getopt
import base64
import getpass
import builtins
import posixpath
import lxml.etree
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
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

#-- PURPOSE: retrieve PO.DAAC Drive WebDAV credentials
def podaac_webdav(USER, PASSWORD, parser):
    #-- https://docs.python.org/3/howto/urllib2.html#id5
    #-- create a password manager
    password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    #-- Add the username and password for NASA Earthdata Login system
    URS = 'https://urs.earthdata.nasa.gov'
    password_mgr.add_password(None,URS,USER,PASSWORD)
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
    base64_string = base64.b64encode('{0}:{1}'.format(USER,PASSWORD).encode())
    authorization_header = "Basic {0}".format(base64_string.decode())
    opener.addheaders = [("Authorization", authorization_header)]
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.
    HOST = posixpath.join('https://podaac-tools.jpl.nasa.gov','drive')
    parameters = {'client_id':'lRY01RPdFZ2BKR77Mv9ivQ', 'response_type':'code',
        'state':base64.b64encode(HOST.encode()),
        'redirect_uri':posixpath.join(HOST,'authenticated'),
        'required_scope': 'country+study_area'}
    #-- retrieve cookies from NASA Earthdata URS
    request = urllib2.Request(url=posixpath.join(URS,'oauth',
        'authorize?{0}'.format(urlencode(parameters))))
    urllib2.urlopen(request)
    #-- read and parse request for webdav password
    request = urllib2.Request(url=HOST)
    tree = lxml.etree.parse(urllib2.urlopen(request, timeout=20), parser)
    WEBDAV, = tree.xpath('//input[@id="password"]/@value')
    #-- return webdav password and cookies jar
    return (WEBDAV,cookie_jar)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -U X, --user=X\t\tUsername for NASA Earthdata Login')
    print(' -N X, --netrc=X\t\tPath to .netrc file for authentication')
    print(' -A, --append\t\tAppend .netrc file instead of printing\n')

#-- Main program that calls podaac_webdav()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','user=','netrc=','append']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hU:N:A',long_options)

    #-- command line parameters
    USER = ''
    NETRC = None
    APPEND = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-U","--user"):
            USER = arg
        elif opt in ("-N","--netrc"):
            NETRC = os.path.expanduser(arg)
        elif opt in ("-A","--append"):
            APPEND = True

    #-- NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    #-- JPL PO.DAAC drive hostname
    HOST = 'podaac-tools.jpl.nasa.gov'
    #-- get NASA Earthdata credentials
    if not USER and not NETRC:
        #-- check that NASA Earthdata credentials were entered
        USER = builtins.input('Username for {0}: '.format(URS))
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,URS))
    elif NETRC:
        USER,LOGIN,PASSWORD = netrc.netrc(NETRC).authenticators(URS)
    else:
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,URS))
    #-- if appending to netrc file and not entered: use default file
    if APPEND and not NETRC:
        NETRC = os.path.join(os.path.expanduser('~'),'.netrc')

    #-- check internet connection before attempting to run program
    if check_connection():
        #-- compile HTML parser for lxml
        parser = lxml.etree.HTMLParser()
        WEBDAV,cookie_jar = podaac_webdav(USER, PASSWORD, parser)
        #-- output to terminal or append to netrc file
        args = (USER,HOST,WEBDAV)
        if APPEND:
            #-- append to netrc file and set permissions level
            with open(NETRC,'a+') as f:
                f.write('machine {1} login {0} password {2}\n'.format(*args))
                os.chmod(NETRC, 0o600)
        else:
            print('\nWebDAV Password for {0}@{1}:\n\t{2}'.format(*args))
        #-- create "opener" (OpenerDirector instance)
        opener = urllib2.build_opener(
            urllib2.HTTPSHandler(context=ssl.SSLContext()),
            urllib2.HTTPCookieProcessor(cookie_jar))
        #-- Encode username/webdav password for request authorization headers
        base64_string = base64.b64encode('{0}:{1}'.format(USER,WEBDAV).encode())
        #-- replace Authorization header for opener
        authorization_header = "Basic {0}".format(base64_string.decode())
        opener.addheaders = [("Authorization", authorization_header)]
        #-- Now all calls to urllib2.urlopen use our opener.
        urllib2.install_opener(opener)
        #-- All calls to urllib2.urlopen will now use handler
        #-- Make sure not to include the protocol in with the URL, or
        #-- HTTPPasswordMgrWithDefaultRealm will be confused.

#-- run main program
if __name__ == '__main__':
    main()
