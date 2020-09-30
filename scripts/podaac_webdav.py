#!/usr/bin/env python
u"""
podaac_webdav.py
Written by Tyler Sutterley (09/2020)

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
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -A, --append: append .netrc file instead of printing

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

UPDATE HISTORY:
    Updated 09/2020: use argparse to set command line parameters
    Written 05/2020 for public release
"""
from __future__ import print_function

import sys
import os
import netrc
import base64
import getpass
import builtins
import argparse
import posixpath
import lxml.etree
import gravity_toolkit.utilities

#-- PURPOSE: retrieve PO.DAAC Drive WebDAV credentials
def podaac_webdav(USER, PASSWORD, parser):
    #-- build opener for retrieving PO.DAAC Drive WebDAV credentials
    #-- Add the username and password for NASA Earthdata Login system
    URS = 'https://urs.earthdata.nasa.gov'
    gravity_toolkit.utilities.build_opener(USER, PASSWORD,
        password_manager=True, authorization_header=True, urs=URS)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.
    HOST = posixpath.join('https://podaac-tools.jpl.nasa.gov','drive')
    parameters = gravity_toolkit.utilities.urlencode(
        {'client_id':'lRY01RPdFZ2BKR77Mv9ivQ', 'response_type':'code',
        'state':base64.b64encode(HOST.encode()),
        'redirect_uri':posixpath.join(HOST,'authenticated'),
        'required_scope': 'country+study_area'}
    )
    #-- retrieve cookies from NASA Earthdata URS
    request = gravity_toolkit.utilities.urllib2.Request(
        url=posixpath.join(URS,'oauth','authorize?{0}'.format(parameters)))
    gravity_toolkit.utilities.urllib2.urlopen(request)
    #-- read and parse request for webdav password
    request = gravity_toolkit.utilities.urllib2.Request(url=HOST)
    response = gravity_toolkit.utilities.urllib2.urlopen(request,timeout=20)
    tree = lxml.etree.parse(response, parser)
    WEBDAV, = tree.xpath('//input[@id="password"]/@value')
    #-- return webdav password
    return WEBDAV

#-- Main program that calls podaac_webdav()
def main():
   #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Retrieves and prints a user's PO.DAAC WebDAV
            credentials
            """
    )
    #-- command line parameters
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default='',
        help='Username for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=os.path.expanduser, default='',
        help='Path to .netrc file for authentication')
    #-- append to netrc
    parser.add_argument('--append','-A',
        default=False, action='store_true',
        help='Append .netrc file instead of printing')
    args = parser.parse_args()

    #-- NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    #-- JPL PO.DAAC drive hostname
    HOST = 'podaac-tools.jpl.nasa.gov'
    #-- get NASA Earthdata credentials
    if not args.user and not args.netrc:
        #-- check that NASA Earthdata credentials were entered
        USER = builtins.input('Username for {0}: '.format(URS))
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,URS))
    elif args.netrc:
        NETRC = args.netrc
        USER,LOGIN,PASSWORD = netrc.netrc(NETRC).authenticators(URS)
    else:
        #-- enter password securely from command-line
        USER = args.user
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,URS))
    #-- if appending to netrc file and not entered: use default file
    if args.append and not args.netrc:
        NETRC = os.path.join(os.path.expanduser('~'),'.netrc')

    #-- check internet connection before attempting to run program
    DRIVE = posixpath.join('https://podaac-tools.jpl.nasa.gov','drive')
    if gravity_toolkit.utilities.check_connection(DRIVE):
        #-- compile HTML parser for lxml
        WEBDAV = podaac_webdav(USER, PASSWORD, lxml.etree.HTMLParser())
        #-- output to terminal or append to netrc file
        a = (USER,HOST,WEBDAV)
        if args.append:
            #-- append to netrc file and set permissions level
            with open(NETRC,'a+') as f:
                f.write('machine {1} login {0} password {2}\n'.format(*a))
                os.chmod(NETRC, 0o600)
        else:
            print('\nWebDAV Password for {0}@{1}:\n\t{2}'.format(*a))

#-- run main program
if __name__ == '__main__':
    main()
