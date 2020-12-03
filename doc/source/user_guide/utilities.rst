============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from PO.DAAC via https when WebDAV credentials are supplied
 - Checks MD5 hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/utilities.py


General Methods
===============


.. method:: gravity_toolkit.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        `relpath`: local relative path as list or string


.. method:: gravity_toolkit.utilities.get_hash(local)

    Get the MD5 hash value from a local file

    Arguments:

        `local`: path to file


.. method:: gravity_toolkit.utilities.url_split(s)

    Recursively split a url path into a list

    Arguments:

        `s`: url string


.. method:: gravity_toolkit.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        `time_string`: formatted time string to parse

    Keyword arguments:

        `format`: format for input time string


.. method:: gravity_toolkit.utilities.ftp_list(HOST,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `output`: list of items in a directory

        `mtimes`: list of last modification times for items in the directory


.. method:: gravity_toolkit.utilities.from_ftp(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `fid`: open file object to print if verbose

        `mode`: permissions mode of output local file


.. method:: gravity_toolkit.utilities.check_connection(HOST)

    Check internet connection

    Arguments:

        `HOST`: remote http host


.. method:: gravity_toolkit.utilities.from_http(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `fid`: open file object to print if verbose

        `mode`: permissions mode of output local file


.. method:: gravity_toolkit.utilities.build_opener(username,password,context=ssl.SSLContext(),password_manager=False,get_ca_certs=False,redirect=False,authorization_header=True,urs=None)

    build urllib opener for NASA Earthdata or JPL PO.DAAC Drive with supplied credentials

    Arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata or JPL PO.DAAC WebDAV password

    Keyword arguments:

        `context`: SSL context for opener object

        `password_manager`: create password manager context using default realm

        `get_ca_certs`: get list of loaded “certification authority” certificates

        `redirect`: create redirect handler object

        `authorization_header`: add base64 encoded authorization header to opener

        `urs`: Earthdata login URS 3 host


.. method:: gravity_toolkit.utilities.check_credentials()

    Check that entered JPL PO.DAAC Drive credentials are valid


.. method:: gravity_toolkit.utilities.podaac_list(HOST,username=None,password=None,build=True,timeout=None,parser=None,pattern='',sort=False)

    Download a file from a PO.DAAC Drive https server

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: JPL PO.DAAC Drive WebDAV password

        `build`: Build opener and check WebDAV credentials

        `timeout`: timeout in seconds for blocking operations

        `parser`: HTML parser for lxml

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list


    Returns:

        `colnames`: list of column names in a directory

        `collastmod`: list of last modification times for items in the directory



.. method:: gravity_toolkit.utilities.from_podaac(HOST,username=None,password=None,build=True,timeout=None,local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a PO.DAAC Drive https server

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: JPL PO.DAAC Drive WebDAV password

        `build`: Build opener and check WebDAV credentials

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `fid`: open file object to print if verbose

        `mode`: permissions mode of output local file


.. method:: gravity_toolkit.utilities.icgem_list(host='http://icgem.gfz-potsdam.de/tom_longtime',timeout=None,parser=lxml.etree.HTMLParser())

    Parse the table of static gravity field models on the GFZ ICGEM server

    Keyword arguments:

        `host`: url for the GFZ ICGEM gravity field table

        `timeout`: timeout in seconds for blocking operations

        `parser`: HTML parser for lxml

    Returns:

        `colfiles`: dictionary of static file urls mapped by field name
