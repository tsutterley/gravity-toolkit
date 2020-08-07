============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from PO.DAAC via https when WebDAV credentials are supplied
 - Checks MD5 hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/master/gravity_toolkit/utilities.py


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


.. method:: gravity_toolkit.utilities.ftp_list(HOST,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list


.. method:: gravity_toolkit.utilities.from_ftp(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: gravity_toolkit.utilities.from_http(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: gravity_toolkit.utilities.build_opener(username,password,urs=None)

    build urllib opener for JPL PO.DAAC Drive with supplied credentials

    Arguments:

        `username`: NASA Earthdata username

        `password`: JPL PO.DAAC Drive WebDAV password

    Keyword arguments:

        urs: Earthdata login URS 3 host


.. method:: gravity_toolkit.utilities.check_credentials()

    Check that entered JPL PO.DAAC Drive credentials are valid


.. method:: gravity_toolkit.utilities.from_podaac(HOST,username=None,password=None,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a PO.DAAC Drive https server

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: JPL PO.DAAC Drive WebDAV password

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file
