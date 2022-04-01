============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from PO.DAAC via https when WebDAV credentials are supplied
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/utilities.py


General Methods
===============


.. method:: gravity_toolkit.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        ``relpath``: local relative path as list or string


.. method:: gravity_toolkit.utilities.get_hash(local, algorithm='MD5')

    Get the hash value from a local file or BytesIO object

    Arguments:

        ``local``: BytesIO object or path to file

    Keyword Arguments:

        ``algorithm``: hashing algorithm for checksum validation

            ``'MD5'``: Message Digest

            ``'sha1'``: Secure Hash Algorithm


.. method:: gravity_toolkit.utilities.url_split(s)

    Recursively split a url path into a list

    Arguments:

        ``s``: url string


.. method:: gravity_toolkit.utilities.convert_arg_line_to_args(arg_line):

    Convert file lines to arguments

    Arguments:

        ``arg_line``: line string containing a single argument and/or comments


.. method:: gravity_toolkit.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        ``time_string``: formatted time string to parse

    Keyword arguments:

        ``format``: format for input time string


.. method:: gravity_toolkit.utilities.isoformat(time_string)

    Reformat a date string to ISO formatting

    Arguments:

        ``time_string``: formatted time string to parse


.. method:: gravity_toolkit.utilities.even(value)

    Rounds a number to an even number less than or equal to original

    Arguments:

        ``value``: number to be rounded


.. method:: gravity_toolkit.utilities.ceil(value)

    Rounds a number upward to its nearest integer

    Arguments:

        ``value``: number to be rounded upward


.. method:: gravity_toolkit.utilities.copy(source, destination, verbose=False, move=False)

    Copy or move a file with all system information

    Arguments:

        ``source``: source file

        ``destination``: copied destination file

    Keyword arguments:

        ``verbose``: print file transfer information

        ``move``: remove the source file


.. method:: gravity_toolkit.utilities.create_unique_file(filename)

    Open a unique file adding a numerical instance if existing

    Arguments:

        ``filename``: full path to output file


.. method:: gravity_toolkit.utilities.check_ftp_connection(HOST,username=None,password=None)

    Check internet connection with ftp host

    Arguments:

        ``HOST``: remote ftp host

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password


.. method:: gravity_toolkit.utilities.ftp_list(HOST,username=None,password=None,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        ``HOST``: remote ftp host path split as list

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password

        ``timeout``: timeout in seconds for blocking operations

        ``basename``: return the file or directory basename instead of the full path

        ``pattern``: regular expression pattern for reducing list

        ``sort``: sort output list

    Returns:

        ``output``: list of items in a directory

        ``mtimes``: list of last modification times for items in the directory


.. method:: gravity_toolkit.utilities.from_ftp(HOST,username=None,password=None,timeout=None,local=None,hash='',chunk=8192,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a ftp host

    Arguments:

        ``HOST``: remote ftp host path split as list

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password

        ``timeout``: timeout in seconds for blocking operations

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file

    Returns:

        ``remote_buffer``: BytesIO representation of file


.. method:: gravity_toolkit.utilities.check_connection(HOST)

    Check internet connection with an http host

    Arguments:

        ``HOST``: remote http host


.. method:: gravity_toolkit.utilities.http_list(HOST,timeout=None,context=ssl.SSLContext(),parser=lxml.etree.HTMLParser(),format='%Y-%m-%d %H:%M',pattern='',sort=False)

    List a directory on an Apache http Server

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``timeout``: timeout in seconds for blocking operations

        ``context``: SSL context for url opener object

        ``parser``: HTML parser for lxml

        ``format``: format for input time string

        ``pattern``: regular expression pattern for reducing list

        ``sort``: sort output list

    Returns:

        ``output``: list of items in a directory

        ``mtimes``: list of last modification times for items in the directory


.. method:: gravity_toolkit.utilities.from_http(HOST,timeout=None,context=ssl.SSLContext(),local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a http host

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``timeout``: timeout in seconds for blocking operations

        ``context``: SSL context for url opener object

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file

    Returns:

        ``remote_buffer``: BytesIO representation of file


.. method:: gravity_toolkit.utilities.attempt_login(urs,context=ssl.SSLContext(),password_manager=False,get_ca_certs=False,redirect=False,authorization_header=True,username=None,password=None,retries=5,netrc=None)

    attempt to build a urllib opener for NASA Earthdata

    Arguments:

        ``urs``: Earthdata login URS 3 host

    Keyword arguments:

        ``context``: SSL context for opener object

        ``password_manager``: create password manager context using default realm

        ``get_ca_certs``: get list of loaded “certification authority” certificates

        ``redirect``: create redirect handler object

        ``authorization_header``: add base64 encoded authorization header to opener

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata password

        ``retries``: number of retry attempts

        ``netrc``: path to .netrc file for authentication


.. method:: gravity_toolkit.utilities.build_opener(username,password,context=ssl.SSLContext(),password_manager=False,get_ca_certs=False,redirect=False,authorization_header=True,urs=None)

    build urllib opener for NASA Earthdata or JPL PO.DAAC Drive with supplied credentials

    Arguments:

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata or JPL PO.DAAC WebDAV password

    Keyword arguments:

        ``context``: SSL context for opener object

        ``password_manager``: create password manager context using default realm

        ``get_ca_certs``: get list of loaded “certification authority” certificates

        ``redirect``: create redirect handler object

        ``authorization_header``: add base64 encoded authorization header to opener

        ``urs``: Earthdata login URS 3 host


.. method:: gravity_toolkit.utilities.s3_client(HOST,timeout=None,region_name='us-west-2')

    get AWS s3 client for PO.DAAC Cumulus

    Keyword arguments:

        ``HOST``: PO.DAAC or ECCO AWS S3 credential host

        ``timeout``: timeout in seconds for blocking operations

        ``region_name``: AWS region name

    Returns:

        ``client``: AWS s3 client for PO.DAAC Cumulus


.. method:: gravity_toolkit.utilities.s3_key(presigned_url)

    get a s3 bucket key from a presigned url

    Arguments:

        ``presigned_url``: s3 presigned url


.. method:: gravity_toolkit.utilities.check_credentials(HOST='https://podaac-tools.jpl.nasa.gov/drive/files')

    Check that entered `JPL PO.DAAC Drive`__ credentials are valid

    Keyword arguments:

        ``HOST``: PO.DAAC or ECCO Drive host

    .. __: https://podaac-tools.jpl.nasa.gov/drive


.. method:: gravity_toolkit.utilities.drive_list(HOST,username=None,password=None,build=True,timeout=None,urs=None,parser=None,pattern='',sort=False)

    List a directory on `JPL PO.DAAC <https://podaac-tools.jpl.nasa.gov/drive>`_  or `ECCO Drive <https://ecco.jpl.nasa.gov/drive/>`_

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``username``: NASA Earthdata username

        ``password``: JPL PO.DAAC Drive WebDAV password

        ``build``: Build opener and check WebDAV credentials

        ``timeout``: timeout in seconds for blocking operations

        ``urs``: JPL PO.DAAC or ECCO login URS 3 host

        ``parser``: HTML parser for lxml

        ``pattern``: regular expression pattern for reducing list

        ``sort``: sort output list

    Returns:

        ``colnames``: list of column names in a directory

        ``collastmod``: list of last modification times for items in the directory



.. method:: gravity_toolkit.utilities.from_drive(HOST,username=None,password=None,build=True,timeout=None,urs=None,local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from `JPL PO.DAAC <https://podaac-tools.jpl.nasa.gov/drive>`_  or `ECCO Drive <https://ecco.jpl.nasa.gov/drive/>`_ https servers

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``username``: NASA Earthdata username

        ``password``: JPL PO.DAAC Drive WebDAV password

        ``build``: Build opener and check WebDAV credentials

        ``timeout``: timeout in seconds for blocking operations

        ``urs``: JPL PO.DAAC or ECCO login URS 3 host

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file

    Returns:

        ``remote_buffer``: BytesIO representation of file


.. method:: gravity_toolkit.utilities.cmr_product_shortname(mission, center, release)

    Create a list of product shortnames for CMR queries

    Arguments:

        ``mission``: GRACE (``'grace'``) or GRACE Follow-On (``'grace-fo'``)

        ``center``: GRACE/GRACE-FO processing center

        ``release``: GRACE/GRACE-FO data release

    Keyword arguments:

        ``level``: GRACE/GRACE-FO product level (``'L1A'``, ``'L1B'``, ``'L2'``)


.. method:: gravity_toolkit.utilities.cmr_readable_granules(product, solution='BA01')

    Create readable granule names pattern for CMR queries

    Arguments:

        ``product``: GRACE/GRACE-FO data product

    Keyword arguments:

        ``level``: GRACE/GRACE-FO product level (``'L1A'``, ``'L1B'``, ``'L2'``)

        ``solution``: monthly gravity field solution for Release-06

            - ``BA01``: unconstrained monthly gravity field solution to d/o 60
            - ``BB01``: unconstrained monthly gravity field solution to d/o 96
            - ``BC01``: computed monthly dealiasing solution to d/o 180


.. method:: gravity_toolkit.utilities.cmr_filter_json(search_results, endpoint='data')

    Filter the CMR json response for desired data files

    Arguments:

        ``search_results``: JSON response from CMR query

    Keyword arguments:

        ``endpoint``: url endpoint type

            - ``data``: PO.DAAC https archive
            - ``s3``: PO.DAAC Cumulus AWS S3 bucket

    Returns:

        ``granule_names``: list of GRACE/GRACE-FO granule names

        ``granule_urls``: list of GRACE/GRACE-FO granule urls

        ``granule_mtimes``: list of GRACE/GRACE-FO granule modification times


.. method:: gravity_toolkit.utilities.cmr(mission=None, center=None, release=None, level='L2', product=None, solution='BA01', start_date=None, end_date=None, provider='POCLOUD', endpoint='data', verbose=False, fid=sys.stdout)

    Query the NASA Common Metadata Repository (CMR) for GRACE/GRACE-FO data

    Keyword arguments:

        ``mission``: GRACE (``'grace'``) or GRACE Follow-On (``'grace-fo'``)

        ``center``: GRACE/GRACE-FO processing center

        ``release``: GRACE/GRACE-FO data release

        ``level``: GRACE/GRACE-FO product level (``'L1A'``, ``'L1B'``, ``'L2'``)

        ``product``: GRACE/GRACE-FO data product

        ``solution``: monthly gravity field solution for Release-06

        ``start_date``: starting date for CMR product query

        ``end_date``: ending date for CMR product query

        ``endpoint``: url endpoint type

            - ``data``: PO.DAAC https archive
            - ``s3``: PO.DAAC Cumulus AWS S3 bucket

        ``verbose``: print CMR query information

        ``fid``: open file object to print if verbose

    Returns:

        ``granule_names``: list of GRACE/GRACE-FO granule names

        ``granule_urls``: list of GRACE/GRACE-FO granule urls

        ``granule_mtimes``: list of GRACE/GRACE-FO granule modification times


.. method:: gravity_toolkit.utilities.from_figshare(directory,article='7388540',timeout=None,context=ssl.SSLContext(),chunk=16384,verbose=False,fid=sys.stdout,pattern='',mode=0o775)

    Download [Sutterley2019]_ geocenter files from `figshare`_

    Arguments:

        ``directory``: download directory

    Keyword arguments:

        ``article``: figshare article number

        ``timeout``: timeout in seconds for blocking operations

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``pattern``: regular expression pattern for reducing list

        ``mode``: permissions mode of output local file

    .. _figshare: https://doi.org/10.6084/m9.figshare.7388540


.. method:: gravity_toolkit.utilities.to_figshare(files,username=None,password=None,directory=None,timeout=None,context=ssl.SSLContext(ssl.PROTOCOL_TLS),get_ca_certs=False,verbose=False,chunk=8192)

    Send files to `figshare`_ using secure `FTP uploader <https://help.figshare.com/article/upload-large-datasets-and-bulk-upload-using-the-ftp-uploader-desktop-uploader-or-api>`_

    Arguments:

      ``files``: list of files to upload

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password

        ``directory``: figshare subdirectory for sending data

        ``timeout``: timeout in seconds for blocking operations

        ``context``: SSL context for ftp connection

        ``get_ca_certs``: get list of loaded “certification authority” certificates

        ``verbose``: print ftp transfer information

        ``chunk``: chunk size for transfer encoding

    .. _figshare: https://doi.org/10.6084/m9.figshare.7388540


.. method:: gravity_toolkit.utilities.from_csr(directory,timeout=None,context=ssl.SSLContext(),chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download `satellite laser ranging (SLR)`__ files from the University of Texas Center for Space Research (UTCSR)

    Arguments:

        ``directory``: download directory

    Keyword arguments

        ``timeout``: timeout in seconds for blocking operations

        ``context``: SSL context for url opener object

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file

    .. __: http://download.csr.utexas.edu/pub/slr/


.. method:: gravity_toolkit.utilities.from_gfz(directory,timeout=None,chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download GravIS and satellite laser ranging (SLR) files from the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ)

    Arguments:

        ``directory``: download directory

    Keyword arguments

        ``timeout``: timeout in seconds for blocking operations

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file


.. method:: gravity_toolkit.utilities.icgem_list(host='http://icgem.gfz-potsdam.de/tom_longtime',timeout=None,parser=lxml.etree.HTMLParser())

    Parse table of gravity field models on the `GFZ International Centre for Global Earth Models (ICGEM)`__ server

    Keyword arguments:

        ``host``: url for the GFZ ICGEM gravity field table

        ``timeout``: timeout in seconds for blocking operations

        ``parser``: HTML parser for lxml

    Returns:

        ``colfiles``: dictionary of static file urls mapped by field name

    .. __: http://icgem.gfz-potsdam.de/

References
##########

.. [Sutterley2019] T. C. Sutterley and I. Velicogna, "Improved Estimates of Geocenter Variability from Time-Variable Gravity and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_
