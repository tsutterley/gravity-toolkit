=========
utilities
=========

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from PO.DAAC via https when WebDAV credentials are supplied
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/utilities.py


General Methods
===============

.. autofunction:: gravity_toolkit.utilities.get_data_path

.. autofunction:: gravity_toolkit.utilities.get_hash

.. autofunction:: gravity_toolkit.utilities.get_git_revision_hash

.. autofunction:: gravity_toolkit.utilities.get_git_status

.. autofunction:: gravity_toolkit.utilities.url_split

.. autofunction:: gravity_toolkit.utilities.convert_arg_line_to_args

.. autofunction:: gravity_toolkit.utilities.get_unix_time

.. autofunction:: gravity_toolkit.utilities.isoformat

.. autofunction:: gravity_toolkit.utilities.even

.. autofunction:: gravity_toolkit.utilities.ceil

.. autofunction:: gravity_toolkit.utilities.copy

.. autofunction:: gravity_toolkit.utilities.create_unique_file

.. autofunction:: gravity_toolkit.utilities.check_ftp_connection

.. autofunction:: gravity_toolkit.utilities.ftp_list

.. autofunction:: gravity_toolkit.utilities.from_ftp

.. autofunction:: gravity_toolkit.utilities.check_connection

.. autofunction:: gravity_toolkit.utilities.http_list

.. autofunction:: gravity_toolkit.utilities.from_http

.. autofunction:: gravity_toolkit.utilities.attempt_login

.. autofunction:: gravity_toolkit.utilities.build_opener

.. autofunction:: gravity_toolkit.utilities.s3_client

.. autofunction:: gravity_toolkit.utilities.s3_bucket

.. autofunction:: gravity_toolkit.utilities.s3_key

.. autofunction:: gravity_toolkit.utilities.check_credentials

.. autofunction:: gravity_toolkit.utilities.drive_list

.. autofunction:: gravity_toolkit.utilities.from_drive

.. autofunction:: gravity_toolkit.utilities.cmr_product_shortname

.. autofunction:: gravity_toolkit.utilities.cmr_readable_granules

.. autofunction:: gravity_toolkit.utilities.cmr_filter_json

.. autofunction:: gravity_toolkit.utilities.cmr_metadata_json

.. autofunction:: gravity_toolkit.utilities.cmr

.. autofunction:: gravity_toolkit.utilities.cmr_metadata

.. autofunction:: gravity_toolkit.utilities.compile_regex_pattern

.. autofunction:: gravity_toolkit.utilities.from_figshare

.. autofunction:: gravity_toolkit.utilities.to_figshare

.. autofunction:: gravity_toolkit.utilities.from_csr

.. autofunction:: gravity_toolkit.utilities.from_gsfc

.. autofunction:: gravity_toolkit.utilities.from_gfz

.. autofunction:: gravity_toolkit.utilities.icgem_list
