==============
NASA Earthdata
==============

NASA Data Distribution Centers
##############################

The NASA Earth Science Data Information Systems Project funds and operates
`12 Distributed Active Archive Centers (DAACs) <https://earthdata.nasa.gov/about/daacs>`_ throughout the United States.
These centers have recently transitioned from ftp to https servers.
The https updates are designed to increase performance and improve security during data retrieval.
NASA Earthdata uses `OAuth2 <https://wiki.earthdata.nasa.gov/pages/viewpage.action?pageId=71700485>`_, an approach to authentication that protects your personal information.

- https://urs.earthdata.nasa.gov/documentation
- https://wiki.earthdata.nasa.gov/display/EL/Knowledge+Base

PO.DAAC
#######
The `Physical Oceanography Distributed Active Archive Center (PO.DAAC) <https://podaac.jpl.nasa.gov/>`_
provides data and related information pertaining to the physical processes and conditions of the global oceans,
including measurements of ocean winds, temperature, topography, salinity, circulation and currents, and sea ice.
PO.DAAC hosts

PO.DAAC is `migrating its data archive to the Earthdata Cloud <https://podaac.jpl.nasa.gov/cloud-datasets/migration>`_,
which is hosted in Amazon Web Services (AWS).
GRACE/GRACE-FO spherical harmonic data is currently scheduled to be "cloud enabled"
during Phase 4 and Phase 5 of the transition, slated for April and July 2022 respectively.

If any problems contact JPL PO.DAAC support at `podaac@podaac.jpl.nasa.gov <mailto:podaac@podaac.jpl.nasa.gov>`_
or the NASA EOSDIS support team `support@earthdata.nasa.gov <mailto:support@earthdata.nasa.gov>`_.

WebDAV
------
JPL PO.DAAC Drive uses passwords generated using the Web Distributed Authoring and Versioning (WebDAV) API.
This password is created at the `PO.DAAC Drive <https://podaac-tools.jpl.nasa.gov/drive/>`_ website.
Use this password rather than your Earthdata password when retrieving model data from PO.DAAC Drive.
`Click here for more information <https://podaac-tools.jpl.nasa.gov/drive/help>`_.

Steps to Sync from PO.DAAC
--------------------------

1. `Register with NASA Earthdata Login system <https://urs.earthdata.nasa.gov/users/new>`_
2. `After registering, login to the system <https://urs.earthdata.nasa.gov/home>`_
3. Add ``PO.DAAC Drive OPS`` `applications to Earthdata <https://wiki.earthdata.nasa.gov/display/EL/How+To+Pre-authorize+an+application>`_
4. Retrieve your `JPL WebDAV password <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_webdav.py>`_ to access PO.DAAC servers
5. `Sync time-variable gravity data <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_grace_sync.py>`_

Can also create a ``.netrc`` file for permanently storing NASA Earthdata and JPL WebDAV credentials:

.. code-block:: bash

    echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> ~/.netrc
    echo "machine podaac-tools.jpl.nasa.gov login <uid> password <webdav>" >> ~/.netrc
    chmod 0600 ~/.netrc

Or set environmental variables for your NASA Earthdata and JPL WebDAV credentials:

.. code-block:: bash

    export EARTHDATA_USERNAME=<uid>
    export EARTHDATA_PASSWORD=<password>
    export PODAAC_PASSWORD=<webdav>

NASA Common Metadata Repository
###############################

The NASA Common Metadata Repository (CMR) is a catalog of all data
and service metadata records contained as part of NASA's Earth
Observing System Data and Information System (EOSDIS).
Querying the CMR system is a way of quickly performing a search
through the NASA Earthdata archive.
Basic queries for the granule names, PO.DAAC URLs and modification times
of GRACE/GRACE-FO data are available through the ``cmr`` routine in the
``utilities`` module.
For AWS instances in ``us-west-2``, CMR queries can access urls for S3 endpoints.

.. code-block:: python

   ids,urls,mtimes = gravity_toolkit.utilities.cmr(mission='grace-fo',
      center='JPL', release='RL06', level='L2', product='GSM',
      solution='BA01', provider='POCLOUD', endpoint='s3', verbose=False)

Other Data Access Examples
##########################
- `Curl and Wget <https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget>`_
- `Python <https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python>`_
