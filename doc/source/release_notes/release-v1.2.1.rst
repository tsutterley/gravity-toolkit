.. _release-v1.2.1:

==================
`Release v1.2.1`__
==================

* ``feat``: add global dealiasing uplift program `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``refactor``: in monthly dealiasing, read data into flattened ``harmonics`` objects `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``feat``:  debug-level logging of member names and header lines `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``feat``: add ``extend_matrix`` function `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``feat``: add ``wrap_longitudes`` to ``tools`` `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``refactor``: convert ``shape`` and ``ndim`` to ``harmonics`` and ``spatial`` class properties `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``refactor``: updated inputs to ``spatial`` ``from_ascii`` function `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``fix``: scale and dimensions `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``fix``: ``spatial`` ``extend_matrix`` `(#113) <https://github.com/tsutterley/gravity-toolkit/pull/113>`_
* ``feat``: add mapping programs from Sutterley et al. 2019 and 2020 `(#114) <https://github.com/tsutterley/gravity-toolkit/pull/114>`_
* ``fix``: output file append in regress_grace_maps.py `(#114) <https://github.com/tsutterley/gravity-toolkit/pull/114>`_
* ``refactor``: merge regional plot programs into single program `(#114) <https://github.com/tsutterley/gravity-toolkit/pull/114>`_
* ``docs``: add documentation for graphing programs `(#114) <https://github.com/tsutterley/gravity-toolkit/pull/114>`_
* ``fix``: place cmap set in try/except `(#114) <https://github.com/tsutterley/gravity-toolkit/pull/114>`_
* ``refactor``: add descriptive file-level attributes to output netCDF4/HDF5 files `(#115) <https://github.com/tsutterley/gravity-toolkit/pull/115>`_
* ``feat``: include option to not compensate for elastic deformation for  `#37 <https://github.com/tsutterley/gravity-toolkit/issues/37>`_ `(#115) <https://github.com/tsutterley/gravity-toolkit/pull/115>`_
* ``feat``: added regex formatting for CNES GRGS harmonics `(#115) <https://github.com/tsutterley/gravity-toolkit/pull/115>`_
* ``feat``: added functions for getting unit attributes for known types `(#116) <https://github.com/tsutterley/gravity-toolkit/pull/116>`_
* ``docs``: improve typing for variables in docstrings `(#116) <https://github.com/tsutterley/gravity-toolkit/pull/116>`_
* ``feat``: add initial geostrophic current program `(#117) <https://github.com/tsutterley/gravity-toolkit/pull/117>`_
* ``feat``: add angular velocity of the Earth to ``units`` `(#117) <https://github.com/tsutterley/gravity-toolkit/pull/117>`_
* ``fix``: set ``case_insensitive_filename`` to ``None`` if filename is empty `(#117) <https://github.com/tsutterley/gravity-toolkit/pull/117>`_
* ``feat``: add solver option to geocenter and mascons `(#118) <https://github.com/tsutterley/gravity-toolkit/pull/118>`_
* ``feat``: allow option ``0`` in ``combine_harmonics.py`` for no unit conversion `(#118) <https://github.com/tsutterley/gravity-toolkit/pull/118>`_
* ``fix``: ``remove_label`` in ``tools.py`` `(#118) <https://github.com/tsutterley/gravity-toolkit/pull/118>`_
* ``fix``: frame in geostrophic currents notebook `(#118) <https://github.com/tsutterley/gravity-toolkit/pull/118>`_
* ``fix``: allow love numbers to be None for custom units case `(#119) <https://github.com/tsutterley/gravity-toolkit/pull/119>`_
* ``fix``: headers in plot script documentation `(#119) <https://github.com/tsutterley/gravity-toolkit/pull/119>`_
* ``feat``: add notes about shida numbers `(#119) <https://github.com/tsutterley/gravity-toolkit/pull/119>`_
* ``fix``: podaac now has different openers for s3 and data endpoints `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``feat``: add piecewise program `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``feat``: add additional fit term options to regression programs `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``docs``: remove binder links :( `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``fix``: add podaac cumulus test and skip drive `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``fix``: bump GFZ GravIS files to Release-03 `(#120) <https://github.com/tsutterley/gravity-toolkit/pull/120>`_
* ``feat``: split S2 tidal aliasing terms into GRACE and GRACE-FO eras `(#121) <https://github.com/tsutterley/gravity-toolkit/pull/121>`_
* ``feat``: add reify decorator for evaluation of properties `(#121) <https://github.com/tsutterley/gravity-toolkit/pull/121>`_
* ``refactor``: use formatting for reading from date file `(#121) <https://github.com/tsutterley/gravity-toolkit/pull/121>`_
* ``refactor``: use ``pathlib`` to define and operate on paths `(#123) <https://github.com/tsutterley/gravity-toolkit/pull/123>`_
* ``feat``: split S2 tidal aliasing terms into GRACE and GRACE-FO eras `(#123) <https://github.com/tsutterley/gravity-toolkit/pull/123>`_
* ``fix``: remove deprecated podaac drive sync programs `(#123) <https://github.com/tsutterley/gravity-toolkit/pull/123>`_
* ``docs``: remove references to deprecated sync programs `(#123) <https://github.com/tsutterley/gravity-toolkit/pull/123>`_
* ``refactor``: use fit module for getting tidal aliasing terms `(#125) <https://github.com/tsutterley/gravity-toolkit/pull/125>`_
* ``feat``: add functions to retrieve and revoke Earthdata tokens `(#125) <https://github.com/tsutterley/gravity-toolkit/pull/125>`_
* ``feat``: more operatations on spatial error if in possible data keys `(#125) <https://github.com/tsutterley/gravity-toolkit/pull/125>`_
* ``docs``: scrub PO.DAAC Drive and WebDAV `(#126) <https://github.com/tsutterley/gravity-toolkit/pull/126>`_
* ``fix``: append amplitude and phase titles when creating flags `(#127) <https://github.com/tsutterley/gravity-toolkit/pull/127>`_
* ``refactor``: modified custom units in spherical caps and disc loads `(#129) <https://github.com/tsutterley/gravity-toolkit/pull/129>`_
* ``fix``: GRACE/GRACE-FO months in drift function for  `#131 <https://github.com/tsutterley/gravity-toolkit/issues/131>`_ `(#132) <https://github.com/tsutterley/gravity-toolkit/pull/132>`_
* ``feat``: can choose different tidal aliasing periods `(#132) <https://github.com/tsutterley/gravity-toolkit/pull/132>`_
* ``feat``: add timescale class for converting between time scales `(#133) <https://github.com/tsutterley/gravity-toolkit/pull/133>`_
* ``feat``: add functions for retrieving leap seconds from NIST or IERF servers `(#133) <https://github.com/tsutterley/gravity-toolkit/pull/133>`_

.. __: https://github.com/tsutterley/gravity-toolkit/releases/tag/1.2.1
