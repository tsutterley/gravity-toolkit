.. _release-v1.0.2.5:

=====================
`Release v1.0.2.5`__
=====================

* ``fix``: add try/except for read_GRACE_geocenter import
* ``feat``: S3 access using PO.DAAC cumulus to address  `#59 <https://github.com/tsutterley/gravity-toolkit/issues/59>`_ `(#63) <https://github.com/tsutterley/gravity-toolkit/pull/63>`_
* ``fix``: update readable granule for L1A/B grav
* ``fix``: for now use podaac drive provider
* ``docs``: updated docstrings to numpy documentation format `(#65) <https://github.com/tsutterley/gravity-toolkit/pull/65>`_
* ``docs``: use autodoc to build documentation `(#65) <https://github.com/tsutterley/gravity-toolkit/pull/65>`_
* ``feat``: new internal ncdf/hdf5 read/write within harmonics and spatial classes `(#65) <https://github.com/tsutterley/gravity-toolkit/pull/65>`_
* ``feat``: add citations and references to read_GIA_model `(#65) <https://github.com/tsutterley/gravity-toolkit/pull/65>`_
* ``docs``: update headers to remove deprecated ncdf/hdf5 read/write modules `(#65) <https://github.com/tsutterley/gravity-toolkit/pull/65>`_
* ``refactor``: moved Load love number wrapper function to within read `(#66) <https://github.com/tsutterley/gravity-toolkit/pull/66>`_
* ``feat``: change badge for pangeo to aws us-west-2 for podaac cloud access `(#66) <https://github.com/tsutterley/gravity-toolkit/pull/66>`_
* ``fix``: pin markupsafe to 2.0.1 to prevent soft_unicode error `(#66) <https://github.com/tsutterley/gravity-toolkit/pull/66>`_
* ``fix``: change function name back to load_love_numbers `(#67) <https://github.com/tsutterley/gravity-toolkit/pull/67>`_
* ``docs``: update history in headers `(#67) <https://github.com/tsutterley/gravity-toolkit/pull/67>`_
* ``test``: try docker build with windows to address  `#64 <https://github.com/tsutterley/gravity-toolkit/issues/64>`_ `(#68) <https://github.com/tsutterley/gravity-toolkit/pull/68>`_
* ``fix``: include utf-8 encoding in reads to be windows compliant `(#68) <https://github.com/tsutterley/gravity-toolkit/pull/68>`_
* ``fix``:  only sync newsletters for mission of interest
* ``fix``: spatial field mapping for output
* ``fix``: mask in sea level equation
* ``feat``: add from_GIA to harmonics class `(#69) <https://github.com/tsutterley/gravity-toolkit/pull/69>`_
* ``feat``: prepare CMR queries for version 1 of RL06 `(#69) <https://github.com/tsutterley/gravity-toolkit/pull/69>`_
* ``fix``: expansion and squeezing of mask variable if None `(#69) <https://github.com/tsutterley/gravity-toolkit/pull/69>`_
* ``feat``: add option for L2 version in sync programs `(#69) <https://github.com/tsutterley/gravity-toolkit/pull/69>`_
* ``fix``: improved passing of filename attribute in harmonic objects `(#70) <https://github.com/tsutterley/gravity-toolkit/pull/70>`_
* ``fix``: include filename when copying spatial objects `(#70) <https://github.com/tsutterley/gravity-toolkit/pull/70>`_
* ``refactor``: always try syncing from both grace and grace-fo missions `(#71) <https://github.com/tsutterley/gravity-toolkit/pull/71>`_
* ``feat``: added AW13 models using IJ05-R2 ice history `(#71) <https://github.com/tsutterley/gravity-toolkit/pull/71>`_
* ``feat``: allow input ascii harmonic files to have additional columns `(#71) <https://github.com/tsutterley/gravity-toolkit/pull/71>`_
* ``docs``: update environment file `(#71) <https://github.com/tsutterley/gravity-toolkit/pull/71>`_
* ``fix``: index using granules `(#71) <https://github.com/tsutterley/gravity-toolkit/pull/71>`_

.. __: https://github.com/tsutterley/gravity-toolkit/releases/tag/v1.0.2.5
