.. _release-v1.2.0:

==================
`Release v1.2.0`__
==================

* ``refactor``: new repo name following  `#84 <https://github.com/tsutterley/gravity-toolkit/issues/84>`_ `(#99) <https://github.com/tsutterley/gravity-toolkit/pull/99>`_
* ``docs``: update docstrings `(#100) <https://github.com/tsutterley/gravity-toolkit/pull/100>`_
* ``fix``: verify warnings have type and filtered after `(#101) <https://github.com/tsutterley/gravity-toolkit/pull/101>`_
* ``feat``: added function to retrieve named units `(#102) <https://github.com/tsutterley/gravity-toolkit/pull/102>`_
* ``fix``: set custom units as top option in if/else statements `(#102) <https://github.com/tsutterley/gravity-toolkit/pull/102>`_
* ``feat``: added wrapper function for smoothing and converting to output units `(#102) <https://github.com/tsutterley/gravity-toolkit/pull/102>`_
* ``test``: add test for getting units `(#102) <https://github.com/tsutterley/gravity-toolkit/pull/102>`_
* ``fix``: expand harmonics case where data is a single degree `(#103) <https://github.com/tsutterley/gravity-toolkit/pull/103>`_
* ``fix``: harmonics case where maximum spherical harmonic degree is 0 `(#103) <https://github.com/tsutterley/gravity-toolkit/pull/103>`_
* ``fix``: units case where maximum spherical harmonic degree is 0 `(#103) <https://github.com/tsutterley/gravity-toolkit/pull/103>`_
* ``fix``: load love number case with maximum spherical harmonic degree is 0 `(#103) <https://github.com/tsutterley/gravity-toolkit/pull/103>`_
* ``feat``: add PREM hard and soft sediment love numbers `(#104) <https://github.com/tsutterley/gravity-toolkit/pull/104>`_
* ``feat``: data class for load love numbers with attributes for model `(#105) <https://github.com/tsutterley/gravity-toolkit/pull/105>`_
* ``feat``: customizable file-level attributes to netCDF4 and HDF5 `(#106) <https://github.com/tsutterley/gravity-toolkit/pull/106>`_
* ``feat``: add root attributes to output netCDF4 and HDF5 files `(#107) <https://github.com/tsutterley/gravity-toolkit/pull/107>`_
* ``fix``: only attempt to squeeze from final dimension in harmonics objects `(#108) <https://github.com/tsutterley/gravity-toolkit/pull/108>`_
* ``feat``: add indexing of filenames to object iterators `(#109) <https://github.com/tsutterley/gravity-toolkit/pull/109>`_
* ``feat``: new ``scaling_factors`` inheritance of ``spatial`` class `(#110) <https://github.com/tsutterley/gravity-toolkit/pull/110>`_
* ``refactor``: simplified unit factors in spherical caps and disc loads `(#111) <https://github.com/tsutterley/gravity-toolkit/pull/111>`_
* ``fix``: case insensitive searching for HDF5 and netCDF4 GIA attributes `(#111) <https://github.com/tsutterley/gravity-toolkit/pull/111>`_
* ``docs``: slimmer builds by placing imports within try/except `(#111) <https://github.com/tsutterley/gravity-toolkit/pull/111>`_
* ``refactor``: remove deprecated functions `(#112) <https://github.com/tsutterley/gravity-toolkit/pull/112>`_

.. __: https://github.com/tsutterley/gravity-toolkit/releases/tag/1.2.0
