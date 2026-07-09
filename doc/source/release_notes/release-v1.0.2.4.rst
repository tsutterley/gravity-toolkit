.. _release-v1.0.2.4:

=====================
`Release v1.0.2.4`__
=====================

* ``fix``: uncertainties for SLR CS21 and CS22
* ``feat``: add option for setting input format of the mascon files `(#46) <https://github.com/tsutterley/gravity-toolkit/pull/46>`_
* ``feat``: add averaging kernel program `(#46) <https://github.com/tsutterley/gravity-toolkit/pull/46>`_
* ``feat``: time-variable gravity data from COST-G, GRAZ, SWARM  `#37 <https://github.com/tsutterley/gravity-toolkit/issues/37>`_ `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``feat``: time-variable gravity data from COST-G GRAZ, SWARM  `#37 <https://github.com/tsutterley/gravity-toolkit/issues/37>`_ `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``feat``: need to add sync programs for all products `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``refactor``: call read ICGEM from read gfc `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``fix``: Swarm strings and titles `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``docs``: update documentation to add additional information `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``feat``: add sync programs for Swarm and GRACE COST-G `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``feat``: add ITSG GRAZ GRACE sync `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``fix``: accidental copy paste `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``fix``: output index in separate loop for COST-G `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``test``: add tests for COST-G, GRAZ and Swarm gfc files `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``test``: use fixture to download geocenter files `(#47) <https://github.com/tsutterley/gravity-toolkit/pull/47>`_
* ``refactor``: merged integration and fourier harmonics programs ``fix``: use fill values for input ascii files in convert_harmonics ``fix``: update grid attributes after allocating for data in combine_harmonics ``ci``: install proj from source for cartopy dependency ``ci``: install devel cartopy from repo `(#48) <https://github.com/tsutterley/gravity-toolkit/pull/48>`_
* ``refactor``: adding new tools module to simplify notebooks `(#49) <https://github.com/tsutterley/gravity-toolkit/pull/49>`_
* ``feat``: adding new time functions for to/from grace months `(#49) <https://github.com/tsutterley/gravity-toolkit/pull/49>`_
* ``docs``: add documentation for tools `(#49) <https://github.com/tsutterley/gravity-toolkit/pull/49>`_
* ``fix``: adjust months if final in time series `(#49) <https://github.com/tsutterley/gravity-toolkit/pull/49>`_
* ``fix``: adjustable minimum in gaussian weights
* ``refactor``: slim requirements `(#50) <https://github.com/tsutterley/gravity-toolkit/pull/50>`_
* ``refactor``: using python logging for handling verbose output `(#51) <https://github.com/tsutterley/gravity-toolkit/pull/51>`_
* ``feat``: add version program and add package __version__ attribute `(#51) <https://github.com/tsutterley/gravity-toolkit/pull/51>`_
* ``feat``: add more remove file options in GRACE maps `(#52) <https://github.com/tsutterley/gravity-toolkit/pull/52>`_
* ``fix``: remove file choices in calc mascon `(#52) <https://github.com/tsutterley/gravity-toolkit/pull/52>`_
* ``fix``: logging ``CRITICAL``
* ``fix``: logging with filenames
* ``feat``: grid conversion routines for publicly available mascon solutions `(#53) <https://github.com/tsutterley/gravity-toolkit/pull/53>`_
* ``test``: attempt to download cost-g from http `(#53) <https://github.com/tsutterley/gravity-toolkit/pull/53>`_
* ``refactor``: use logging to print multiprocessing exceptions
* ``refactor``: netCDF4 and HDF5 input programs `(#54) <https://github.com/tsutterley/gravity-toolkit/pull/54>`_
* ``ci``: update for geoid-toolkit dependencies `(#54) <https://github.com/tsutterley/gravity-toolkit/pull/54>`_
* ``docs``: pin docutils to 0.18 `(#54) <https://github.com/tsutterley/gravity-toolkit/pull/54>`_
* ``fix``: modify legendre normalization to prevent high degree overflows `(#55) <https://github.com/tsutterley/gravity-toolkit/pull/55>`_
* ``test``: add legendre test `(#55) <https://github.com/tsutterley/gravity-toolkit/pull/55>`_
* ``test``: add associated legendre polynomial test `(#56) <https://github.com/tsutterley/gravity-toolkit/pull/56>`_
* ``fix``: HDF5 and netCDF4 io for single time case `(#56) <https://github.com/tsutterley/gravity-toolkit/pull/56>`_
* ``fix``: format for index in notebooks `(#56) <https://github.com/tsutterley/gravity-toolkit/pull/56>`_
* ``feat``: use new GSFC weekly 5x5s for CS2 and C50 `(#57) <https://github.com/tsutterley/gravity-toolkit/pull/57>`_
* ``refactor``: create geocenter program for reading and operating `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``test``: include new geocenter class `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``feat``: add load love numbers for unit conversion `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``feat``: add mean function to geocenter `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``refactor``: rename SLF geocenter reader `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``feat``: add more geocenter operations `(#58) <https://github.com/tsutterley/gravity-toolkit/pull/58>`_
* ``feat``: add geocenter figure programs from Sutterley and Velicogna (2019)
* ``feat``: update geocenter plot programs for AGU `(#60) <https://github.com/tsutterley/gravity-toolkit/pull/60>`_
* ``feat``: added netCDF4 reader for UCI iteration files `(#60) <https://github.com/tsutterley/gravity-toolkit/pull/60>`_
* ``feat``: added custom colormap function for some common scales `(#60) <https://github.com/tsutterley/gravity-toolkit/pull/60>`_
* ``feat``: added UNITS list option for converting from custom units `(#60) <https://github.com/tsutterley/gravity-toolkit/pull/60>`_
* ``feat``: option to specify a specific geocenter correction file `(#61) <https://github.com/tsutterley/gravity-toolkit/pull/61>`_
* ``feat``: can use variable loglevels for verbose output `(#61) <https://github.com/tsutterley/gravity-toolkit/pull/61>`_
* ``fix``: fix default file prefix to include center and release information `(#61) <https://github.com/tsutterley/gravity-toolkit/pull/61>`_

.. __: https://github.com/tsutterley/gravity-toolkit/releases/tag/v1.0.2.4
