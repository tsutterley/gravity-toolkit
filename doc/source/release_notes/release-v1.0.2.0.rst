.. _release-v1.0.2.0:

=====================
`Release v1.0.2.0`__
=====================

* ``docs``: Update documentation to use sphinx and readthedocs
* ``docs``: update ``readthedocs.yml`` to use conda
* ``fix``: update ``environment.yml`` to include pip
* ``docs``: different environment for docs
* ``feat``: add AOD1B oblateness and isomorphic parameters
* ``feat``: add netCDF4 and HDF5 options to read_GIA_model.py
* ``feat``: add GIA step to getting started
* ``feat``: added ``harmonics`` class for correcting GRACE/GRACE-FO data
* ``feat``: add ``mean`` and ``from_dict`` to ``harmonics`` class
* ``feat``: use ``harmonics`` class to ``index``, ``add``, ``subtract``, ``convolve`` and ``destripe``
* ``refactor``: separate ``units`` into its own class
* ``feat``: ``units`` factors for ``harmonics`` and ``spatial``
* ``docs``: add header notes to notebook describing GRACE/GRACE-FO measurements
* ``fix``: podaac program default release to RL06
* ``fix``: enumeration in ``GRACE-Data-File-Formats.md``
* ``feat``: add gravity model read from GFZ ICGEM
* ``feat``: add more functionality and mathematical functions to ``harmonics`` class
* ``feat``: check dimensions of ``harmonics`` objects if using mathematics functions
* ``feat``: add viscoelastic crustal uplift to ``units`` class for converting GIA rates
* ``docs``: include source code links to documents
* ``fix``: include degree and order in read GIA program
* ``docs``: include level-2 handbooks and processing standard documents
* ``docs``: include notes about pole tide and atmospheric corrections
* ``feat``: add ocean redistribution to ipynb, add spherical harmonic calculations
* ``feat``: add ocean redistribution to ipynb, add spherical harmonic calculations
* ``feat``: add date option to ncdf/hdf5 write programs
* ``fix``: in ``harmonics`` ``from_list`` separate date and sort
* ``feat``: get ``harmonics`` third dimension from shape
* ``feat``: add options to ``flatten`` and ``expand`` ``harmonics`` matrices or arrays
* ``docs``: updated ``README.md`` to have data repositories
* ``docs``: add more notes about spatial units and conversion from harmonics
* ``feat``: include file of load love numbers `(Han and Wahr, 1995) <http://doi.wiley.com/10.1111/j.1365-246X.1995.tb01819.x>`_
* ``docs``: update html link to `Martin Mohlenkamp's uguide <http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf>`_
* ``docs``: add link to love numbers to getting started doc
* ``feat``: add more legendre and harmonics programs
* ``docs``: updated readme and documentation
* ``docs``: update html links to https
* ``feat``: add webdav program for retrieving PO.DAAC credentials
* ``feat``: add ``netrc`` option to PO.DAAC sync program
* ``feat``: ``netrc`` file can be appended from webdav program
* ``docs``: updated README, getting started and program documentation
* ``feat``: add additional date and conversion programs
* ``feat``: output list of filenames if using ``from_list()`` in ``harmonics``
* ``fix``: subset and index can output the harmonics filename if set
* ``fix``: increase timeout to 2 minutes in ``podaac_grace_sync.py``
* ``docs``: Add install with pip from git to readme
* ``docs``: update readme to note CC4 license for non-code content
* ``feat``: Add ``spatial`` class for reading, writing and processing grids
* ``feat``: additional capablities within ``spatial`` class
* ``docs``: update ``spatial`` class documentation
* ``refactor``: reorganize code structure
* ``docs``: update documentation for new code structure
* ``fix``: use dependencies from ``requirements.txt`` in ``setup.py``
* ``fix``: include files within scripts directory in ``setup.py``
* ``feat``: add CLI spatial and regression programs
* ``feat``: add ``zeros_like`` to ``harmonics`` class
* ``docs``: update documentation for added modules
* ``docs``: update readme for added modules
* ``feat``: add GRACE spatial error program
* ``docs``: update documentation and readme for changes
* ``feat``: update setup to mark new version
* ``docs``: add note about John to citations
* ``docs``: add function docstrings
* ``feat``: add back level-1b dealiasing sync programs
* ``docs``: update documentation and readme
* ``feat``: add case insensitive file search
* ``feat``: update mask if no ``fill_value`` in ``spatial``
* ``feat``: update jupyter notebook to use ``spatial`` class in plot
* ``ci``: add github actions for continuous integration
* ``ci``: ``flake8`` linter updates for CI
* ``fix``: add ``scipy`` to ``environment.yml`` and ``requirements.txt``
* ``feat``: add github dependency for ``geocenter`` to requirements
* ``fix``: remove dependency links in lieu of requirements
* ``ci``: add github actions for continuous integration
* ``test``: add spherical harmonic conversion test
* ``fix``: update regular expressions for ``flake8`` compat
* ``test``: add spherical harmonic conversion test
* ``test``: add test for downloading and reading GRACE data
* ``refactor``: move build opener to ``utilities`` routines
* ``docs``: update readme and documentation for utilities
* ``test``: add more tests for downloading and reading GRACE data
* ``feat``: use ``podaac_list()`` within ``podaac_grace_sync`` program
* ``feat``: ``read_GRACE_harmonics()`` can read ``bytesIO`` objects
* ``feat``: added GFZ ftp download and read test
* ``feat``: added compression options to ``harmonic`` and ``spatial`` file input
* ``fix``: ``flake8`` updates for ``python3``
* ``fix``: update legendre polynomial programs for divide by zero in differentials
* ``fix``: add ``KeyError`` to ``from_dict``
* ``fix``: ``flake8`` updates for python3
* ``feat``: use ``utilities`` to define path to load love numbers file
* ``feat``: include data in package
* ``feat``: Update ``MANIFEST.in`` for included data
* ``ci``: add ``macos-latest`` to testing strategy
* ``ci``: will use homebrew package manager to install dependencies
* ``feat``: add podaac sync within jupyter notebook with magics
* ``refactor``: reorganize base directory: .binder and notebooks
* ``test``: calculate test coverage
* ``test``: upload coverage file in github actions
* ``feat``: include GSFC GRACE mascons in dates
* ``fix``: update python language support
* ``fix``: use ``urllib`` from ``gravity_toolkit`` ``utilities``
* ``feat``: generalize build opener for different earthdata instances
* ``refactor``: switching to main branch as primary
* ``chore``: update links to main branch in readme and docs
* ``feat``: use ``argparse`` to set parameters
* ``feat``: ``abspath`` and ``expanduser`` in ``argparse`` paths
* ``feat``: add ``spatial`` ascii header option
* ``fix``: update ``spatial`` ``mean`` to catch more exceptions
* ``feat``: add more routines to ``spatial`` class `(#17) <https://github.com/tsutterley/gravity-toolkit/pull/17>`_
* ``refactor``: update podaac programs to simplify args `(#17) <https://github.com/tsutterley/gravity-toolkit/pull/17>`_
* ``feat``: add updated CNES sync program `(#18) <https://github.com/tsutterley/gravity-toolkit/pull/18>`_
* ``feat``: add GFZ ICGEM list for static models `(#18) <https://github.com/tsutterley/gravity-toolkit/pull/18>`_
* ``docs``: update documentation `(#18) <https://github.com/tsutterley/gravity-toolkit/pull/18>`_
* ``feat``: added more love number options and from gfc for mean files `(#19) <https://github.com/tsutterley/gravity-toolkit/pull/19>`_
* ``feat``: add `Sutterley and Velicogna geocenter <https://doi.org/10.6084/m9.figshare.7388540>`_ download `(#19) <https://github.com/tsutterley/gravity-toolkit/pull/19>`_
* ``feat``: updated SLR geocenter for new solutions from Minkang Cheng `(#19) <https://github.com/tsutterley/gravity-toolkit/pull/19>`_
* ``feat``: added download for satellite laser ranging (SLR) files from UTCSR `(#19) <https://github.com/tsutterley/gravity-toolkit/pull/19>`_
* ``feat``: add first public versions of mascon programs `(#20) <https://github.com/tsutterley/gravity-toolkit/pull/20>`_
* ``docs``: add documentation outlining programs `(#20) <https://github.com/tsutterley/gravity-toolkit/pull/20>`_
* ``docs``: add documentation outlining grace/grace-fo processing `(#20) <https://github.com/tsutterley/gravity-toolkit/pull/20>`_
* ``docs``: add blurbs to add Yara's comments `(#20) <https://github.com/tsutterley/gravity-toolkit/pull/20>`_
* ``docs``: use restructuredtext for background `(#20) <https://github.com/tsutterley/gravity-toolkit/pull/20>`_
* ``refactor``: generalize utilities for downloading from JPL drive (PO.DAAC/ECCO) `(#21) <https://github.com/tsutterley/gravity-toolkit/pull/21>`_
* ``feat``: can calculate means (``spatial`` and ``harmonic``) for a subset `(#21) <https://github.com/tsutterley/gravity-toolkit/pull/21>`_
* ``feat``: add pressure harmonics routines for OBP/surface pressure `(#21) <https://github.com/tsutterley/gravity-toolkit/pull/21>`_
* ``fix``: update requirements `(#21) <https://github.com/tsutterley/gravity-toolkit/pull/21>`_
* ``feat``: added ``time`` module to be able to convert delta times `(#22) <https://github.com/tsutterley/gravity-toolkit/pull/22>`_
* ``refactor``: merged ``convert_calendar_decimal`` and ``convert_julian`` with ``time`` module `(#22) <https://github.com/tsutterley/gravity-toolkit/pull/22>`_
* ``feat``: update netCDF4 and HDF5 programs for attributes and references `(#22) <https://github.com/tsutterley/gravity-toolkit/pull/22>`_
* ``docs``: update documentation `(#22) <https://github.com/tsutterley/gravity-toolkit/pull/22>`_
* ``test``: add test module for time programs `(#22) <https://github.com/tsutterley/gravity-toolkit/pull/22>`_
* ``feat``: update netCDF and HDF5 programs to read from memory `(#23) <https://github.com/tsutterley/gravity-toolkit/pull/23>`_
* ``feat``: update ``ftp_list`` and read for protected ftp `(#23) <https://github.com/tsutterley/gravity-toolkit/pull/23>`_
* ``test``: add ftp connection check `(#23) <https://github.com/tsutterley/gravity-toolkit/pull/23>`_
* ``refactor``: update ftp programs to use ``utilities`` `(#24) <https://github.com/tsutterley/gravity-toolkit/pull/24>`_
* ``feat``: add even rounding utility `(#24) <https://github.com/tsutterley/gravity-toolkit/pull/24>`_
* ``refactor``: moved pressure harmonics function to ``model_harmonics`` `(#24) <https://github.com/tsutterley/gravity-toolkit/pull/24>`_
* ``feat``: use ``harmonics`` class as output from SH generators `(#25) <https://github.com/tsutterley/gravity-toolkit/pull/25>`_
* ``feat``: add piecewise regression routine for breakpoint analysis `(#26) <https://github.com/tsutterley/gravity-toolkit/pull/26>`_
* ``docs``: add harmonic triangle plot notebook
* ``docs``: add regression plots to spatial map notebook
* ``docs``: use ``sphinx_rtd_theme`` for documentation
* ``docs``: change some markdown docs to rst
* ``feat``: add date parser for cases when only a date and no units
* ``docs``: add badges to examples documentation
* ``feat``: add kfactor calculation program to ``spatial`` class
* ``feat``: add degree amplitude function to ``harmonics`` class
* ``fix``: prevent warnings with python3 compatible regex strings in nc/hdf5 read
* ``test``: add point mass test
* ``fix``: modify legendre case with underflow
* ``docs``: update references in point harmonics programs
* ``feat``: added ``replace_masked`` to replace masked values in ``spatial`` data
* ``fix``: in ``spatial`` broadcast mask over third dimension
* ``refactor``: changed remove index to files with specified formats `(#27) <https://github.com/tsutterley/gravity-toolkit/pull/27>`_
* ``feat``: added generic reader, generic writer and write to list functions `(#27) <https://github.com/tsutterley/gravity-toolkit/pull/27>`_
* ``feat``: added ``adjust_months`` function to fix "special" months cases `(#27) <https://github.com/tsutterley/gravity-toolkit/pull/27>`_
* ``feat``: include geocenter read program for coefficents from Sean `(#27) <https://github.com/tsutterley/gravity-toolkit/pull/27>`_
* ``fix``: replaced ``numpy`` bool to prevent deprecation warning `(#27) <https://github.com/tsutterley/gravity-toolkit/pull/27>`_
* ``refactor``: generalize kwargs to ascii, netCDF4 and HDF5 readers and writers
* ``refactor``: moved model mascon programs to ``model_harmonics``
* ``refactor``: merged read ICGEM harmonics with ``geoid_toolkit`` reader
* ``docs``: add contribution guidelines `(#28) <https://github.com/tsutterley/gravity-toolkit/pull/28>`_
* ``docs``: more documentation standardization `(#28) <https://github.com/tsutterley/gravity-toolkit/pull/28>`_
* ``ci``: remove python 3.5 from tests `(#28) <https://github.com/tsutterley/gravity-toolkit/pull/28>`_
* ``docs``: documentation standardization
* ``docs``: update documentation `(#29) <https://github.com/tsutterley/gravity-toolkit/pull/29>`_
* ``feat``: set a default netrc file and check access `(#29) <https://github.com/tsutterley/gravity-toolkit/pull/29>`_
* ``feat``: default credentials from environmental variables `(#29) <https://github.com/tsutterley/gravity-toolkit/pull/29>`_
* ``docs``: add rst format citations to documentation `(#30) <https://github.com/tsutterley/gravity-toolkit/pull/30>`_
* ``fix``: update CSR SLR function (thanks @hulecom for pointing out the file format change) `(#30) <https://github.com/tsutterley/gravity-toolkit/pull/30>`_
* ``fix``: update setup file to check if ``readthedocs`` `(#30) <https://github.com/tsutterley/gravity-toolkit/pull/30>`_
* ``feat``: adding more SLR low-degree replacements `(#31) <https://github.com/tsutterley/gravity-toolkit/pull/31>`_
* ``docs``: update documentation for SLR harmonics `(#31) <https://github.com/tsutterley/gravity-toolkit/pull/31>`_
* ``feat``: add parser object for removing commented or empty lines `(#32) <https://github.com/tsutterley/gravity-toolkit/pull/32>`_
* ``feat``: add GFZ SLR solutions for C20/C21+S21/C30 `(#33) <https://github.com/tsutterley/gravity-toolkit/pull/33>`_
* ``feat``: add GFZ GravIS geocenter solutions `(#33) <https://github.com/tsutterley/gravity-toolkit/pull/33>`_
* ``docs``: update documentation for GFZ solutions `(#33) <https://github.com/tsutterley/gravity-toolkit/pull/33>`_
* ``fix``: update grace input months for GFZ SLR `(#33) <https://github.com/tsutterley/gravity-toolkit/pull/33>`_
* ``feat``: added option for connection timeout to sync programs `(#34) <https://github.com/tsutterley/gravity-toolkit/pull/34>`_
* ``fix``: define int/float precision to prevent deprecation warning `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``feat``: use try/except for retrieving netrc credentials `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``feat``: add figshare secure FTP uploader to utilities `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: use cartopy no-binary in build `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``fix``: use first value in requirements in setup `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: brew install ``pkg-config`` `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: use older proj7 in brew install for cartopy `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: add LD and CPP flags for proj7 `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: add cython to installations `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: ``ACCEPT_USE_OF_DEPRECATED_PROJ_API_H`` `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``ci``: set pkg-config path `(#35) <https://github.com/tsutterley/gravity-toolkit/pull/35>`_
* ``refactor``: switch from parameter files to argparse arguments `(#36) <https://github.com/tsutterley/gravity-toolkit/pull/36>`_
* ``fix``: degree spacing in spatial programs `(#36) <https://github.com/tsutterley/gravity-toolkit/pull/36>`_
* ``fix``: cycles in regression program `(#38) <https://github.com/tsutterley/gravity-toolkit/pull/38>`_
* ``fix``: documentation for spatial programs `(#38) <https://github.com/tsutterley/gravity-toolkit/pull/38>`_
* ``refactor``: simplified file exports using wrappers in harmonics `(#39) <https://github.com/tsutterley/gravity-toolkit/pull/39>`_
* ``fix``: gfc format in ``from_file`` wrapper in harmonics `(#39) <https://github.com/tsutterley/gravity-toolkit/pull/39>`_
* ``fix``: format for mean files `(#40) <https://github.com/tsutterley/gravity-toolkit/pull/40>`_
* ``docs``: clenshaw summation citations `(#40) <https://github.com/tsutterley/gravity-toolkit/pull/40>`_
* ``feat``: Add 3-hour AOD interval for RL06  `#37 <https://github.com/tsutterley/gravity-toolkit/issues/37>`_ `(#41) <https://github.com/tsutterley/gravity-toolkit/pull/41>`_
* ``feat``: release monthly dealiasing (for CSR GAA etc) `(#41) <https://github.com/tsutterley/gravity-toolkit/pull/41>`_
* ``fix``: inputs to AOD-corrected SLR geocenter coefficients `(#41) <https://github.com/tsutterley/gravity-toolkit/pull/41>`_
* ``feat``: output index file for monthly dealiasing SHM files `(#42) <https://github.com/tsutterley/gravity-toolkit/pull/42>`_
* ``feat``: added check if needing to interpolate love numbers `(#42) <https://github.com/tsutterley/gravity-toolkit/pull/42>`_
* ``feat``: added path to default land-sea mask for mass redistribution `(#42) <https://github.com/tsutterley/gravity-toolkit/pull/42>`_
* ``feat``: added option to output mean harmonics in gfc format `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``refactor``: rename monthly mean dealiasing program `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``feat``: output uncalibrated spherical harmonic errors (eclm and eslm) `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``fix``: remove choices for argparse processing centers `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``fix``: remove defaults in monthly dealiasing `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``refactor``: no default processing center `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``fix``: require processing center argument `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``feat``: add months option to gfz dealiasing sync `(#43) <https://github.com/tsutterley/gravity-toolkit/pull/43>`_
* ``feat``: add harmonic resolution calculator `(#44) <https://github.com/tsutterley/gravity-toolkit/pull/44>`_

.. __: https://github.com/tsutterley/gravity-toolkit/releases/tag/v1.0.2.0
