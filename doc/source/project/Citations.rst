====================
Citation Information
====================

References
##########

The programs included in this software have contributed
most recently to the following work:

    I. Velicogna, Y. Mohajerani, G. A, F. Landerer, J. Mouginot, B. No\ |euml|\ l,
    E. Rignot, T. C. Sutterley, M. van den Broeke, J. M. van Wessem, and D. Wiese,
    "Continuity of ice sheet mass loss in Greenland and Antarctica from the GRACE
    and GRACE Follow-On missions", *Geophysical Research Letters*, 47,
    (2020). `doi: 10.1029/2020GL087291 <https://doi.org/10.1029/2020GL087291>`_

    T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self-Consistent Ice Mass Balance
    and Regional Sea Level From Time-Variable Gravity", *Earth and Space Science*, 7,
    (2020). `doi: 10.1029/2019EA000860 <https://doi.org/10.1029/2019EA000860>`_

    T. C. Sutterley and I. Velicogna, "Improved estimates of geocenter variability
    from time-variable gravity and ocean model outputs", *Remote Sensing*, 11(18),
    2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

Some of the text within this documentation are modifications from these
original works, which is allowed under the
`publisher's permissions policies for authors <https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Policies/Permission-policy>`_
or through the `Creative Commons licensing <http://creativecommons.org/licenses/by/4.0/>`_ of the work.

.. admonition:: Please consider citing our library

    T. C. Sutterley, et al., "gravity-toolkit: Python tools for obtaining and
    working with data from the GRACE and GRACE Follow-On missions", (2020).
    `doi: 10.5281/zenodo.5156864 <https://doi.org/10.5281/zenodo.5156864>`_

Dependencies
############

This software is also dependent on other commonly used Python packages:

- `boto3: Amazon Web Services (AWS) SDK for Python <https://boto3.amazonaws.com/v1/documentation/api/latest/index.html>`_
- `future: Compatibility layer between Python 2 and Python 3 <https://python-future.org/>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_
- `matplotlib: Python 2D plotting library <https://matplotlib.org/>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/>`_
- `numpy: Scientific Computing Tools For Python <https://numpy.org>`_
- `platformdirs: Python module for determining platform-specific directories <https://pypi.org/project/platformdirs/>`_
- `python-dateutil: powerful extensions to datetime <https://dateutil.readthedocs.io/en/stable/>`_
- `PyYAML: YAML parser and emitter for Python <https://github.com/yaml/pyyaml>`_
- `scipy: Scientific Tools for Python <https://docs.scipy.org/doc/>`_

Optional Dependencies
---------------------

- `cartopy: Python package designed for geospatial data processing <https://scitools.org.uk/cartopy/docs/latest/>`_
- `dask: Parallel computing with task scheduling <https://www.dask.org/>`_
- `geoid-toolkit: Python utilities for calculating geoid heights from static gravity field coefficients <https://github.com/tsutterley/geoid-toolkit/>`_
- `gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL) <https://pypi.python.org/pypi/GDAL>`_
- `h5py: Python interface for Hierarchal Data Format 5 (HDF5) <https://www.h5py.org/>`_
- `ipywidgets: interactive HTML widgets for Jupyter notebooks and IPython <https://ipywidgets.readthedocs.io/en/latest/>`_
- `obstore: Simple, high-throughput Python interface for object storage <https://developmentseed.org/obstore>`_
- `pyarrow: Apache Arrow Python bindings <https://arrow.apache.org/docs/python/>`_
- `pyshp: Python read/write support for ESRI Shapefile format <https://github.com/GeospatialPython/pyshp>`_
- `s3fs: Pythonic file interface to S3 built on top of botocore <https://s3fs.readthedocs.io/en/latest/>`_
- `shapely: PostGIS-ish operations outside a database context for Python <http://toblerity.org/shapely/index.html>`_
- `tkinter: Python interface to the Tcl/Tk GUI toolkit <https://docs.python.org/3/library/tkinter.html>`_
- `zarr: Chunked, compressed, N-dimensional arrays in Python <https://zarr.readthedocs.io/en/stable/>`_

Disclaimer
##########

This work is currently supported by the NASA GRACE-FO Science Team (Grant Number 80NSSC24K1153).
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR),
the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.

.. caution::
    This software is provided here for your convenience but *with no guarantees whatsoever*.

This product includes software developed at:

- University of California, Irvine, Department of Earth System Science
- Jet Propulsion Laboratory, California Institute of Technology
- National Aeronautics and Space Administration, Goddard Space Flight Center
- University of Washington, Applied Physics Laboratory, Polar Science Center

Acknowledgements
################

Much of this software stems from the work of `John Wahr (1951--2015) <http://www.johnwahr.com/>`_,
who made immeasurable contributions towards time-variable gravity research and was an inspiration.

.. |euml|    unicode:: U+00EB .. LATIN SMALL LETTER E WITH DIAERESIS
