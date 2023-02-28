=======
spatial
=======

Spatial data class for reading, writing and processing spatial data

 - Can read ascii, netCDF4, HDF5 files
 - Can read from an index of the above file types
 - Can merge a list of ``spatial`` objects into a single object
 - Can subset to a list of GRACE/GRACE-FO months
 - Can calculate the mean field of a ``spatial`` object
 - Can output ``spatial`` objects to ascii, netCDF4 or HDF5 files

Calling Sequence
================

Reading a netCDF4 file

.. code-block:: python

    from gravity_toolkit.spatial import spatial
    grid = spatial().from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    from gravity_toolkit.spatial import spatial
    grid = spatial().from_HDF5(path_to_HDF5_file)

Reading an index file of netCDF4 files and then outputting as a single file

.. code-block:: python

    from gravity_toolkit.spatial import spatial
    grid = spatial().from_index(path_to_index_file,'netCDF4')
    grid.to_netCDF4(path_to_netCDF4_file)

Reading an index file of HDF5 files and subsetting to specific months

.. code-block:: python

    from gravity_toolkit.spatial import spatial
    grid = spatial().from_index(path_to_index_file,'HDF5').subset(months)

Converting a dictionary object to a spatial object and removing the mean field

.. code-block:: python

    from gravity_toolkit.spatial import spatial
    grid = spatial().from_dict(grid_dict)
    grid.mean(apply=True)


`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/spatial.py

General Attributes and Methods
==============================

.. autoclass:: gravity_toolkit.spatial
   :members:
