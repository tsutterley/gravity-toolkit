==========
spatial.py
==========

Spatial data class for reading, writing and processing spatial data

 - Can read ascii, netCDF4, HDF5 files
 - Can read from an index of the above file types
 - Can merge a list of spatial objects into a single object
 - Can subset to a list of GRACE/GRACE-FO months
 - Can calculate the mean field of a spatial object
 - Can output spatial objects to ascii, netCDF4 or HDF5 files

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

.. __: https://github.com/tsutterley/read-GRACE-spatial/blob/master/gravity_toolkit/spatial.py

General Attributes and Methods
==============================


.. attribute:: object.data

    spatial grid data


.. attribute:: object.mask

    spatial grid mask


.. attribute:: object.lat

    grid latitudes


.. attribute:: object.lon

    grid longitudes


.. attribute:: object.time

    time variable of the spatial data


.. attribute:: object.month

    GRACE/GRACE-FO months variable of the spatial data


.. attribute:: object.fill_value

    invalid value for spatial grid data


.. attribute:: object.extent

    invalid value for spatial grid data


.. attribute:: object.shape

    grid dimensions


.. attribute:: object.ndim

    number of grid dimensions


.. attribute:: object.spacing

    grid step size [longitudes,latitudes]


.. method:: object.from_ascii(filename, date=True)

    Read a spatial object from an ascii file

    Inputs: full path of input ascii file

    Options: ascii file contains date information


.. method:: object.from_netCDF4(filename, date=True)

    Read a spatial object from a netCDF4 file

    Inputs: full path of input netCDF4 file

    Options: netCDF4 file contains date information


.. method:: object.from_HDF5(filename, date=True)

    Read a spatial object from a HDF5 file

    Inputs: full path of input HDF5 file

    Options: HDF5 file contains date information


.. method:: object.from_index(filename, format=None, date=True, sort=True)

    Read a spatial object from an index of ascii, netCDF4 or HDF5 files

    Inputs: full path of index file to be read into a spatial object

    Options:
        format of files in index (ascii, netCDF4 or HDF5)

        ascii, netCDF4, or HDF5 contains date information

        sort spatial objects by date information


.. method:: object.from_list(object_list, date=True, sort=True)

    Build a sorted spatial object from a list of other spatial objects

    Inputs: list of spatial object to be merged

    Options:
        spatial objects contain date information

        sort spatial objects by date information


.. method:: object.from_dict(dict_object)

    Convert a dict object to a spatial object

    Inputs: dictionary object to be converted


.. method:: object.to_ascii(filename)

    Write a spatial object to ascii file

    Inputs: full path of output ascii file

    Options: spatial objects contain date information


.. method:: object.to_netCDF4(filename)

    Write a spatial object to netCDF4 file

    Inputs: full path of output netCDF4 file

    Options: spatial objects contain date information


.. method:: object.to_HDF5(filename)

    Write a spatial object to HDF5 file

    Inputs: full path of output HDF5 file

    Options: spatial objects contain date information


.. method:: object.update_spacing()

    Calculate the step size of spatial object


.. method:: object.update_extent()

    Calculate the bounds of spatial object


.. method:: object.update_dimensions()

    Update the dimensions of the spatial object with new extents


.. method:: object.update_mask()

    Update the mask of the spatial object


.. method:: object.copy()

    Copy a spatial object to a new spatial object


.. method:: object.expand_dims()

    Add a singleton dimension to a spatial object if non-existent


.. method:: object.squeeze()

    Remove singleton dimensions from a spatial object


.. method:: object.index(indice, date=True)

    Subset a spatial object to specific index

    Inputs: `indice` in matrix to subset

    Options: spatial objects contain date information


.. method:: object.subset(months)

    Subset a spatial object to specific GRACE/GRACE-FO months

    Inputs: GRACE/GRACE-FO months


.. method:: object.mean(apply=False)

    Compute mean gravitational field from the spatial object

    Option: `apply` to remove the mean field from the input spatial
