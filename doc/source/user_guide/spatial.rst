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

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/spatial.py

General Attributes and Methods
==============================

.. class:: spatial(object)


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

        spatial grid bounds ``[minimum longitude, maximum longitude, minimum latitude, maximum latitude]``


    .. attribute:: object.spacing

        grid step size ``[longitude,latitude]``


    .. attribute:: object.shape

        grid dimensions


    .. attribute:: object.ndim

        number of grid dimensions


    .. method:: object.case_insensitive_filename(filename)

        Searches a directory for a filename without case dependence


    .. method:: object.from_ascii(filename, date=True, compression=None, verbose=False, columns=['lon','lat','data','time'], header=0)

        Read a spatial object from an ascii file

        Arguments:

            full path of input ascii file

        Keyword arguments:

            ``date`` ascii file contains date information

            ``compression`` ascii file is compressed or streamed from memory

            ``verbose`` print ascii filename

            ``columns`` variable names for each column

            ``header`` rows of header lines to skip


    .. method:: object.from_netCDF4(filename, date=True, compression=None, verbose=False, varname='z', lonname='lon', latname='lat', timename='time')

        Read a spatial object from a netCDF4 file

        Arguments:

            full path of input netCDF4 file

        Keyword arguments:

            ``date`` netCDF4 file contains date information

            ``compression`` netCDF4 file is compressed or streamed from memory

            ``verbose`` print netCDF4 file information

            ``varname`` input variable name in netCDF4 file

            ``lonname`` input longitude variable name in netCDF4 file

            ``latname`` input latitude variable name in netCDF4 file

            ``timename`` input time variable name in netCDF4 file


    .. method:: object.from_HDF5(filename, date=True, compression=None, verbose=False, varname='z', lonname='lon', latname='lat', timename='time')

        Read a spatial object from a HDF5 file

        Arguments:

            full path of input HDF5 file

        Keyword arguments:

            ``date`` HDF5 file contains date information

            ``compression`` HDF5 file is compressed or streamed from memory

            ``verbose`` print HDF5 file information

            ``varname`` input variable name in HDF5 file

            ``lonname`` input longitude variable name in HDF5 file

            ``latname`` input latitude variable name in HDF5 file

            ``timename`` input time variable name in HDF5 file


    .. method:: object.from_index(filename, format=None, date=True, sort=True)

        Read a spatial object from an index of ascii, netCDF4 or HDF5 files

        Arguments:

            full path of index file to be read into a spatial object

        Keyword arguments:

            format of files in index (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            ascii, netCDF4, or HDF5 files contain date information

            sort spatial objects by date information


    .. method:: object.from_list(object_list, date=True, sort=True, clear=False)

        Build a sorted spatial object from a list of other spatial objects

        Arguments:

            list of spatial object to be merged

        Keyword arguments:

            spatial objects contain date information

            sort spatial objects by date information

            clear the list of objects from memory


    .. method:: object.from_file(filename, format=None, date=True, **kwargs)

        Read a spatial object from a specified format

        Arguments:

            full path of input file

        Keyword arguments:

            file format (``'ascii'``, ``'netCDF4'``, ``'HDF5'``)

            file contains date information

            keyword arguments for input readers


    .. method:: object.from_dict(dict_object)

        Convert a dict object to a spatial object

        Arguments: dictionary object to be converted


    .. method:: object.to_ascii(filename, date=True, verbose=False)

        Write a spatial object to ascii file

        Arguments:

            full path of output ascii file

        Keyword arguments:

            ``date`` spatial objects contain date information

            ``verbose`` print ascii file name


    .. method:: object.to_netCDF4(filename, date=True, varname='z', units=None, longname=None, title=None, verbose=False)

        Write a spatial object to netCDF4 file

        Arguments:

            full path of output netCDF4 file

        Keyword arguments:

            ``date`` spatial objects contain date information

            ``varname`` output variable name in netCDF4 file

            ``units`` output variable units in netCDF4 file

            ``longname`` output variable unit longname in netCDF4 file

            ``title`` output netCDF4 file title

            ``verbose`` print netCDF4 file information


    .. method:: object.to_HDF5(filename, date=True, varname='z', units=None, longname=None, title=None, verbose=False)

        Write a spatial object to HDF5 file

        Arguments:

            full path of output HDF5 file

        Keyword arguments:

            ``date`` spatial objects contain date information

            ``varname`` output variable name in HDF5 file

            ``units`` output variable units in HDF5 file

            ``longname`` output variable unit longname in HDF5 file

            ``title`` output HDF5 file title

            ``verbose`` print HDF5 file information


    .. method:: object.to_index(filename, file_list, format=None, date=True, **kwargs)

        Write a spatial object to index of ascii, netCDF4 or HDF5 files

        Arguments:

            full path of output HDF5 file

            list of filenames for each output file

        Keyword arguments:

            format of files in index (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            spatial object contains date information

            keyword arguments for output writers


    .. method:: object.to_file(filename, format=None, date=True, **kwargs)

        Write a spatial object to a specified format

        Arguments:

            full path of output file

        Keyword arguments:

            file format (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            spatial object contains date information

            keyword arguments for output writers


    .. method:: object.to_masked_array()

        Convert a spatial object to a masked numpy array


    .. method:: object.update_spacing()

        Calculate the step size of spatial object


    .. method:: object.update_extents()

        Calculate the bounds of spatial object


    .. method:: object.update_dimensions()

        Update the dimensions of the spatial object


    .. method:: object.update_mask()

        Update the mask of the spatial object


    .. method:: object.copy()

        Copy a spatial object to a new spatial object


    .. method:: object.zeros_like()

        Create a spatial object using the dimensions of another


    .. method:: object.expand_dims()

        Add a singleton dimension to a spatial object if non-existent


    .. method:: object.squeeze()

        Remove singleton dimensions from a spatial object


    .. method:: object.index(indice, date=True)

        Subset a spatial object to specific index

        Arguments:

            ``indice`` in matrix to subset

        Keyword arguments:

            spatial objects contain date information


    .. method:: object.subset(months)

        Subset a spatial object to specific GRACE/GRACE-FO months

        Arguments: GRACE/GRACE-FO months


    .. method:: object.offset(var)

        Offset a spatial object by a constant

        Arguments: scalar value to which the spatial object will be offset


    .. method:: object.scale(var)

        Multiply a spatial object by a constant

        Arguments: scalar value to which the spatial object will be multiplied


    .. method:: object.kfactor(var)

        Calculate the scaling factor and scaling factor errors from two spatial objects following `Landerer and Swenson (2012) <https://doi.org/10.1029/2011WR011453>`_

        Arguments:

            spatial object to used for scaling

        Returns:

            scaling factor and scaling factor error


    .. method:: object.mean(apply=False, indices=Ellipsis)

        Compute mean spatial field and remove from data if specified

        Keyword arguments:

            ``apply`` to remove the mean field from the input spatial

            ``indices`` of spatial object to compute mean


    .. method:: object.reverse(axis=0)

        Reverse the order of data and dimensions along an axis

        Keyword arguments:

            ``axis`` to reorder


    .. method:: object.transpose(axes=None)

        Reverse or permute the axes of a spatial object

        Keyword arguments:

            order of the output axes


    .. method:: object.sum(power=1)

        Compute summation of spatial field

        Keyword arguments:

            apply a power before calculating summation


    .. method:: object.power(pow)

        Raise a spatial object to a power

        Arguments:

            power to which the spatial object will be raised


    .. method:: object.max()

        Compute maximum value of spatial field


    .. method:: object.min()

        Compute minimum value of spatial field


    .. method:: object.replace_invalid(fill_value, mask=None)

        Replace the masked values with a new fill_value

        Keyword arguments:

            ``mask`` to update the current mask


    .. method:: object.replace_masked()

        Replace the masked values with fill_value
