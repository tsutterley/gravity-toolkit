============
harmonics.py
============

Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data

 - Can read ascii, netCDF4, HDF5 files
 - Can read from an index of the above file types
 - Can merge a list of harmonics objects into a single object
 - Can subset to a list of GRACE/GRACE-FO months
 - Can calculate the mean field of a harmonics object
 - Can filter harmonics for correlated "striping" errors
 - Can output harmonics objects to ascii, netCDF4 or HDF5 files

Calling Sequence
================

Reading a netCDF4 file

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_HDF5(path_to_HDF5_file)

Reading an index file of netCDF4 files and then outputting as a single file

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_index(path_to_index_file,'netCDF4')
    Ylms.to_netCDF4(path_to_netCDF4_file)

Reading an index file of HDF5 files and subsetting to specific months

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_index(path_to_index_file,'HDF5').subset(months)

Reading an index file of ascii files and truncating to a new degree and order

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_index(path_to_index_file,'ascii').truncate(lmax)

Converting a dictionary object to a harmonics object and removing the mean field

.. code-block:: python

    from gravity_toolkit.harmonics import harmonics
    Ylms = harmonics().from_dict(Ylms_dict)
    Ylms.mean(apply=True)


`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/harmonics.py

General Attributes and Methods
==============================

.. class:: harmonics(object)


    .. attribute:: object.lmax

        maximum degree of the spherical harmonic field


    .. attribute:: object.mmax

        maximum order of the spherical harmonic field


    .. attribute:: object.clm

        cosine spherical harmonics


    .. attribute:: object.slm

        sine spherical harmonics


    .. attribute:: object.time

        time variable of the spherical harmonics


    .. attribute:: object.month

        GRACE/GRACE-FO months variable of the spherical harmonics


    .. attribute:: object.shape

        dimensions of harmonics object


    .. attribute:: object.ndim

        number of dimensions of harmonics object


    .. method:: object.case_insensitive_filename(filename)

        Searches a directory for a filename without case dependence


    .. method:: object.from_ascii(filename, date=True, compression=None, verbose=False)

        Read a harmonics object from an ascii file

        Inputs: full path of input ascii file

        Options:
            `date` ascii file contains date information

            `compression` ascii file is compressed or streamed from memory

            `verbose` print ascii file information


    .. method:: object.from_netCDF4(filename, date=True, compression=None, verbose=False)

        Read a harmonics object from a netCDF4 file

        Inputs: full path of input netCDF4 file

        Options:
            `date` netCDF4 file contains date information

            `compression` netCDF4 file is compressed or streamed from memory

            `verbose` print netCDF4 file information


    .. method:: object.from_HDF5(filename, date=True, compression=None, verbose=False)

        Read a harmonics object from a HDF5 file

        Inputs: full path of input HDF5 file

        Options:
            `date` HDF5 file contains date information

            `compression` HDF5 file is compressed or streamed from memory

            `verbose` print HDF5 file information


    .. method:: object.from_gfc(filename, verbose=False)

        Read a harmonics object from a gfc gravity model file from the `GFZ ICGEM`__.

        Inputs: full path of input gfc file

        Options:
            `verbose` print gfc file information

    .. __: http://icgem.gfz-potsdam.de/


    .. method:: object.from_index(filename, format=None, date=True, sort=True)

        Read a harmonics object from an index of ascii, netCDF4 or HDF5 files

        Inputs: full path of index file to be read into a harmonics object

        Options:
            format of files in index (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            ascii, netCDF4, or HDF5 files contain date information

            sort harmonics objects by date information


    .. method:: object.from_list(object_list, date=True, sort=True, clear=False)

        Build a sorted harmonics object from a list of other harmonics objects

        Inputs: list of harmonics object to be merged

        Options:
            harmonics objects contain date information

            sort harmonics objects by date information

            clear the list of objects from memory


    .. method:: object.from_file(filename, format=None, date=True, **kwargs)

        Read a harmonics object from a specified format

        Inputs: full path of input file

        Options:
            file format (``'ascii'``, ``'netCDF4'``, ``'HDF5'`` or ``'gfc'``)

            file contains date information

            keyword arguments for input readers


    .. method:: object.from_dict(dict_object)

        Convert a dict object to a harmonics object

        Inputs: dictionary object to be converted


    .. method:: object.to_ascii(filename, date=True)

        Write a harmonics object to ascii file

        Inputs: full path of output ascii file

        Options: harmonics objects contain date information


    .. method:: object.to_netCDF4(filename, date=True)

        Write a harmonics object to netCDF4 file

        Inputs: full path of output netCDF4 file

        Options: harmonics objects contain date information


    .. method:: object.to_HDF5(filename, date=True)

        Write a harmonics object to HDF5 file

        Inputs: full path of output HDF5 file

        Options: harmonics objects contain date information


    .. method:: object.to_index(filename, file_list, format=None, date=True, **kwargs)

        Write a harmonics object to index of ascii, netCDF4 or HDF5 files

        Inputs:
            full path of output file

            list of filenames for each output file

        Options:

            format of files in index (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            harmonics object contains date information

            keyword arguments for output writers


    .. method:: object.to_file(filename, format=None, date=True, **kwargs)

        Write a harmonics object to a specified format

        Inputs: full path of output file

        Options:

            file format (``'ascii'``, ``'netCDF4'`` or ``'HDF5'``)

            harmonics object contains date information

            keyword arguments for output writers


    .. method:: object.to_masked_array()

        Convert a harmonics object to a masked numpy array


    .. method:: object.update_dimensions()

        Update the dimensions of the harmonics object


    .. method:: object.add(temp)

        Add two harmonics objects

        Inputs: harmonic object to be added


    .. method:: object.subtract(temp)

        Subtract one harmonics object from another

        Inputs: harmonic object to be subtracted


    .. method:: object.multiply(temp)

        Multiply two harmonics objects

        Inputs: harmonic object to be multiplied


    .. method:: object.divide(temp)

        Divide one harmonics object from another

        Inputs: harmonic object to be divided


    .. method:: object.copy()

        Copy a harmonics object to a new harmonics object


    .. method:: object.zeros_like()

        Create a harmonics object using the dimensions of another


    .. method:: object.expand_dims()

        Add a singleton dimension to a harmonics object if non-existent


    .. method:: object.squeeze()

        Remove singleton dimensions from a harmonics object


    .. method:: object.flatten(date=True)

        Flatten harmonics matrices into arrays

        Options: harmonics objects contain date information


    .. method:: expand.expand(date=True)

        Expand flattened harmonics into matrices

        Options: harmonics objects contain date information


    .. method:: object.index(indice, date=True)

        Subset a harmonics object to specific index

        Inputs: `indice` in matrix to subset

        Options: harmonics objects contain date information


    .. method:: object.subset(months)

        Subset a harmonics object to specific GRACE/GRACE-FO months

        Inputs: GRACE/GRACE-FO months


    .. method:: object.truncate(lmax, lmin=0, mmax=None)

        Truncate a harmonics object to a new degree and order

        Inputs: `lmax` maximum degree of spherical harmonics

        Options:
            `lmin` minimum degree of spherical harmonics

            `mmax` maximum order of spherical harmonics


    .. method:: object.mean(apply=False, indices=Ellipsis)

        Compute mean gravitational field from the harmonics object

        Option:
            `apply` to remove the mean field from the input harmonics

            `indices` of input harmonics object to compute mean

    .. method:: object.scale(var)

        Multiply a harmonics object by a constant

        Inputs: scalar value to which the harmonics object will be multiplied


    .. method:: object.power(pow)

        Raise a harmonics object to a power

        Inputs: power to which the harmonics object will be raised


    .. method:: object.convolve(var)

        Convolve spherical harmonics with a degree-dependent array

        Inputs: degree dependent array for convolution


    .. method:: object.destripe()

        Filters spherical harmonic coefficients for correlated "striping" errors following `Swenson and Wahr (2006)`__.

    .. __: https://doi.org/10.1029/2005GL025285


    .. method:: object.amplitude(mmax=None)

        Calculates the degree amplitude of a harmonics object

        Options: mmax maximum order of spherical harmonics
