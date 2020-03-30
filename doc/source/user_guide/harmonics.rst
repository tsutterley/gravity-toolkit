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
 - Can output harmonics objects to netCDF4 or HDF5 files

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


General Attributes and Methods
==============================

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


.. attribute:: object.from_ascii(filename)

    read a harmonics object from an ascii file

    Inputs: full path of input ascii file


.. attribute:: object.from_netCDF4(filename)

    read a harmonics object from a netCDF4 file

    Inputs: full path of input netCDF4 file


.. attribute:: object.from_HDF5(filename)

    read a harmonics object from a HDF5 file

    Inputs: full path of input HDF5 file


.. attribute:: object.from_index(filename, format=None)

    read a harmonics object from an index of ascii, netCDF4 or HDF5 files

    Inputs: full path of index file to be read into a harmonics object

    Options: format of files in index (ascii, netCDF4 or HDF5)


.. attribute:: object.from_list(object_list)

    build a sorted harmonics object from a list of other harmonics objects

    Inputs: list of harmonics object to be merged


.. attribute:: object.from_dict(dict_object)

    convert a dict object to a harmonics object

    Inputs: dictionary object to be converted


.. attribute:: object.to_netCDF4(filename)

    write a harmonics object to netCDF4 file

    Inputs: full path of output netCDF4 file


.. attribute:: object.to_HDF5(filename)

    write a harmonics object to HDF5 file

    Inputs: full path of output HDF5 file


.. attribute:: object.add(temp)

    add two harmonics objects

    Inputs: harmonic object to be added


.. attribute:: object.subtract(temp)

    subtract one harmonics object from another

    Inputs: harmonic object to be subtracted


.. attribute:: object.index(indice)

    subset a harmonics object to specific index

    Inputs: `indice` in matrix to subset

.. attribute:: object.subset(months)

    subset a harmonics object to specific GRACE/GRACE-FO months

    Inputs: GRACE/GRACE-FO months


.. attribute:: object.truncate(lmax, mmax=None)

    truncate a harmonics object to a new degree and order

    Inputs: `lmax` maximum degree of spherical harmonics

    Options: `mmax` maximum order of spherical harmonics


.. attribute:: object.mean(apply=False)

    Compute mean gravitational field from the harmonics object

    Option: `apply` to remove the mean field from the input harmonics


.. attribute:: object.convolve(var)

    Convolve spherical harmonics with a degree-dependent array

    Inputs: degree dependent array for convolution


.. attribute:: object.destripe(self)

    Filters spherical harmonic coefficients for correlated "striping" errors following `Swenson and Wahr (2006)`__.

.. __: https://doi.org/10.1029/2005GL025285
