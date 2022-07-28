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

.. autoclass:: gravity_toolkit.harmonics
   :members:
