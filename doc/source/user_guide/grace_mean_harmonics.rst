=======================
grace_mean_harmonics.py
=======================

- Calculates the temporal mean of the GRACE/GRACE-FO spherical harmonics for a specified date range
- Used to estimate the static gravitational field over a given date rage

Calling Sequence
################
.. code-block:: bash

     python grace_mean_harmonics.py @default_arguments_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_mean_harmonics.py

Command Line Options
####################

- ``-D X``, ``--directory X``: Working data directory
- ``-c X``, ``--center X``: GRACE/GRACE-FO processing center
- ``-r X``, ``--release X``: GRACE/GRACE-FO data release
- ``-p X``, ``--product X``: GRACE/GRACE-FO Level-2 data product
- ``-S X``, ``--start X``: starting GRACE/GRACE-FO month
- ``-E X``, ``--end X``: ending GRACE/GRACE-FO month
- ``-N X``, ``--missing X``: Missing GRACE/GRACE-FO months
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``--atm-correction``: Apply atmospheric jump correction coefficients
- ``--pole-tide``: Correct for pole tide drift
- ``--geocenter X``: Update Degree 1 coefficients with SLR or derived values

    * ``None``
    * ``'Tellus'``: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
    * ``'SLR'``: satellite laser ranging coefficients from CSR
    * ``'SLF'``: Sutterley and Velicogna coefficients, Remote Sensing (2019)
    * ``'Swenson'``: GRACE-derived coefficients from Sean Swenson
    * ``'GFZ'``: GRACE/SLR derived coefficients from GFZ GravIS
- ``--geocenter-file X``: Specific geocenter file if not default
- ``--interpolate-geocenter``: Least-squares model missing Degree 1 coefficients
- ``--slr-c20 X``: Replace *C*\ :sub:`20` coefficients with SLR values

    * ``None``: use original values
    * ``'CSR'``: use values from CSR (TN-07, TN-09, TN-11)
    * ``'GFZ'``: use values from GFZ
    * ``'GSFC'``: use values from GSFC (TN-14)
- ``--slr-21 X``: Replace *C*\ :sub:`21` and *S*\ :sub:`21` coefficients with SLR values

    * ``None``: use original values
    * ``'CSR'``: use values from CSR
    * ``'GFZ'``: use values from GFZ GravIS
- ``--slr-22 X``: Replace *C*\ :sub:`22` and *S*\ :sub:`22` coefficients with SLR values

    * ``None``: use original values
    * ``'CSR'``: use values from CSR
- ``--slr-c30 X``: Replace *C*\ :sub:`30` coefficients with SLR values

    * ``None``: use original values
    * ``'CSR'``: use values from CSR (5x5 with 6,1)
    * ``'GFZ'``: use values from GFZ GravIS
    * ``'GSFC'``: use values from GSFC (TN-14)
    * ``'LARES'``: use filtered values from CSR
- ``--slr-c50 X``: Replace *C*\ :sub:`50` coefficients with SLR values

    * ``None``: use original values
    * ``'CSR'``: use values from CSR (5x5 with 6,1)
    * ``'GSFC'``: use values from GSFC
    * ``'LARES'``: use filtered values from CSR
- ``--mean-file X``: GRACE/GRACE-FO mean file to remove from the harmonic data
- ``--mean-format X``: Output data format for GRACE/GRACE-FO mean file

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
    * ``'gfc'``
- ``--log``: Output log file for job
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created
