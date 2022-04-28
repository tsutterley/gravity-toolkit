=====================
mascon_reconstruct.py
=====================

- Calculates the equivalent spherical harmonics from a mascon time series

Calling Sequence
################

.. code-block:: bash

    python mascon_reconstruct.py @mascon_reconstruct

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/mascon_reconstruct.py

Command Line Options
####################

- ``-O X``, ``--output-directory X``: output directory for mascon files
- ``-p X``, ``--product X``: GRACE/GRACE-FO Level-2 data product
- ``-S X``, ``--start X``: starting GRACE/GRACE-FO month
- ``-E X``, ``--end X``: ending GRACE/GRACE-FO month
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``-R X``, ``--radius X``: Gaussian smoothing radius (km)
- ``-d``, ``--destripe``: use decorrelation filter (destriping filter)
- ``-n X``, ``--love X``: Load Love numbers dataset

     * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
     * ``1``: Gegout (2005) values from PREM [Gegout2010]_
     * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``--reference X``: Reference frame for load love numbers

     * ``'CF'``: Center of Surface Figure (default)
     * ``'CM'``: Center of Mass of Earth System
     * ``'CE'``: Center of Mass of Solid Earth
- ``-F X``, ``--format X``: input data format for auxiliary files

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
- ``-G X``, ``--gia X``: GIA model type to read

    * ``'IJ05-R2'``: `Ivins R2 GIA Models <https://doi.org/10.1002/jgrb.50208>`_
    * ``'W12a'``: `Whitehouse GIA Models <https://doi.org/10.1111/j.1365-246X.2012.05557.x>`_
    * ``'SM09'``: `Simpson/Milne GIA Models <https://doi.org/10.1029/2010JB007776>`_
    * ``'ICE6G'``: `ICE-6G GIA Models <https://doi.org/10.1002/2014JB011176>`_
    * ``'Wu10'``: `Wu (2010) GIA Correction <https://doi.org/10.1038/ngeo938>`_
    * ``'AW13-ICE6G'``: `Geruo A ICE-6G GIA Models <https://doi.org/10.1093/gji/ggs030>`_
    * ``'AW13-IJ05'``: `Geruo A IJ05-R2 GIA Models <https://doi.org/10.1093/gji/ggs030>`_
    * ``'Caron'``: `Caron JPL GIA Assimilation <https://doi.org/10.1002/2017GL076644>`_
    * ``'ICE6G-D'``: `ICE-6G Version-D GIA Models <https://doi.org/10.1002/2016JB013844>`_
    * ``'ascii'``: reformatted GIA in ascii format
    * ``'netCDF4'``: reformatted GIA in netCDF4 format
    * ``'HDF5'``: reformatted GIA in HDF5 format
- ``--gia-file X``: GIA file to read
- ``--atm-correction``: Apply atmospheric jump correction coefficients
- ``--mask X``: Land-sea mask for redistributing mascon mass and land water flux
- ``--mascon-file X``: index file of mascons spherical harmonics
- ``--redistribute-mascons``: redistribute mascon mass over the ocean
- ``--fit-method X``: method for fitting sensitivity kernel to harmonics

    * ``1``: mass coefficients
    * ``2``: geoid coefficients
- ``--reconstruct-file X``: reconstructed mascon time series file
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created
