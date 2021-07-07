=====================
regress_grace_maps.py
=====================

- Reads in GRACE/GRACE-FO spatial files and fits a regression model at each grid point

Calling Sequence
################

.. code-block:: bash

    python regress_grace_maps.py --order 1 --cycles=0.5,1.0 @default_arguments_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py

Command Line Options
####################

- ``-O X``, ``--output-directory X``: output directory for spatial files
- ``-P X``, ``--file-prefix X``: prefix string for input and output files
- ``-S X``, ``--start X``: starting GRACE/GRACE-FO month for time series regression
- ``-E X``, ``--end X``: ending GRACE/GRACE-FO month for time series regression
- ``-N X``, ``--missing X``: Missing GRACE/GRACE-FO months
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``-R X``, ``--radius X``: Gaussian smoothing radius (km)
- ``-d``, ``--destripe``: use decorrelation filter (destriping filter)
- ``-U X``, ``--units X``: output units
    * ``1``: cm of water thickness
    * ``2``: mm of geoid height
    * ``3``: mm of elastic crustal deformation [Davis 2004]
    * ``4``: microGal gravitational perturbation
    * ``5``: mbar equivalent surface pressure
- ``--spacing X``: spatial resolution of output data (dlon,dlat)
- ``--interval X``: output grid interval
    * ``1``: (0:360, 90:-90)
    * ``2``: (degree spacing/2)
    * ``3``: non-global grid (set with defined bounds)
- ``--bounds X``: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
- ``-F X``, ``--format X``: input/output data format

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
- ``--redistribute-removed``: redistribute removed mass fields over the ocean
- `--order X`: regression fit polynomial order
- `--cycles X`: regression fit cyclical terms as wavelength in decimal years
- ``--log``: Output log file for job
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created
