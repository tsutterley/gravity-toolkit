==============
calc_mascon.py
==============

- Reads in GRACE/GRACE-FO spherical harmonic coefficients
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms [Jekeli1981]_ [Swenson2006]_
- Calculates a time-series of regional mass anomalies through a least-squares mascon procedure following [Tiwari2009]_ [Jacob2012]_
- Calculates the regional mascon errors following [Wahr2006]_

Calling Sequence
################

.. code-block:: bash

     python calc_mascon.py @default_arguments_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/calc_mascon.py

Command Line Options
####################

- ``-D X``, ``--directory X``: Working data directory
- ``-O X``, ``--output-directory X``: output directory for mascon files
- ``-c X``, ``--center X``: GRACE/GRACE-FO processing center
- ``-r X``, ``--release X``: GRACE/GRACE-FO data release
- ``-p X``, ``--product X``: GRACE/GRACE-FO Level-2 data product
- ``-S X``, ``--start X``: starting GRACE/GRACE-FO month
- ``-E X``, ``--end X``: ending GRACE/GRACE-FO month
- ``-N X``, ``--missing X``: Missing GRACE/GRACE-FO months
- ``--lmin X``: minimum spherical harmonic degree
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
    * ``'Caron'``: `Caron JPL GIA Assimilation <https://doi.org/10.1002/2017GL076644>`_
    * ``'ICE6G-D'``: `ICE-6G Version-D GIA Models <https://doi.org/10.1002/2016JB013844>`_
    * ``'ascii'``: reformatted GIA in ascii format
    * ``'netCDF4'``: reformatted GIA in netCDF4 format
    * ``'HDF5'``: reformatted GIA in HDF5 format
- ``--gia-file X``: GIA file to read
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
- ``--mean-format X``: Input data format for GRACE/GRACE-FO mean file

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
    * ``'gfc'``
- ``--mask X``: Land-sea mask for redistributing mascon mass and land water flux
- ``--mascon-file X``: index file of mascons spherical harmonics
- ``--mascon-format X``: Input data format for mascon files

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
- ``--redistribute-mascons``: redistribute mascon mass over the ocean
- ``--fit-method X``: method for fitting sensitivity kernel to harmonics

    * ``1``: mass coefficients
    * ``2``: geoid coefficients
- ``--remove-file X``: Monthly files to be removed from the GRACE/GRACE-FO data
- ``--remove-format X``: Input data format for files to be removed

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
    * ``'index-ascii'``
    * ``'index-netCDF4'``
    * ``'index-HDF5'``
- ``--redistribute-removed``: redistribute removed mass fields over the ocean
- ``--reconstruct-file X``: reconstructed mascon time series file to be removed
- ``--remove-reconstruct``: remove reconstructed mascon time series fields
- ``--log``: Output log file for job
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Jacob2012] T. Jacob, J. Wahr, W. T. Pfeffer, and S. Swenson, "Recent contributions of glaciers and ice caps to sea level rise", *Nature*, 482, 514--518, (2012). `doi: 10.1038/nature10847 <https://doi.org/10.1038/nature10847>`_

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Tiwari2009] V. M. Tiwari, J. Wahr, and S. Swenson, "Dwindling groundwater resources in northern India, from satellite gravity observations", *Geophysical Research Letters*, 36(L18401), (2009). `doi: 10.1029/2009GL039401 <https://doi.org/10.1029/2009GL039401>`_

.. [Wahr2006] J. Wahr, S. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates", Geophysical Research Letters, 33(L06401), (2006). `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
