====================
convert_harmonics.py
====================

- Converts a file from the spatial domain into the spherical harmonic domain [Wahr1998]_

Calling Sequence
################

.. code-block:: bash

     python convert_harmonics.py --format netCDF4 --lmax 60 --units 1 input_file output_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/convert_harmonics.py

Command Line Options
####################

- ``--help``: list the command line options
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``-n X``, ``--love X``: Load Love numbers dataset

     * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
     * ``1``: Gegout (2005) values from PREM [Gegout2010]_
     * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``-r X``, ``--reference X``: Reference frame for load love numbers

     * ``'CF'``: Center of Surface Figure (default)
     * ``'CM'``: Center of Mass of Earth System
     * ``'CE'``: Center of Mass of Solid Earth
- ``-U X``, ``--units X``: output units

     * ``1``: cm of water thickness (cm.w.e., g/cm\ :sup:`2`)
     * ``2``: Gigatonnes (Gt)
     * ``3``: mm of water thickness (kg/m\ :sup:`2`)
- ``-S X``, ``--spacing X``: spatial resolution of output data (dlon,dlat)
- ``-I X``, ``--interval X``: output grid interval

     * ``1``: (0:360, 90:-90)
     * ``2``: (degree spacing/2)
- ``--missing``: input spatial fields have missing values
- ``--fill-value X``: set fill_value for input spatial fields
- ``--header X``: number of header rows to skip in input ascii files
- ``--delimiter X``: delimiter in input ascii files
- ``-F X``, ``--format X``: input and output data format

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
