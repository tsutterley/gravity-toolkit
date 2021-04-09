===================
scale_grace_maps.py
===================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports scaled spatial fields
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms [Jekeli1981]_ [Swenson2006]_
- Converts data to centimeters water equivalent, performs a spherical harmonic summation to convert to the spatial domain [Wahr1998]_
- Scales the spatial fields following [Landerer2012]_
- Calculates the scaled spatial error field following [Wahr2006]_

Calling Sequence
################

.. code-block:: bash

     python scale_grace_maps.py --mode 0o775 parameter_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/scale_grace_maps.py

Inputs
######

- parameter file containing specific variables for the analysis

Command Line Options
####################

- ``-P X``, ``--np X``: Run in parallel with X number of processes
- ``-D X``, ``--directory X``: Working data directory
- ``-n X``, ``--love X``: Load Love numbers dataset

     * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
     * ``1``: Gegout (2005) values from PREM [Gegout2010]_
     * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``-r X``, ``--reference X``: Reference frame for load love numbers

     * ``'CF'``: Center of Surface Figure (default)
     * ``'CM'``: Center of Mass of Earth System
     * ``'CE'``: Center of Mass of Solid Earth
- ``-V``, ``--verbose``: verbose output of processing run
- ``-l``, ``--log``: output log file for each job
- ``-M X``, ``--mode X``: permissions mode of output files

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Landerer2012] F. W. Landerer and S. C. Swenson, "Accuracy of scaled GRACE terrestrial water storage estimates", *Water Resources Research*, 48(W04531), (2012). `doi: 10.1029/2011WR011453 <https://doi.org/10.1029/2011WR011453>`_

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wahr2006] J. Wahr, S. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates", Geophysical Research Letters, 33(L06401), (2006). `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
