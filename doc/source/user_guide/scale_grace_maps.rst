===================
scale_grace_maps.py
===================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports scaled spatial fields
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms [Jekeli1981]_ [Swenson2006]_
- Converts data to centimeters water equivalent, performs a spherical harmonic summation to convert to the spatial domain [Wahr1998]_
- Scales the spatial fields following [Landerer2012]_
- Calculates the scaled spatial error field following [Wahr2006]_

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/scale_grace_maps.py

Calling Sequence
################

.. argparse::
    :filename: ../scripts/scale_grace_maps.py
    :func: arguments
    :prog: scale_grace_maps.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
        * ``1``: Gegout (2005) values from PREM [Gegout2010]_
        * ``2``: Wang et al. (2012) values from PREM [Wang2012]_

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth

    --gia -G : @after
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

    --geocenter : @after
        * ``None``
        * ``'Tellus'``: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        * ``'SLR'``: satellite laser ranging coefficients from CSR
        * ``'SLF'``: Sutterley and Velicogna coefficients, Remote Sensing (2019)
        * ``'Swenson'``: GRACE-derived coefficients from Sean Swenson
        * ``'GFZ'``: GRACE/SLR derived coefficients from GFZ GravIS

    --slr-c20 : @replace
        Replace *C*\ :sub:`20` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (TN-07, TN-09, TN-11)
        * ``'GFZ'``: use values from GFZ
        * ``'GSFC'``: use values from GSFC (TN-14)

    --slr-21 X : @replace
        Replace *C*\ :sub:`21` and *S*\ :sub:`21` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR
        * ``'GFZ'``: use values from GFZ GravIS
        * ``'GSFC'``: use values from GSFC

    --slr-22 : @replace
        Replace *C*\ :sub:`22` and *S*\ :sub:`22` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR
        * ``'GSFC'``: use values from GSFC

    --slr-c30 : @replace
        Replace *C*\ :sub:`30` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (5x5 with 6,1)
        * ``'GFZ'``: use values from GFZ GravIS
        * ``'GSFC'``: use values from GSFC (TN-14)
        * ``'LARES'``: use filtered values from CSR

    --slr-c50 : @replace
        Replace *C*\ :sub:`50` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (5x5 with 6,1)
        * ``'GSFC'``: use values from GSFC
        * ``'LARES'``: use filtered values from CSR

    --interval : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)

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
