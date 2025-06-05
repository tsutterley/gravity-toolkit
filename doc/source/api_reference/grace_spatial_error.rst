======================
grace_spatial_error.py
======================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports spatial error field following :cite:t:`Wahr:2006bx`
- Filters and smooths data with specified processing algorithms :cite:p:`Jekeli:1981vj,Swenson:2006hu`
- Converts data to specified units and performs a spherical harmonic summation to convert error field to the spatial domain :cite:p:`Wahr:1998hy`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/grace_spatial_error.py

Calling Sequence
################

.. argparse::
    :filename: grace_spatial_error.py
    :func: arguments
    :prog: grace_spatial_error.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM :cite:p:`Han:1995go`
        * ``1``: Gegout (2005) values from PREM :cite:p:`Gegout:2010gc`
        * ``2``: Wang et al. (2012) values from PREM :cite:p:`Wang:2012gc`
        * ``3``: Wang et al. (2012) values from PREM with hard sediment :cite:p:`Wang:2012gc`
        * ``4``: Wang et al. (2012) values from PREM with soft sediment :cite:p:`Wang:2012gc`

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth

    --geocenter : @after
        * ``None``
        * ``'Tellus'``: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        * ``'SLR'``: satellite laser ranging coefficients from CSR
        * ``'UCI'``: Sutterley and Velicogna coefficients, Remote Sensing (2019)
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
        * ``'CSR'``: use values from CSR (5\ |times|\ 5 with 6,1)
        * ``'GFZ'``: use values from GFZ GravIS
        * ``'GSFC'``: use values from GSFC (TN-14)
        * ``'LARES'``: use filtered values from CSR

    --slr-c50 : @replace
        Replace *C*\ :sub:`50` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (\ |times|\ 5 with 6,1)
        * ``'GSFC'``: use values from GSFC
        * ``'LARES'``: use filtered values from CSR

    --units -U : @after
        * ``1``: cm of water thickness
        * ``2``: mm of geoid height
        * ``3``: mm of elastic crustal deformation
        * ``4``: |mu|\ Gal gravitational perturbation
        * ``5``: mbar equivalent surface pressure

    --interval : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)

.. |mu|      unicode:: U+03BC .. GREEK SMALL LETTER MU

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
