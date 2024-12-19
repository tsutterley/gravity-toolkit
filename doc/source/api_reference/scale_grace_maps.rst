===================
scale_grace_maps.py
===================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports scaled spatial fields
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms :cite:p:`Jekeli:1981vj` :cite:p:`Swenson:2006hu`
- Converts data to centimeters water equivalent, performs a spherical harmonic summation to convert to the spatial domain :cite:p:`Wahr:1998hy`
- Scales the spatial fields following :cite:p:`Landerer:2012kf`
- Calculates the scaled spatial error field following :cite:p:`Wahr:2006bx`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/scale_grace_maps.py

Calling Sequence
################

.. argparse::
    :filename: scale_grace_maps.py
    :func: arguments
    :prog: scale_grace_maps.py
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

    --gia -G : @after
        * ``'IJ05-R2'``: Ivins R2 GIA Models :cite:p:`Ivins:2013cq`
        * ``'W12a'``: Whitehouse GIA Models :cite:p:`Whitehouse:2012jj`
        * ``'SM09'``: Simpson/Milne GIA Models :cite:p:`Simpson:2009hg`
        * ``'ICE6G'``: ICE-6G GIA Models :cite:p:`Peltier:2015bo`
        * ``'Wu10'``: Wu (2010) GIA Correction :cite:p:`Wu:2010dq`
        * ``'AW13-ICE6G'``: Geruo A ICE-6G GIA Models :cite:p:`A:2013kh`
        * ``'AW13-IJ05'``: Geruo A IJ05-R2 GIA Models :cite:p:`A:2013kh`
        * ``'Caron'``: Caron JPL GIA Assimilation :cite:p:`Caron:2018ba`
        * ``'ICE6G-D'``: ICE-6G Version-D GIA Models :cite:p:`Peltier:2018dp`
        * ``'ascii'``: reformatted GIA in ascii format
        * ``'netCDF4'``: reformatted GIA in netCDF4 format
        * ``'HDF5'``: reformatted GIA in HDF5 format

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

    --slr-c40 : @replace
        Replace *C*\ :sub:`40` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (5\ |times|\ 5 with 6,1)
        * ``'GSFC'``: use values from GSFC
        * ``'LARES'``: use filtered values from CSR

    --slr-c50 : @replace
        Replace *C*\ :sub:`50` coefficients with SLR values

        * ``None``: use original values
        * ``'CSR'``: use values from CSR (5\ |times|\ 5 with 6,1)
        * ``'GSFC'``: use values from GSFC
        * ``'LARES'``: use filtered values from CSR

    --interval : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
