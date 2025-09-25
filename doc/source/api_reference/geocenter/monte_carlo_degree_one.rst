=========================
monte_carlo_degree_one.py
=========================

- Estimates uncertainties in degree 1 using GRACE/GRACE-FO coefficients of degree 2 and greater, and modeled ocean bottom pressure variations in a Monte Carlo scheme :cite:p:`Swenson:2008cr,Sutterley:2019bx`.
- Calculates the estimated spherical harmonic errors following :cite:t:`Wahr:2006bx`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/geocenter/monte_carlo_degree_one.py

Calling Sequence
################

.. argparse::
    :filename: monte_carlo_degree_one.py
    :func: arguments
    :prog: monte_carlo_degree_one.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM :cite:p:`Han:1995go`
        * ``1``: Gegout (2005) values from PREM :cite:p:`Gegout:2010gc`
        * ``2``: Wang et al. (2012) values from PREM :cite:p:`Wang:2012gc`
        * ``3``: Wang et al. (2012) values from PREM with hard sediment :cite:p:`Wang:2012gc`
        * ``4``: Wang et al. (2012) values from PREM with soft sediment :cite:p:`Wang:2012gc`

    --kl -k : @after
        * ``None``: use derived values from :cite:p:`Trupin:1992kp,Blewitt:2003bz`.

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

    --solver -s : @replace
        Least squares solver for degree one solutions

        * ``'inv'``: matrix inversion
        * ``'lstsq'``: least squares solution
        * ``'gelsy'``: complete orthogonal factorization solution
        * ``'gelss'``: singular value decomposition (SVD) solution
        * ``'gelsd'``: singular value decomposition (SVD) solution with a divide and conquer method

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
