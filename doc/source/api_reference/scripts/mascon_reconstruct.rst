=====================
mascon_reconstruct.py
=====================

- Calculates the equivalent spherical harmonics from a mascon time series

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/mascon_reconstruct.py

Calling Sequence
################

.. argparse::
    :filename: mascon_reconstruct.py
    :func: arguments
    :prog: mascon_reconstruct.py
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

    --fit-method : @after
        * ``1``: mass coefficients
        * ``2``: geoid coefficients

