====================
combine_harmonics.py
====================

- Converts a file from the spherical harmonic domain into the spatial domain :cite:p:`Wahr:1998hy`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/combine_harmonics.py

Calling Sequence
################

.. argparse::
    :filename: combine_harmonics.py
    :func: arguments
    :prog: combine_harmonics.py
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

    --units -U : @after
        * ``0``: no unit conversion
        * ``1``: cm of water thickness
        * ``2``: mm of geoid height
        * ``3``: mm of elastic crustal deformation
        * ``4``: |mu|\ Gal gravitational perturbation
        * ``5``: mbar equivalent surface pressure

    --interval -I : @after
        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)

.. |mu|      unicode:: U+03BC .. GREEK SMALL LETTER MU
