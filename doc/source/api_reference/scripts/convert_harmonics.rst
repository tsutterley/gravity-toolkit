====================
convert_harmonics.py
====================

- Converts a file from the spatial domain into the spherical harmonic domain :cite:p:`Wahr:1998hy`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/convert_harmonics.py

Calling Sequence
################

.. argparse::
    :filename: convert_harmonics.py
    :func: arguments
    :prog: convert_harmonics.py
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
        * ``1``: cm of water thickness (cm.w.e., g/cm\ :sup:`2`)
        * ``2``: Gigatonnes (Gt)
        * ``3``: mm of water thickness (kg/m\ :sup:`2`)

    --interval -I : @after
        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
