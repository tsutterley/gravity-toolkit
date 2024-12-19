=========================
run_sea_level_equation.py
=========================

- Solves the sea level equation with the option of including polar motion feedback :cite:p:`Farrell:1976hm` :cite:p:`Kendall:2005ds` :cite:p:`Mitrovica:2003cq`
- Uses a Clenshaw summation to calculate the spherical harmonic summation :cite:p:`Holmes:2002ff` :cite:p:`Tscherning:1982tu`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/run_sea_level_equation.py


Calling Sequence
################

.. argparse::
    :filename: run_sea_level_equation.py
    :func: arguments
    :prog: run_sea_level_equation.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM :cite:p:`Han:1995go`
        * ``1``: Gegout (2005) values from PREM :cite:p:`Gegout:2010gc`
        * ``2``: Wang et al. (2012) values from PREM :cite:p:`Wang:2012gc`
        * ``3``: Wang et al. (2012) values from PREM with hard sediment :cite:p:`Wang:2012gc`
        * ``4``: Wang et al. (2012) values from PREM with soft sediment :cite:p:`Wang:2012gc`

    --body -b : @after
        * ``0``: :cite:p:`Wahr:1981ea` and :cite:p:`Wahr:1985gr` values from PREM
        * ``1``: :cite:p:`Farrell:1972cm` values from Gutenberg-Bullen oceanic mantle model

    --fluid -f : @after
        * ``0``: :cite:p:`Han:1989kj` fluid love number
        * ``1``: :cite:p:`Munk:1960uk` secular love number
        * ``2``: :cite:p:`Munk:1960uk` fluid love number
        * ``3``: :cite:p:`Lambeck:1980um` fluid love number

    --polar-feedback : @replace
        Include polar feedback :cite:p:`Wahr:1985gr`

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth
