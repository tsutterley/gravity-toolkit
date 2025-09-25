===========================
dealiasing_global_uplift.py
===========================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for global atmospheric and oceanic loading and estimates anomalies in elastic crustal uplift :cite:p:`Davis:2004il,Wahr:1998hy`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/dealiasing/dealiasing_global_uplift.py

Calling Sequence
################

.. argparse::
    :filename: dealiasing_global_uplift.py
    :func: arguments
    :prog: dealiasing_global_uplift.py
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

    --interval -I : @after
        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)
