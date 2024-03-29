=======================
piecewise_grace_maps.py
=======================

- Reads in GRACE/GRACE-FO spatial files and fits a piecewise regression model at each grid point for breakpoint analysis

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/piecewise_grace_maps.py

Calling Sequence
################

.. argparse::
    :filename: piecewise_grace_maps.py
    :func: arguments
    :prog: piecewise_grace_maps.py
    :nodescription:
    :nodefault:

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
