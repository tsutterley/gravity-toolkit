==================
aod1b_geocenter.py
==================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product

    * ``'atm'``: atmospheric loading from ECMWF
    * ``'ocn'``: oceanic loading from OMCT/MPIOM
    * ``'glo'``: global atmospheric and oceanic loading
    * ``'oba'``: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of geocenter variations at 3 or 6-hour intervals

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/aod1b_geocenter.py

Calling Sequence
################

.. argparse::
    :filename: aod1b_geocenter.py
    :func: arguments
    :prog: aod1b_geocenter.py
    :nodescription:
    :nodefault:
