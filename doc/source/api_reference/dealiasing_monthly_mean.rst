==========================
dealiasing_monthly_mean.py
==========================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product and outputs monthly the mean for a specific GRACE/GRACE-FO processing center and data release

    * ``'GAA'``: atmospheric loading from ECMWF
    * ``'GAB'``: oceanic loading from OMCT/MPIOM
    * ``'GAC'``: global atmospheric and oceanic loading
    * ``'GAD'``: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of oblateness variations at 3 or 6-hour intervals

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/dealiasing_monthly_mean.py

Calling Sequence
################

.. argparse::
    :filename: dealiasing_monthly_mean.py
    :func: arguments
    :prog: dealiasing_monthly_mean.py
    :nodescription:
    :nodefault:
