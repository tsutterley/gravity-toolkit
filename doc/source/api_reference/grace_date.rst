==========
grace_date
==========

- Reads GRACE/GRACE-FO index file from `podaac_cumulus.py` or `gfz_isdc_grace_ftp.py`
- Parses dates of each GRACE/GRACE-FO file and assigns the month number
- Creates an index of dates for GRACE/GRACE-FO files

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.grace_date import grace_date
    grace_files = grace_date(base_dir, PROC=PROC, DREL=DREL, DSET=DSET)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/grace_date.py

.. autofunction:: gravity_toolkit.grace_date

