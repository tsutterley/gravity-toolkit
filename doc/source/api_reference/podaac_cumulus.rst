=================
podaac_cumulus.py
=================

- Syncs GRACE/GRACE-FO data from `NASA JPL PO.DAAC Cumulus AWS S3 bucket <https://podaac.jpl.nasa.gov/cloud-datasets/about>`_
- S3 Cumulus syncs are only available in AWS instances in ``us-west-2``
- Creates an index file for each data product

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/podaac_cumulus.py

Calling Sequence
################

.. argparse::
    :filename: podaac_cumulus.py
    :func: arguments
    :prog: podaac_cumulus.py
    :nodescription:
    :nodefault:
