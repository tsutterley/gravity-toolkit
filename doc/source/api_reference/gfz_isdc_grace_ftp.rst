=====================
gfz_isdc_grace_ftp.py
=====================

- Syncs GRACE/GRACE-FO and auxiliary data from the `GFZ Information System and Data Center (ISDC) <http://isdc.gfz-potsdam.de/grace-isdc/>`_
- Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
- Gets the latest technical note (TN) files
- Gets the monthly GRACE/GRACE-FO newsletters
- Creates an index file for each data product

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/gfz_isdc_grace_ftp.py

Calling Sequence
################

.. argparse::
    :filename: gfz_isdc_grace_ftp.py
    :func: arguments
    :prog: gfz_isdc_grace_ftp.py
    :nodescription:
    :nodefault:
