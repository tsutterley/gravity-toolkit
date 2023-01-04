====================
podaac_grace_sync.py
====================

- Syncs GRACE/GRACE-FO and auxiliary data from the `NASA JPL PO.DAAC Drive Server <https://podaac-tools.jpl.nasa.gov/drive>`_
- Syncs CSR/GFZ/JPL files for RL04/RL05/RL06 GAA/GAB/GAC/GAD/GSM (GAA and GAB are GFZ/JPL only)
- Gets the latest technical note (TN) files
- Gets the monthly GRACE/GRACE-FO newsletters
- Creates an index file for each data product

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/podaac_grace_sync.py

Calling Sequence
################

.. argparse::
    :filename: podaac_grace_sync.py
    :func: arguments
    :prog: podaac_grace_sync.py
    :nodescription:
    :nodefault:
