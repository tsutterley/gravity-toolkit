==================
read_gfc_harmonics
==================

- Reads gfc files and extracts spherical harmonics for Swarm and GRAZ GRACE/GRACE-FO data
- Parses date of GRACE/GRACE-FO data from filename

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
    Ylms = read_gfc_harmonics(input_file, TIDE='tide_free')

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/read_gfc_harmonics.py

.. autofunction:: gravity_toolkit.read_gfc_harmonics
