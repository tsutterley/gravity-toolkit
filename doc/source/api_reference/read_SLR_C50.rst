===============
read_SLR_C50.py
===============

- Reads monthly degree 5 zonal spherical harmonic data files from satellite laser ranging (SLR)

    * CSR: ``CSR_Monthly_5x5_Gravity_Harmonics.txt``
    * GSFC: ``gsfc_slr_5x5c61s61.txt``
    * LARES: ``C50_LARES_filtered.txt``

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_C50 import read_SLR_C50
    SLR_C50 = read_SLR_C50(SLR_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C50.py

.. autofunction:: gravity_toolkit.read_SLR_C50
