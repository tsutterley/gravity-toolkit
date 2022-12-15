============
read_SLR_C40
============

- Reads monthly degree 4 zonal spherical harmonic data files from satellite laser ranging (SLR)

    * CSR: ``CSR_Monthly_5x5_Gravity_Harmonics.txt``
    * GSFC: ``gsfc_slr_5x5c61s61.txt``
    * LARES: ``C40_LARES_filtered.txt``

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_C40 import read_SLR_C40
    SLR_C40 = read_SLR_C40(SLR_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C40.py

.. autofunction:: gravity_toolkit.read_SLR_C40
