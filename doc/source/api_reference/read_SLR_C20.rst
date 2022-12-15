============
read_SLR_C20
============

- Reads monthly oblateness (degree 2 zonal) spherical harmonic data files from satellite laser ranging (SLR)

    * RL04: ``TN-05_C20_SLR.txt``
    * RL05: ``TN-07_C20_SLR.txt``
    * RL06: ``TN-11_C20_SLR.txt``
    * CSR: ``C20_RL05.txt``
    * GSFC: ``TN-14_C30_C30_GSFC_SLR.txt``

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_C20 import read_SLR_C20
    SLR_C20 = read_SLR_C20(SLR_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C20.py

.. autofunction:: gravity_toolkit.read_SLR_C20
