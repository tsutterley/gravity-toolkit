===============
read_SLR_C30.py
===============

- Reads monthly degree 3 zonal spherical harmonic data files from satellite laser ranging (SLR)

    * CSR: ``CSR_Monthly_5x5_Gravity_Harmonics.txt``
    * GFZ: ``GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat``
    * GSFC: ``TN-14_C30_C30_GSFC_SLR.txt``
    * LARES: ``C30_LARES_filtered.txt``

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_C30 import read_SLR_C30
    SLR_C30 = read_SLR_C30(SLR_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C30.py


.. autofunction:: gravity_toolkit.read_SLR_C30
