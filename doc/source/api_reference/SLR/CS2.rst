=======
SLR.CS2
=======

- Reads monthly degree 2,m (figure axis and azimuthal dependence) spherical harmonic data files from satellite laser ranging (SLR)

    * CSR 2,1: ``C21_S21_RL06.txt``
    * CSR 2,2: ``C22_S22_RL06.txt``
    * GFZ: ``GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat``
    * GSFC: ``gsfc_slr_5x5c61s61.txt``

Calling Sequence
################

.. code-block:: python

    import gravity_toolkit.SLR
    SLR_CS2 = gravity_toolkit.SLR.CS2(SLR_file)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/SLR/CS2.py

.. autofunction:: gravity_toolkit.SLR.CS2
