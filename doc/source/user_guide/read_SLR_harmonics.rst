=====================
read_SLR_harmonics.py
=====================

- Reads 5\ |times|\ 5 spherical harmonic coefficients with 1 coefficient from degree 6 all calculated from satellite laser ranging (SLR) measurements
- Calculated by the University of Texas Center for Space Research (CSR) and NASA Goddard Space Flight Center (GSFC)

    * CSR: ``CSR_Monthly_5x5_Gravity_Harmonics.txt``
    * GSFC: ``gsfc_slr_5x5c61s61.txt``

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_harmonics import read_SLR_harmonics
    Ylms = read_SLR_harmonics(input_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_harmonics.py

.. autofunction:: gravity_toolkit.read_SLR_harmonics

.. autofunction:: gravity_toolkit.read_SLR_harmonics.read_CSR_monthly_6x1

.. autofunction:: gravity_toolkit.read_SLR_harmonics.read_GSFC_weekly_6x1

.. autofunction:: gravity_toolkit.read_SLR_harmonics.convert_weekly

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
