=======================
read_GRACE_harmonics.py
=======================

- Reads GRACE/GRACE-FO files and extracts spherical harmonic data and drift rates (RL04)
- Adds drift rates to clm and slm for Release-4 harmonics
- Correct Release-5 GSM data for drift in pole tide
- Extracts start and end date of GRACE/GRACE-FO files and calculates mean of range

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
    CSR_L2_input = read_GRACE_harmonics('GSM-2_2002095-2002120_0021_UTCSR_0060_0005.gz',60)
    GFZ_L2_input = read_GRACE_harmonics('GSM-2_2002094-2002120_0024_EIGEN_G---_005a.gz',90)
    JPL_L2_input = read_GRACE_harmonics('GSM-2_2002091-2002120_0018_JPLEM_0001_0005.gz',60)
    JPLMSC_input = read_GRACE_harmonics('GSM-2_2003001-2003031_0029_JPLMSC_0719_0005',719)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_GRACE_harmonics.py

.. autofunction:: gravity_toolkit.read_GRACE_harmonics

.. autofunction:: gravity_toolkit.read_GRACE_harmonics.parse_file

.. autofunction:: gravity_toolkit.read_GRACE_harmonics.extract_file
