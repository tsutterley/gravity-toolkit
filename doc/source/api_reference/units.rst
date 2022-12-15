=====
units
=====

Class for converting GRACE/GRACE-FO Level-2 data to specific units

Calling Sequence
================

Calculating the degree dependent factors for converting harmonic units

.. code-block:: python

    from gravity_toolkit.units import units
    from gravity_toolkit.read_love_numbers import read_love_numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, REFERENCE='CF')
    dfactor = units(lmax=lmax).harmonic(hl,kl,ll)
    to_cmwe = dfactor.cmwe

Calculating the degree dependent factors for converting spatial units

.. code-block:: python

    from gravity_toolkit.units import units
    from gravity_toolkit.read_love_numbers import read_love_numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, REFERENCE='CF')
    dfactor = units(lmax=lmax).spatial(hl,kl,ll)
    from_cmwe = dfactor.cmwe

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/units.py

General Attributes and Methods
==============================

.. autoclass:: gravity_toolkit.units
   :members:

References
----------
.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan,
    "Time variability of the Earth's gravity field:
    Hydrological and oceanic effects and their possible
    detection using GRACE", *Journal of Geophysical Research*,
    103(B12), 30205--30229, (1998).
    `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_
.. [Wahr2000] J. Wahr, D. Wingham, and C. Bentley,
    "A method of combining ICESat and GRACE satellite data
    to constrain Antarctic mass balance",
    *Journal of Geophysical Research: Solid Earth*,
    105(B7), 16279--16294, (2000).
    `doi: 10.1029/2000JB900113 <https://doi.org/10.1029/2000JB900113>`_
