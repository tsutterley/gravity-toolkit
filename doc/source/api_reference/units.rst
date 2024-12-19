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

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/units.py

General Attributes and Methods
==============================

.. autoclass:: gravity_toolkit.units
   :members:

