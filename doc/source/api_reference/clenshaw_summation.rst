==================
clenshaw_summation
==================

- Returns the spatial field for a series of spherical harmonics at a sequence of ungridded points
- Uses a Clenshaw summation to calculate the spherical harmonic summation

Calling Sequence
################

.. code-block:: python

   from gravity_toolkit.clenshaw_summation import clenshaw_summation
   spatial = clenshaw_summation(clm,slm,lon,lat,UNITS=1,LMAX=60,LOVE=LOVE)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/clenshaw_summation.py

.. autofunction:: gravity_toolkit.clenshaw_summation
