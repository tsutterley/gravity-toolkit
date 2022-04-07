=================
gen_point_load.py
=================

- Calculates gravitational spherical harmonic coefficients for point masses

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_point_load import gen_point_load
    Ylms = gen_point_load(data, lon, lat, UNITS=1, LMAX=LMAX, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_point_load.py

.. autofunction:: gravity_toolkit.gen_point_load
