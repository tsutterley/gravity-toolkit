=============
gen_disc_load
=============

- Calculates gravitational spherical harmonic coefficients for a uniform disc load

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_disc_load import gen_disc_load
    from gravity_toolkit.associated_legendre import plm_holmes
    PLM, dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_disc_load(data, lon, lat, area, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/gen_disc_load.py

.. autofunction:: gravity_toolkit.gen_disc_load
