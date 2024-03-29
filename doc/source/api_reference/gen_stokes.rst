==========
gen_stokes
==========

- Converts data from the spatial domain to spherical harmonic coefficients

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_stokes import gen_stokes
    from gravity_toolkit.associated_legendre import plm_holmes
    PLM, dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_stokes(data, lon, lat, UNITS=1, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/gen_stokes.py

.. autofunction:: gravity_toolkit.gen_stokes
