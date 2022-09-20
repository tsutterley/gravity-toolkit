====================
gen_spherical_cap.py
====================

- Calculates gravitational spherical harmonic coefficients for a spherical cap

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_spherical_cap import gen_spherical_cap
    from gravity_toolkit.plm_holmes import plm_holmes
    PLM, dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_spherical_cap(data, lon, lat, UNITS=1, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_spherical_cap.py

.. autofunction:: gravity_toolkit.gen_spherical_cap
