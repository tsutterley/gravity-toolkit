================
gen_harmonics.py
================

- Converts data from the spatial domain to spherical harmonic coefficients
- Does not compute the solid Earth elastic response or convert units

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_harmonics import gen_harmonics
    from gravity_toolkit.plm_holmes import plm_holmes
    PLM,dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_harmonics(data, lon, lat, LMAX=LMAX, PLM=PLM)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_harmonics.py

.. autofunction:: gravity_toolkit.gen_harmonics

.. autofunction:: gravity_toolkit.gen_harmonics.integration

.. autofunction:: gravity_toolkit.gen_harmonics.fourier
