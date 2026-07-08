==================
sea_level_equation
==================

- Solves the sea level equation with the option of including polar motion feedback

Calling Sequence
################

.. code-block:: python

     from gravity_toolkit.associated_legendre import plm_holmes
     from gravity_toolkit.sea_level_equation import sea_level_equation
     PLM, dPLM = plm_holmes(LMAX, np.cos(th))
     Ylms = sea_level_equation(loadClm, loadSlm, lon, lat, land_function,
         LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll), ITERATIONS=6, POLAR=True)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/sea_level_equation.py

.. autofunction:: gravity_toolkit.sea_level_equation
