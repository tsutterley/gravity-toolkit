===================
associated_legendre
===================

- Computes fully-normalized associated Legendre Polynomials and their first derivative for a vector of ``x`` values

Calling Sequence
################

.. code-block:: python

    import gravity_toolkit.associated_legendre
    PLM, dPLM = gravtk.associated_legendre.polynomials(LMAX, x, method='holmes')

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/associated_legendre.py

.. autofunction:: gravity_toolkit.associated_legendre

.. autofunction:: gravity_toolkit.plm_colombo

.. autofunction:: gravity_toolkit.plm_holmes

.. autofunction:: gravity_toolkit.plm_mohlenkamp
