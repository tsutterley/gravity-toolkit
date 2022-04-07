=============
plm_holmes.py
=============

- Computes fully-normalized associated Legendre Polynomials and their first derivative for a vector of ``x`` values using the Holmes and Featherstone recursion relation
- Recursion relation is stable up to very high degree and order

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.plm_holmes import plm_holmes
    plm,dplm = plm_holmes(LMAX, x)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/plm_holmes.py

.. autofunction:: gravity_toolkit.plm_holmes
