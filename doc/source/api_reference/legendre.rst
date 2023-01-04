========
legendre
========

- Computes associated Legendre functions of degree ``l`` evaluated for elements ``x``
- ``l`` must be a scalar integer and ``x`` must contain real values ranging -1 <= ``x`` <= 1
- Unnormalized associated Legendre function values will overflow for ``l`` > 150

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.legendre import legendre
    Pl = legendre(l, x)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/legendre.py

.. autofunction:: gravity_toolkit.legendre
