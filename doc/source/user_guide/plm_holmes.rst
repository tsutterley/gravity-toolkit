=============
plm_holmes.py
=============

- Computes fully-normalized associated Legendre Polynomials and their first derivative for a vector of x values using the [Holmes2002]_ recursion relation
- Recursion relation is stable up to very high degree and order

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.plm_holmes import plm_holmes
    plm,dplm = plm_holmes(LMAX, x)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/plm_holmes.py

Arguments
#########

- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``x``: elements ranging from -1 to 1. Typically ``cos(theta)``, where ``theta`` is the colatitude in radians

Keyword arguments
#################

- ``ASTYPE``: output variable type. Default is 64-bit floating point

Returns
#######

- ``plms``: Legendre polynomials of ``x`` (geodesy normalization)
- ``dplms``: first differentials of Legendre polynomials of `x`

References
##########

.. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002). `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_
