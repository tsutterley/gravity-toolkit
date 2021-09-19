===================
fourier_legendre.py
===================

- Computes Fourier coefficients of the associated Legendre functions

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.fourier_legendre import fourier_legendre
    plm = fourier_legendre(lmax, mmax)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/fourier_legendre.py

Arguments
#########

- ``lmax``: Upper bound of Spherical Harmonic Degrees
- ``mmax``: Upper bound of Spherical Harmonic Orders

Returns
#######

- ``plm``: Fourier coefficients
