=======================
legendre_polynomials.py
=======================

- Computes fully normalized Legendre polynomials for an array of x values and their first derivative following [HofmannWellenhof2006]_
- Calculates Legendre polynomials for zonal harmonics (order 0)

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.legendre_polynomials import legendre_polynomials
    pl,dpl = legendre_polynomials(LMAX, x)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/legendre_polynomials.py

Arguments
#########

- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``x``: elements ranging from -1 to 1. Typically ``cos(theta)``, where ``theta`` is the colatitude in radians

Keyword arguments
#################

- ``ASTYPE``: output variable type. Default is 64-bit floating point

Returns
#######

- ``pl``: Legendre polynomials of ``x`` (geodesy normalization)
- ``dpl``: first differentials of Legendre polynomials of ``x``

References
##########

.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
