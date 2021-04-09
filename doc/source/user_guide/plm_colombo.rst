==============
plm_colombo.py
==============

- Computes fully-normalized associated Legendre Polynomials and their first derivative for a vector of x values using a standard forward column method
- Uses the [Colombo1981]_ recursion relation listed in [Losch2003]_ and Holmes2002]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.plm_colombo import plm_colombo
    plm,dplm = plm_colombo(LMAX, x)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/plm_colombo.py

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
- ``dplms``: first differentials of Legendre polynomials of ``x``

References
##########

.. [Colombo1981] O. L. Colombo, "Numerical Methods for Harmonic Analysis on the Sphere", Air Force Contract No. F19628-79-C-0027, OSURF Proj. No. 711664, 140 pp., (1981).

.. [Losch2003] M. Losch and V. Seufer, "How to Compute Geoid Undulations (Geoid Height Relative to a Given Reference Ellipsoid) from Spherical Harmonic Coefficients for Satellite Altimetry Applications", (2003). `eprint ID: 11802 <http://mitgcm.org/~mlosch/geoidcookbook.pdf>`_

.. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002). `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_
