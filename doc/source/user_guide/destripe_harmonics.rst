=====================
destripe_harmonics.py
=====================

- Filters spherical harmonic coefficients for correlated "striping" errors following [Swenson2006]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.destripe_harmonics import destripe_harmonics
    Ylms = destripe_harmonics(clm,slm,LMAX=60)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/destripe_harmonics.py

Arguments
#########

1. ``clm``: cosine spherical harmonic coefficients
2. ``slm``: sine spherical harmonic coefficients

Keyword arguments
#################

- ``LMIN``: Lower bound of Spherical Harmonic Degrees
- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``MMAX``: Upper bound of Spherical Harmonic Orders
- ``ROUND``: use round to find nearest even (``True``) or use floor (``False``)
- ``NARROW``: sets harmonics to zero if number of points is less than window size (``False``)

Returns
#######

- ``Wclm``: filtered cosine spherical harmonic coefficients
- ``Wslm``: filtered sine spherical harmonic coefficients

References
##########

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_
