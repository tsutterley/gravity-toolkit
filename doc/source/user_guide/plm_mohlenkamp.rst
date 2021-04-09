=================
plm_mohlenkamp.py
=================

- Computes fully-normalized associated Legendre Polynomials for a vector of x values using Martin Mohlenkamp's recursion relation [Mohlenkamp2016]_
- Derived from [Szego1939]_ recurrence formula for Jacobi Polynomials


Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp
    plm = plm_mohlenkamp(LMAX, x)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/plm_mohlenkamp.py

Arguments
#########

- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``x``: elements ranging from -1 to 1. Typically ``cos(theta)``, where ``theta`` is the colatitude in radians

Keyword arguments
#################

- ``MMAX``: Upper bound of Spherical Harmonic Orders (default: ``LMAX``)

Returns
#######

- ``plms``: Legendre polynomials of ``x`` (geodesy normalization)

References
##########

.. [Mohlenkamp2016] M. J. Mohlenkamp, "A User’s Guide to Spherical Harmonics", (2016). `[pdf] <http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf>`_

.. [Szego1939] Gabor Szeg\ |ouml|\ , "Orthogonal Polynomials", 440 pp., (1939) `[pdf] <https://people.math.osu.edu/nevai.1/AT/SZEGO/szego=szego1975=ops=OCR.pdf>`_

.. |ouml|    unicode:: U+00F6 .. LATIN SMALL LETTER O WITH DIAERESIS
