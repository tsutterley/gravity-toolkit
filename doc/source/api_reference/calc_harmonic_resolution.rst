===========================
calc_harmonic_resolution.py
===========================

- Calculates the spatial resolution that can be resolved by the spherical harmonics of a certain degree [Barthelmes2013]_ [HofmannWellenhof2006]_
- Default method uses the smallest half-wavelength that can be resolved
- Secondary method calculates the smallest possible bump that can be resolved

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/calc_harmonic_resolution.py

Calling Sequence
################

.. argparse::
    :filename: ../scripts/calc_harmonic_resolution.py
    :func: arguments
    :prog: calc_harmonic_resolution.py
    :nodescription:
    :nodefault:

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
