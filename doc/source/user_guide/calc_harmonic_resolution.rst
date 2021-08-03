===========================
calc_harmonic_resolution.py
===========================

- Calculates the spatial resolution that can be resolved by the spherical harmonics of a certain degree [Barthelmes2013]_ [HofmannWellenhof2006]_
- Default method uses the smallest half-wavelength that can be resolved
- Secondary method calculates the smallest possible bump that can be resolved

Calling Sequence
################

.. code-block:: bash

     python calc_harmonic_resolution --lmax 60 --cap

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/calc_harmonic_resolution.py

Command Line Options
####################

- ``--help``: list the command line options
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-R X``, ``--radius X``: Average radius of the Earth (km)
- ``-C``, ``--cap``: Calculate the smallest possible bump that can be resolved

References
##########

.. [Barthelmes2013] F. Barthelmes, "Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models", GeoForschungsZentrum Scientific Technical Report, STR09/02, (2013). `doi: 10.2312/GFZ.b103-0902-26 <https://doi.org/10.2312/GFZ.b103-0902-26>`_

.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
