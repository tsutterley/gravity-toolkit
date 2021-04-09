=========================
read_swenson_geocenter.py
=========================

- Reads `monthly geocenter coefficients <https://github.com/swensosc/GRACE_Tiles/blob/master/ancillary_data/gad_gsm.rl05.txt>`_ in units mm w.e
- Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_swenson_geocenter import read_swenson_geocenter
    deg1_input = read_swenson_geocenter(geocenter_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_swenson_geocenter.py

Arguments
#########

- ``geocenter_file``: degree 1 file

Keyword arguments
#################

- ``HEADER``: file contains header text to be skipped (default: True)

Returns
#######

- ``C10``: Cosine degree 1, order 0 spherical harmonic coefficients
- ``C11``: Cosine degree 1, order 1 spherical harmonic coefficients
- ``S11``: Sine degree 1, order 1 spherical harmonic coefficients
- ``month``: GRACE/GRACE-FO month (April 2002 = 004)
- ``time``: date of GRACE/GRACE-FO month in decimal format

References
##########

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
