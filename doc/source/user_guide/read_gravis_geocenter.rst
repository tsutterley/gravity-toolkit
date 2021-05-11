========================
read_gravis_geocenter.py
========================

- Reads `monthly geocenter coefficients <ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/GRAVIS-2B_GFZOP_GEOCENTER_0002.dat>`_ from GFZ GravIS
- Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_ [Sun2016]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_gravis_geocenter import read_gravis_geocenter
    deg1_input = read_gravis_geocenter(geocenter_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_gravis_geocenter.py

Arguments
#########

- ``geocenter_file``: degree 1 file

Keyword arguments
#################

- ``HEADER``: file contains header text to be skipped

Returns
#######

- ``C10``: Cosine degree 1, order 0 spherical harmonic coefficients
- ``C11``: Cosine degree 1, order 1 spherical harmonic coefficients
- ``S11``: Sine degree 1, order 1 spherical harmonic coefficients
- ``month``: GRACE/GRACE-FO month (April 2002 = 004)
- ``time``: date of GRACE/GRACE-FO month in decimal format

References
##########

.. [Sun2016] Y. Sun, P. Ditmar, and R. Riva, "Observed changes in the Earth's dynamic oblateness from GRACE data and geophysical models", *Journal of Geodesy*, 90(1), 81--89, (2016). `doi: 10.1007/s00190-015-0852-y <https://doi.org/10.1007/s00190-015-0852-y>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
