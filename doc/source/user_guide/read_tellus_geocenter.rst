========================
read_tellus_geocenter.py
========================

- Reads monthly geocenter spherical harmonic data files from `GRACE Tellus Technical Notes <https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1>`_
- Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_ [Sun2016]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_tellus_geocenter import read_tellus_geocenter
    deg1_input = read_tellus_geocenter(geocenter_file, JPL=True)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_tellus_geocenter.py

Arguments
#########

- ``geocenter_file``: degree 1 file

    * CSR: ``TN-13_GEOC_CSR_RL06.txt``
    * GFZ: ``TN-13_GEOC_GFZ_RL06.txt``
    * JPL: ``TN-13_GEOC_JPL_RL06.txt``

Keyword arguments
#################

- ``HEADER``: file contains header text to be skipped (default: True)
- ``JPL``: use JPL TN-13 geocenter files calculated following [Sun2016]_

Returns
#######
- ``C10``: Cosine degree 1, order 0 spherical harmonic coefficients
- ``C11``: Cosine degree 1, order 1 spherical harmonic coefficients
- ``S11``: Sine degree 1, order 1 spherical harmonic coefficients
- ``eC10``: Cosine degree 1, order 0 spherical harmonic coefficients Error
- ``eC11``: Cosine degree 1, order 1 spherical harmonic coefficients Error
- ``eS11``: Sine degree 1, order 1 spherical harmonic coefficients Error
- ``month``: GRACE/GRACE-FO month (April 2002 = 004)
- ``time``: date of GRACE/GRACE-FO month in decimal format

References
##########

.. [Sun2016] Y. Sun, P. Ditmar, and R. Riva, "Observed changes in the Earth's dynamic oblateness from GRACE data and geophysical models", *Journal of Geodesy*, 90(1), 81--89, (2016). `doi: 10.1007/s00190-015-0852-y <https://doi.org/10.1007/s00190-015-0852-y>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
