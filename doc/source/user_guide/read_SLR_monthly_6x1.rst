=======================
read_SLR_monthly_6x1.py
=======================

- Reads in monthly 5x5 spherical harmonic coefficients with 1 coefficient from degree 6 all calculated from satellite laser ranging (SLR) measurements
- Calculated by the University of Texas Center for Space Research (CSR) [Cheng2010]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_monthly_6x1 import read_SLR_monthly_6x1
    Ylms = read_SLR_monthly_6x1(input_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_monthly_6x1.py

Arguments
#########

- ``input_file``: input satellite laser ranging file

Keyword arguments
#################

- ``HEADER``: file contains header text to be skipped

Returns
#######

- ``clm``: Cosine spherical harmonic coefficients
- ``slm``: Sine spherical harmonic coefficients
- ``error/clm``: Cosine spherical harmonic coefficient uncertainty
- ``error/slm``: Sine spherical harmonic coefficients uncertainty
- ``MJD``: output date as Modified Julian Day
- ``date``: output date in year-decimal

References
##########

.. [Cheng2010] M. Cheng, J. C. Ries, and B. D. Tapley, "Variations of the Earth's figure axis from satellite laser ranging and GRACE", *Journal of Geophysical Research*, 116(B01409), `doi: 10.1029/2010JB000850 <https://doi.org/10.1029/2010JB000850>`_
