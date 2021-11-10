=====================
read_SLR_harmonics.py
=====================

- Reads 5x5 spherical harmonic coefficients with 1 coefficient from degree 6 all calculated from satellite laser ranging (SLR) measurements
- Calculated by the University of Texas Center for Space Research (CSR) [Cheng2010]_ and NASA Goddard Space Flight Center (GSFC) [Loomis2020]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_SLR_harmonics import read_SLR_harmonics
    Ylms = read_SLR_harmonics(input_file)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_harmonics.py

Arguments
#########

- ``input_file``: input satellite laser ranging file
    * CSR: `CSR_Monthly_5x5_Gravity_Harmonics.txt`
    * GSFC: `gsfc_slr_5x5c61s61.txt`

Keyword arguments
#################

- ``SCALE``: scale factor for converting to fully-normalized spherical harmonics
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

.. [Loomis2020] B. D. Loomis, K. E. Rachlin, D. N. Wiese, F. W. Landerer, and S. B. Luthcke, "Replacing GRACE/GRACE‐FO *C*\ :sub:`30` with satellite laser ranging: Impacts on Antarctic Ice Sheet mass change". *Geophysical Research Letters*, 47, (2020). `doi: 10.1029/2019GL085488 <https://doi.org/10.1029/2019GL085488>`_
