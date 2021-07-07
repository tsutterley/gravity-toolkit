=====================
grace_input_months.py
=====================

- Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range
- Includes Degree 1 coefficients with input values (if specified)
- Replaces C20 with SLR values (if specified)
- Replaces C21/S21/C22/S22/C30/C50 with SLR values for months 179+ (if specified)
- Corrects for ECMWF atmospheric "jumps" using the GAE, GAF and GAG files following [Fagiolini2015]_
- Corrects for Pole Tide drift following [Wahr2015]_

Calling Sequence
################

.. code-block:: python

      from gravity_toolkit.grace_input_months import grace_input_months
      GRACE_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        start_mon, end_mon, missing, SLR_C20, DEG1, SLR_C30=SLR_C30)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/grace_input_months.py

Arguments
#########

1. ``base_dir``: Working data directory for GRACE/GRACE-FO data
2. ``PROC``: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)

      * ``'CSR'``: University of Texas Center for Space Research
      * ``'GFZ'``: German Research Centre for Geosciences (GeoForschungsZentrum)
      * ``'JPL'``: Jet Propulsion Laboratory
      * ``'CNES'``: French Centre National D'Etudes Spatiales
3. ``DREL``: GRACE/GRACE-FO data release (RL04, RL05, RL06)
4. ``DSET``: GRACE/GRACE-FO data product (GAA, GAB, GAC, GAD, GSM)

      * ``'GAA'``: non-tidal atmospheric correction
      * ``'GAB'``: non-tidal oceanic correction
      * ``'GAC'``: combined non-tidal atmospheric and oceanic correction
      * ``'GAD'``: GRACE/GRACE-FO ocean bottom pressure product
      * ``'GSM'``: corrected monthly GRACE/GRACE-FO static field product
5. ``LMAX``: Upper bound of Spherical Harmonic Degrees
6. ``start_mon``: starting month to consider in analysis
7. ``end_mon``: ending month to consider in analysis
8. ``missing``: missing months to not consider in analysis
9. ``SLR_C20``: Replaces C20 with values from Satellite Laser Ranging (SLR)

      * ``None``: use original values
      * ``'CSR'``: use values from CSR (TN-07, TN-09, TN-11)
      * ``'GFZ'``: use values from GFZ
      * ``'GSFC'``: use values from GSFC (TN-14)
10. ``DEG1``: Use Degree 1 coefficients

      * ``None``: No degree 1 replacement
      * ``'Tellus'``: `GRACE/GRACE-FO TN-13 coefficients from PO.DAAC <https://grace.jpl.nasa.gov/data/get-data/geocenter/>`_ [Sun2016]_
      * ``'SLR'``: `Satellite laser ranging coefficients from CSR <ftp://ftp.csr.utexas.edu/pub/slr/geocenter/>`_ [Cheng2013]_
      * ``'SLF'``: `GRACE/GRACE-FO coefficients from Sutterley and Velicogna <https://doi.org/10.6084/m9.figshare.7388540>`_ [Sutterley2019]_
      * ``'Swenson'``: GRACE-derived coefficients from Sean Swenson [Swenson2008]_
      * ``'GFZ'``: `GRACE/GRACE-FO coefficients from GFZ GravIS <http://gravis.gfz-potsdam.de/corrections>`_

Keyword arguments
#################

- ``MMAX``: Upper bound of Spherical Harmonic Orders
- ``SLR_21``: Replaces C21/S21 with values from Satellite Laser Ranging (SLR)

    * ``None``: use original values
    * ``'CSR'``: use values from CSR
    * ``'GFZ'``: use values from GFZ GravIS
    * ``'GSFC'``: use values from GSFC
- ``SLR_22``: Replaces C22/S22 with values from Satellite Laser Ranging (SLR)

    * ``None``: use original values
    * ``'CSR'``: use values from CSR
- ``SLR_C30``: Replaces C30 with values from Satellite Laser Ranging (SLR)

    * ``None``: use original values
    * ``'CSR'``: use values from CSR (5x5 with 6,1)
    * ``'GFZ'``: use values from GFZ GravIS
    * ``'GSFC'``: use values from GSFC (TN-14)
- ``SLR_C50``: Replaces C50 with values from Satellite Laser Ranging (SLR)

    * ``None``: use original values
    * ``'CSR'``: use values from CSR (5x5 with 6,1)
    * ``'GSFC'``: use values from GSFC
- ``POLE_TIDE``: correct GSM data for pole tide drift
- ``ATM``: correct data with ECMWF "jump" corrections GAE, GAF and GAG
- ``MODEL_DEG1``: least-squares model missing degree 1 coefficients
- ``DEG1_GIA``: GIA-correction used when calculating degree 1 coefficients

Returns
#######

- ``clm``: GRACE/GRACE-FO cosine spherical harmonics to degree/order ``LMAX`` and ``MMAX``
- ``slm``: GRACE/GRACE-FO sine spherical harmonics to degree/order ``LMAX`` and ``MMAX``
- ``time``: time of each GRACE/GRACE-FO measurement (mid-month)
- ``month``: GRACE/GRACE-FO months of input datasets
- ``l``: spherical harmonic degree to ``LMAX``
- ``m``: spherical harmonic order to ``MMAX``
- ``title``: string denoting low degree zonals replacement, geocenter usage and corrections
- ``directory``: directory of exact GRACE/GRACE-FO product

References
##########

.. [Cheng2013] M. Cheng, "Geocenter Variations from Analysis of SLR Data", *Reference Frames for Applications in Geosciences*, 19--25, (2013). `doi: 10.1007/978-3-642-32998-2_4 <https://doi.org/10.1007/978-3-642-32998-2_4>`_

.. [Fagiolini2015] E. Fagiolini, F. Flechtner, M. Horwath, and H. Dobslaw, "Correction of inconsistencies in ECMWF's operational analysis data during de-aliasing of GRACE gravity models", *Geophysical Journal International*, 202(3), 2150--2158, (2015). `doi: 10.1093/gji/ggv276 <https://doi.org/10.1093/gji/ggv276>`_

.. [Sun2016] Y. Sun, P. Ditmar, and R. Riva, "Observed changes in the Earth's dynamic oblateness from GRACE data and geophysical models", *Journal of Geodesy*, 90(1), 81--89, (2016). `doi: 10.1007/s00190-015-0852-y <https://doi.org/10.1007/s00190-015-0852-y>`_

.. [Sutterley2019] T. C. Sutterley and I. Velicogna, "Improved Estimates of Geocenter Variability from Time-Variable Gravity and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_

.. [Wahr2015] J. Wahr, R. S. Nerem, and S. V. Bettadpur, "The pole tide and its effect on GRACE time‐variable gravity measurements: Implications for estimates of surface mass variations". *Journal of Geophysical Research: Solid Earth*, 120, 4597--4615. `doi: 10.1002/2015JB011986 <https://doi.org/10.1002/2015JB011986>`_
