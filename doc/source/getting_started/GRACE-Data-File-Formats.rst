=================
Data File Formats
=================

Product Identifier
##################

GRACE Level-2 products consist of spherical harmonic coefficients of the Earth's gravitational field.
The data files are typically gzipped ascii files with names formatted as the following:
``PID-2_YYYYDOY-yyyydoy_ndays_center_flag_rrrr`` or ``PID-2_YYYYDOY-yyyydoy_mssn_center_flag_rrrr``

 - ``PID`` is a product identification string (for standard products: GSM, GAD, GAC, GAA, GAB)
 - ``-2`` denotes that the data is a GRACE Level-2 product
 - ``YYYYDOY`` denotes the start date (year and day-of-year) of the measurement range
 - ``yyyydoy`` denotes the end date (year and day-of-year) of the measurement range
 - ``ndays`` is the number of calendar days used to produce the monthly estimate
 - ``mssn`` is the mission GRAC for GRACE and GRFO for GRACE Follow-On
 - ``center`` is an institution specific string (UTCSR for CSR, JPLEM for JPL Spherical Harmonics, JPLMSC for JPL Mascons, EIGEN or GFZOP for GFZ)
 - ``flag`` is a 4 character processing center dependent string (CSR denotes the maximum degree and possibly maximum order of the solutions, JPL denotes if the data is an intermediate release, GFZ denotes constrained versus unconstrained solutions).  For Release-6 and beyond this flag denotes the processing.
 - ``rrrr`` is a 4 character release string which is typically a 4 digit number (GFZ datasets can denote intermediate releases in the 4\ :sup:`th` character)

Character Description
#####################

.. table::
    :widths: 10 90

    +-------+-------------------------------------------------+
    |  1st  |                   Description                   |
    +=======+=================================================+
    | **G** | Geopotential coefficients                       |
    +-------+-------------------------------------------------+

.. table::
    :widths: 10 90

    +-------+-------------------------------------------------+
    |  2nd  |                   Description                   |
    +=======+=================================================+
    | **S** | Estimate made from only GRACE/GRACE-FO data     |
    +-------+-------------------------------------------------+
    | **C** | Combination estimate from GRACE/GRACE-FO and    |
    |       | terrestrial gravity information                 |
    +-------+-------------------------------------------------+
    | **E** | Any background model specified as a time-series |
    +-------+-------------------------------------------------+
    | **A** | Average of any background model over a time     |
    |       | period                                          |
    +-------+-------------------------------------------------+

.. table::
    :widths: 10 90

    +-------+-------------------------------------------------+
    |  3rd  |                   Description                   |
    +=======+=================================================+
    | **M** | Estimate of the Static field                    |
    +-------+-------------------------------------------------+
    | **U** | Geopotential estimate relative to the           |
    |       | background gravity model [1]_                   |
    +-------+-------------------------------------------------+
    | **T** | Total background gravity model except for       |
    |       | background static model                         |
    +-------+-------------------------------------------------+
    | **A** | Non-tidal atmosphere                            |
    +-------+-------------------------------------------------+
    | **B** | Non-tidal Oceans                                |
    +-------+-------------------------------------------------+
    | **C** | Combination of non-tidal atmosphere and         |
    |       | ocean                                           |
    +-------+-------------------------------------------------+
    | **D** | Ocean bottom pressure product [2]_              |
    +-------+-------------------------------------------------+


.. [1] Data files for this product also contain records with the epochs and rates used to model secular changes in the background gravity model
.. [2] Summation of non-tidal ocean and atmosphere over the ocean.  Atmospheric pressure values are equal to zero over the land

Summary
#######

- **GSM:** Geopotential coefficients of the static gravity field estimated from GRACE satellite data (produced by all centers).
- **GAA:** Non-tidal atmosphere geopotential coefficients averaged over certain time period (produced by GFZ \& JPL).
- **GAB:** Non-tidal ocean geopotential coefficients averaged over certai n time period (produced by GFZ \& JPL).
- **GAC:** Combination of non-tidal atmosphere and ocean averaged over certain time period (produced by all centers).
- **GAD:** Ocean bottom pressure product.  Combination of surface pressure and ocean pressure over the oceans, zero over the land (produced by all centers).

References
##########

- `Product Description Document for AOD1B Release 06, Technical Document GRACE 327-750, Revision 6.1, GFZ German Research Centre for Geosciences (2017). <https://podaac-tools.jpl.nasa.gov/drive/files/allData/gracefo/docs/AOD1B_PDD_RL06_v6.1.pdf>`_
- `GRACE Follow-On Level-2 Gravity Field Product User Handbook, Technical Document GRACE-FO D-103922, Revision 1.1, NASA Jet Propulsion Laboratory (2019). <https://podaac-tools.jpl.nasa.gov/drive/files/allData/gracefo/docs/GRACE-FO_L2_UserHandbook.pdf>`_
