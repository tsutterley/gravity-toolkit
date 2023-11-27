====================
Geocenter Variations
====================

Variations in the Earth's geocenter reflect the largest scale
variability of mass within the Earth system, and are essential
inclusions for the complete recovery of surface mass change from
time-variable gravity.
The Earth's geocenter is the translation between the Earth's
center of mass (CM) and center of figure (CF) reference frames.
Measurements of time-variable gravity from GRACE and GRACE Follow-On
(GRACE-FO) are set in an instantaneous center of mass (CM) reference frame.
For most science applications of time-variable gravity, a coordinate
system with an origin coinciding with the Earth's center of figure
(CF) is required.
Geocenter variations are represented by the degree one spherical harmonic terms.
The exclusion of degree one terms can have a significant impact on estimates
of ocean mass, ice sheet mass change, and terrestrial hydrology due to
far-field signals leaking into each regional estimate [Velicogna2009]_.

``calc_degree_one.py`` calculates coefficients of degree one by combining
GRACE/GRACE-FO spherical harmonic products with estimates of
ocean bottom pressure (OBP) following [Swenson2008]_ and [Sutterley2019]_.
The method assumes that the change in global surface mass density,
:math:`\Delta\sigma(\theta,\phi)`, can be separated into individual
land and ocean components using a land-function
:math:`\vartheta(\theta,\phi)` [Swenson2008]_.

.. math::
    :label: 4

		\Delta\sigma(\theta,\phi) &= \Delta\sigma_{land}(\theta,\phi) + \Delta\sigma_{ocean}(\theta,\phi)\\
		\Delta\sigma_{ocean}(\theta,\phi) &= \vartheta(\theta,\phi)~\Delta\sigma(\theta,\phi)

The oceanic components of the change in degree one spherical harmonics
(:math:`\Delta C^{ocean}_{10}`, :math:`\Delta C^{ocean}_{11}`, and :math:`\Delta S^{ocean}_{11}`)
can then be calculated from the changes in ocean mass,
:math:`\Delta\sigma_{ocean}(\theta,\phi)` [Swenson2008]_ [Wahr1998]_.
If the oceanic contributions to degree one variability
(:math:`\Delta C^{ocean}_{10}`, :math:`\Delta C^{ocean}_{11}`, and :math:`\Delta S^{ocean}_{11}`)
can be estimated from an ocean model, then the unknown complete degree one terms
(:math:`\Delta C_{10}`, :math:`\Delta C_{11}`, and :math:`\Delta S_{11}`) can be
calculated from the residual between the oceanic degree one terms and the
measured mass change over the ocean calculated using all other degrees of
the global spherical harmonics from GRACE/GRACE-FO [Swenson2008]_ [Sutterley2019]_.

The ``calc_degree_one.py`` program will output geocenter files in ascii format
for each GRACE/GRACE-FO month following [Sutterley2019]_.
Uncertainties in geocenter due to a combination of error sources can be
estimated using the  ``monte_carlo_degree_one.py`` program.

Load Love Numbers
#################

The degree one Love number of gravitational potential :math:`k_1` is defined so
that the degree one terms describe the offset between the center of mass (CM)
of the combined surface mass and deformed solid Earth, and the center of figure (CF)
of the deformed solid Earth surface [Trupin1992]_ [Blewett2003]_.
For the CF coordinate system, this means

.. math::
    :label: 5

	k_1 = -(h_1 + 2\ell_1)/3

where :math:`h_1` and :math:`\ell_1` are the degree one vertical and
horizontal displacement Love numbers.

Geocenter and Degree One
########################

Fully-normalized degree one variations can be converted to
cartesian geocenter variations using the following relation:

.. math::
    :label: 6

	\Delta X &= a\sqrt{3}~\Delta C_{11} \\
	\Delta Y &= a\sqrt{3}~\Delta S_{11} \\
	\Delta Z &= a\sqrt{3}~\Delta C_{10}


The ``geocenter`` class has utilities for converting between
spherical harmonics and geocenter variation along with
readers for different geocenter datasets.

References
##########

.. [Blewett2003] G. Blewitt, "Self-consistency in reference frames, geocenter definition, and surface loading of the solid Earth", *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003). `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

.. [Cheng2013] M. Cheng, "Geocenter Variations from Analysis of SLR Data", *Reference Frames for Applications in Geosciences*, 19--25, (2013). `doi: 10.1007/978-3-642-32998-2_4 <https://doi.org/10.1007/978-3-642-32998-2_4>`_

.. [Dziewonski1981] A. M. Dziewonski and D. L. Anderson, "Preliminary reference Earth model", *Physics of the Earth and Planetary Interiors*, 25(4), 297--356, (1981). `doi: 10.1016/0031-9201(81)90046-7 <https://doi.org/10.1016/0031-9201(81)90046-7>`_

.. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads", *Reviews of Geophysics*, 10(3), 761--797, (1972). `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_

.. [Sutterley2019] T. C. Sutterley and I. Velicogna, "Improved Estimates of Geocenter Variability from Time-Variable Gravity and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_

.. [Trupin1992] A. S. Trupin, M. F. Meier, and J. Wahr, "Effect of melting glaciers on the Earth's rotation and gravitational field: 1965--1984", *Geophysical Journal International*, 108(1), (1992). `doi: 10.1111/j.1365-246X.1992.tb00835.x <https://doi.org/10.1111/j.1365-246X.1992.tb00835.x>`_

.. [Velicogna2009] I. Velicogna, "Increasing rates of ice mass loss from the Greenland and Antarctic ice sheets revealed by GRACE", *Geophysical Research Letters*, 36(L19503), (2009). `doi: 10.1029/2009GL040222 <https://doi.org/10.1029/2009GL040222>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wahr2006] J. Wahr, S. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates", Geophysical Research Letters, 33(L06401), (2006). `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_
