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
far-field signals leaking into each regional estimate :cite:p:`Velicogna:2009ft`.

``calc_degree_one.py`` calculates coefficients of degree one by combining
GRACE/GRACE-FO spherical harmonic products with estimates of
ocean bottom pressure (OBP) following :cite:p:`Swenson:2008cr` and :cite:p:`Sutterley:2019bx`.
The method assumes that the change in global surface mass density,
:math:`\Delta\sigma(\theta,\phi)`, can be separated into individual
land and ocean components using a land-function
:math:`\vartheta(\theta,\phi)` :cite:p:`Swenson:2008cr`.

.. math::
    :label: 4

		\Delta\sigma(\theta,\phi) &= \Delta\sigma_{land}(\theta,\phi) + \Delta\sigma_{ocean}(\theta,\phi)\\
		\Delta\sigma_{ocean}(\theta,\phi) &= \vartheta(\theta,\phi)~\Delta\sigma(\theta,\phi)

The oceanic components of the change in degree one spherical harmonics
(:math:`\Delta C^{ocean}_{10}`, :math:`\Delta C^{ocean}_{11}`, and :math:`\Delta S^{ocean}_{11}`)
can then be calculated from the changes in ocean mass,
:math:`\Delta\sigma_{ocean}(\theta,\phi)` :cite:p:`Swenson:2008cr` :cite:p:`Wahr:1998hy`.
If the oceanic contributions to degree one variability
(:math:`\Delta C^{ocean}_{10}`, :math:`\Delta C^{ocean}_{11}`, and :math:`\Delta S^{ocean}_{11}`)
can be estimated from an ocean model, then the unknown complete degree one terms
(:math:`\Delta C_{10}`, :math:`\Delta C_{11}`, and :math:`\Delta S_{11}`) can be
calculated from the residual between the oceanic degree one terms and the
measured mass change over the ocean calculated using all other degrees of
the global spherical harmonics from GRACE/GRACE-FO :cite:p:`Swenson:2008cr` :cite:p:`Sutterley:2019bx`.

The ``calc_degree_one.py`` program will output geocenter files in ascii format
for each GRACE/GRACE-FO month following :cite:p:`Sutterley:2019bx`.
Uncertainties in geocenter due to a combination of error sources can be
estimated using the  ``monte_carlo_degree_one.py`` program.

Load Love Numbers
#################

The degree one Love number of gravitational potential :math:`k_1` is defined so
that the degree one terms describe the offset between the center of mass (CM)
of the combined surface mass and deformed solid Earth, and the center of figure (CF)
of the deformed solid Earth surface :cite:p:`Trupin:1992kp` :cite:p:`Blewitt:2003bz`.
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
