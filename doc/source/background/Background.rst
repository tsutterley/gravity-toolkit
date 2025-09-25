==========
Background
==========

Measurement Principle
#####################

GRACE and the GRACE Follow-on (GRACE-FO) missions each consist of twin satellites in similar low Earth orbits :cite:p:`Tapley:2019cm`.
The primary and secondary instrumentation onboard the GRACE/GRACE-FO satellites are the ranging instrument
(GRACE has a microwave ranging instrument, GRACE-FO has both a microwave ranging instrument and a laser interferometer),
the global positioning system (GPS), the accelerometers and the star cameras.
Data from these instruments are combined to estimate the distance between the two satellites,
the positions of the satellites in space, the pointing vector of the satellites and any non-gravitational
accelerations the satellites experience.

.. admonition:: The Big Idea

    GRACE/GRACE-FO senses changes in gravity by measuring the change in distance between the two twin satellites:

    1) As the satellites approach a mass anomaly: leading satellite "feels" a greater gravitational attraction and accelerates |rarr| **distance increases**
    2) As the trailing satellite approaches: greater gravitational attraction |rarr| accelerated by the mass anomaly |rarr| **distance decreases**
    3) Leading satellite passes the anomaly: gravitational attraction pulls backwards |rarr| decelerated by the mass anomaly |rarr| **distance decreases**
    4) When the trailing satellite passes the anomaly and leading satellite is far from the anomaly: trailing satellite decelerated by mass anomaly |rarr| **distance increases back to standard separation**

All the onboard measurements are combined with estimates of the background gravity field, atmospheric and oceanic variability,
and tides to create the `Level-2 spherical harmonic product of GRACE and GRACE-FO`__.

.. __: https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/GRACE-FO_L2_UserHandbook.pdf

Data Processing
###############

There are three main processing centers that create the Level-2 spherical harmonic data as part of the GRACE/GRACE-FO Science Data System (SDS):
the `University of Texas Center for Space Research (CSR) <http://www2.csr.utexas.edu/grace/>`_,
the `German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) <https://www.gfz-potsdam.de/en/grace/>`_ and
the `Jet Propulsion Laboratory (JPL) <https://grace.jpl.nasa.gov/>`_.

GRACE/GRACE-FO data is freely available in the US from
the `NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC) <https://podaac.jpl.nasa.gov/grace>`_ and
internationally from the `GFZ Information System and Data Center (ISDC) <http://isdc.gfz-potsdam.de/grace-isdc/>`_.

.. tip::
    There are programs within this repository that can sync with both of these data archives:
    ``podaac_cumulus.py`` for `PO.DAAC AWS <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/podaac_cumulus.py>`_ and
    ``gfz_isdc_grace_ftp.py`` for the `GFZ ISDC <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/gfz_isdc_grace_ftp.py>`_.

Geoid Height
############

The Level-2 spherical harmonic product of GRACE and GRACE-FO provides monthly
estimates of the Earth's gravitational field.
The Earth's gravitational field varies in time as masses on and within the
Earth move and are exchanged between components of the Earth system :cite:p:`Wahr:1998hy`.
The instantaneous shape of the Earth's gravitational field can be described
in terms of an equipotential surface, a surface of constant potential energy
where the gravitational potential is constant :cite:p:`HofmannWellenhof:2006hy`.
The Earth's geoid is the equipotential surface that coincides with global mean
sea level if the oceans were at rest :cite:p:`HofmannWellenhof:2006hy,Wahr:1998hy`.
The distance between the geoid and an Earth reference ellipsoid is the
geoid height (:math:`N`), or the geoidal undulation :cite:p:`HofmannWellenhof:2006hy`.

.. figure:: ../_assets/geoid_height.svg
    :width: 400
    :align: center

    Relationship between ellipsoid height, geoid height, and topographic height :cite:p:`NRC:1997ea`

In spherical coordinates, the change in the height of the geoid,
:math:`\Delta N(\theta,\phi)`, at colatitude :math:`\theta` and longitude :math:`\phi`,
can be estimated from a series of spherical harmonics as:

.. math::
    :label: 1

    \Delta N(\theta,\phi) = a\sum_{l=1}^{l_{max}}\sum_{m=0}^lP_{lm}(\cos\theta)\left[\Delta C_{lm}\cos{m\phi} + \Delta S_{lm}\sin{m\phi}\right]

where :math:`a` is the average radius of the Earth,
:math:`P_{lm}(\cos\theta)` are the fully-normalized Legendre polynomials of degree :math:`l` and order :math:`m` for the cosine of colatitude :math:`\theta`, and
:math:`\Delta C_{lm}`, :math:`\Delta S_{lm}` are the changes in the cosine and sine spherical harmonics of degree :math:`l` and order :math:`m` :cite:p:`Chao:1987fq`.

Surface Mass Density
####################

The radial component of a density change within the Earth cannot be uniquely
determined using satellite gravity observations alone :cite:p:`Wahr:1998hy`.
However, fluctuations in water storage and transport can be assumed to be largely
concentrated within a thin layer near the Earth's surface :cite:p:`Wahr:1998hy`.
With this assumption, the Earth's surface mass density
(:math:`\Delta\sigma(\theta,\phi)`), the integral of the density change
(:math:`\Delta\rho(r,\theta,\phi)`) through the thin surface layer,
can be estimated as the following:

.. math::
    :label: 2

    \Delta\sigma(\theta,\phi) = \frac{a\rho_{ave}}{3}\sum_{l=0}^{l_{max}}\sum_{m=0}^l\frac{2l+1}{1+k_l}P_{lm}(\cos\theta)\left[\Delta C_{lm}\cos{m\phi} + \Delta S_{lm}\sin{m\phi}\right]

where :math:`\rho_{ave}` is the average density of the Earth, and
:math:`k_l` is the gravitational potential load Love number of degree :math:`l`.
Using this assumption, solid Earth variations occurring outside of this
thin layer, such as Glacial Isostatic Adjustment (GIA) effects,
must be independently estimated and removed.

.. |rarr|    unicode:: U+2192 .. RIGHTWARDS ARROW
