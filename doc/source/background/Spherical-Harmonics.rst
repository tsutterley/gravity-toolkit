.. _spherical-harmonics:

Spherical Harmonics
===================

Geoid Height
------------

The Level-2 spherical harmonic product of GRACE and GRACE-FO provides monthly estimates of the Earth's gravitational field [see :ref:`fig-sphharm`].
The Earth's gravitational field varies in time as masses on and within the Earth move and are exchanged between components of the Earth system :cite:p:`Wahr:1998hy`.
The instantaneous shape of the Earth's gravitational field can be described in terms of an equipotential surface, a surface of constant potential energy where the gravitational potential is constant :cite:p:`HofmannWellenhof:2006hy`.
The Earth's geoid is the equipotential surface that coincides with global mean sea level if the oceans were at rest :cite:p:`HofmannWellenhof:2006hy,Wahr:1998hy`.
The distance between the geoid and an Earth reference ellipsoid is the geoid height (:math:`N`), or the geoidal undulation :cite:p:`HofmannWellenhof:2006hy`.

.. figure:: ../_assets/geoid_height.svg
    :width: 400
    :align: center

    Relationship between ellipsoid height, geoid height, and topographic height :cite:p:`NRC:1997ea`

In spherical coordinates, the change in the height of the geoid, :math:`\Delta N(\theta,\phi)`, at colatitude :math:`\theta` and longitude :math:`\phi`, can be estimated from a series of spherical harmonics as:

.. math::
    :label: 1

    \Delta N(\theta,\phi) = a\sum_{l=1}^{l_{max}}\sum_{m=0}^lP_{lm}(\cos\theta)\left[\Delta C_{lm}\cos{m\phi} + \Delta S_{lm}\sin{m\phi}\right]

where :math:`a` is the average radius of the Earth, :math:`P_{lm}(\cos\theta)` are the fully-normalized Legendre polynomials of degree :math:`l` and order :math:`m` for the cosine of colatitude :math:`\theta`, and :math:`\Delta C_{lm}`, :math:`\Delta S_{lm}` are the changes in the cosine and sine spherical harmonics of degree :math:`l` and order :math:`m` :cite:p:`Chao:1987fq`.

Surface Mass Density
--------------------

The radial component of a density change within the Earth cannot be uniquely determined using satellite gravity observations alone :cite:p:`Wahr:1998hy`.
However, fluctuations in water storage and transport can be assumed to be largely concentrated within a thin layer near the Earth's surface :cite:p:`Wahr:1998hy`.
With this assumption, the Earth's surface mass density (:math:`\Delta\sigma(\theta,\phi)`), the integral of the density change (:math:`\Delta\rho(r,\theta,\phi)`) through the thin surface layer, can be estimated as the following:

.. math::
    :label: 2

    \Delta\sigma(\theta,\phi) = \frac{a\rho_{ave}}{3}\sum_{l=0}^{l_{max}}\sum_{m=0}^l\frac{2l+1}{1+k_l}P_{lm}(\cos\theta)\left[\Delta C_{lm}\cos{m\phi} + \Delta S_{lm}\sin{m\phi}\right]

where :math:`\rho_{ave}` is the average density of the Earth, and :math:`k_l` is the gravitational potential load Love number of degree :math:`l`.
Using this assumption, solid Earth variations occurring outside of this thin layer, such as Glacial Isostatic Adjustment (GIA) effects, must be independently estimated and removed.

.. _fig-sphharm:

Low-Degree Harmonics
--------------------

.. plot:: ./background/sphharm.py
    :show-source-link: False
    :caption: Spherical harmonics for degrees 1 through 4
    :align: center
