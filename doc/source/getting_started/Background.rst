==========
Background
==========


This documentation is intended to explain how to compute spatial and time series
estimates using GRACE/GRACE-FO time-variable gravity measurements.
``gravity-toolkit`` is a Python-based geophysical software that reads
GRACE/GRACE-FO time-variable gravity solutions for estimating regional mass change.
A suite of geophysical corrections can be applied to the gravity solutions to
optimize the GRACE/GRACE-FO data for particular applications.
This software was developed with the goal of supporting science applications for
time-variable gravity.
``gravity-toolkit`` provides data access utilities for ascii, netCDF4, HDF5 and gfc file formats.
``gravity-toolkit`` also provides some very high-level plotting programs through the
use of `Jupyter Notebooks <../user_guide/Examples.html>`_.

.. graphviz::
    :caption: Data Processing Framework
    :align: center

    digraph {
        G [label="GRACE/GRACE-FO\ntime-variable gravity"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        A [label="Non-tidal Ocean and\nAtmospheric Variation"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        I [label="Glacial Isostatic\nAdjustment"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        W [label="Terrestrial Water\nStorage"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        R [label="gravity-toolkit"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Spatial-Maps.html"]
        T [label="Time Series Analysis"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Time-Series-Analysis.html"]
        D [label="Geocenter Variation"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Geocenter-Variations.html"]
        G -> R [arrowsize=0.8]
        A -> R [arrowsize=0.8]
        I -> R [arrowsize=0.8]
        W -> R [arrowsize=0.8]
        R -> S [arrowsize=0.8]
        R -> T [arrowsize=0.8]
        R -> D [arrowsize=0.8]
    }

Measurement Principle
#####################

GRACE and the GRACE Follow-on (GRACE-FO) missions each consist of twin satellites in similar low Earth orbits [Tapley2019]_.
The primary and secondary instrumentation onboard the GRACE/GRACE-FO satellites are the ranging instrument
(GRACE has a microwave ranging instrument, GRACE-FO has both a microwave ranging instrument and a laser interferometer),
the global positioning system (GPS), the accelerometers and the star cameras.
Data from these instruments are combined to estimate the distance between the two satellites,
the positions of the satellites in space, the pointing vector of the satellites and any non-gravitational
accelerations the satellites experience.

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
There are programs within this repository that can sync with both of these data archives:
``podaac_cumulus.py`` for `PO.DAAC AWS <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/podaac_cumulus.py>`_ and
``gfz_isdc_grace_ftp.py`` for the `GFZ ISDC <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/gfz_isdc_grace_ftp.py>`_.

Geoid Height
############

The Level-2 spherical harmonic product of GRACE and GRACE-FO provides monthly
estimates of the Earth's gravitational field.
The Earth's gravitational field varies in time as masses on and within the
Earth move and are exchanged between components of the Earth system [Wahr1998]_.
The instantaneous shape of the Earth's gravitational field can be described
in terms of an equipotential surface, a surface of constant potential energy
where the gravitational potential is constant [HofmannWellenhof2006]_.
The Earth's geoid is the equipotential surface that coincides with global mean
sea level if the oceans were at rest [HofmannWellenhof2006]_ [Wahr1998]_.
The distance between the geoid and an Earth reference ellipsoid is the
geoid height (:math:`N`), or the geoidal undulation [HofmannWellenhof2006]_.

.. figure:: ../_assets/geoid_height.svg
    :width: 400
    :align: center

    Relationship between ellipsoid height, geoid height, and topographic height [NRC2010]_

In spherical coordinates, the change in the height of the geoid,
:math:`\Delta N(\theta,\phi)`, at colatitude :math:`\theta` and longitude :math:`\phi`,
can be estimated from a series of spherical harmonics as:

.. math::
    :label: 1

    \Delta N(\theta,\phi) = a\sum_{l=1}^{l_{max}}\sum_{m=0}^lP_{lm}(\cos\theta)\left[\Delta C_{lm}\cos{m\phi} + \Delta S_{lm}\sin{m\phi}\right]

where :math:`a` is the average radius of the Earth,
:math:`P_{lm}(\cos\theta)` are the fully-normalized Legendre polynomials of degree :math:`l` and order :math:`m` for the cosine of colatitude :math:`\theta`, and
:math:`\Delta C_{lm}`, :math:`\Delta S_{lm}` are the changes in the cosine and sine spherical harmonics of degree :math:`l` and order :math:`m` [Chao1987]_.

Surface Mass Density
####################

The radial component of a density change within the Earth cannot be uniquely
determined using satellite gravity observations alone [Wahr1998]_.
However, fluctuations in water storage and transport can be assumed to be largely
concentrated within a thin layer near the Earth's surface [Wahr1998]_.
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

References
##########

.. [Chao1987] B. F. Chao and R. S. Gross, "Changes in the Earth's rotation and low-degree gravitational field induced by earthquakes", *Geophysical Journal International*, 91(3), 569--596 (1987). `doi: 10.1111/j.1365-246X.1987.tb01659.x <https://doi.org/10.1111/j.1365-246X.1987.tb01659.x>`_

.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_

.. [NRC2010] National Research Council. *Satellite Gravity and the Geosphere: Contributions to the Study of the Solid Earth and Its Fluid Envelopes*. The National Academies Press, Washington, DC, 1997. ISBN 978-0-309-05792-9. `doi: 10.17226/5767 <https://doi.org/10.17226/5767>`_

.. [Tapley2019] B. D. Tapley, M. M. Watkins, F. Flechtner et al. "Contributions of GRACE to understanding climate change", *Nature Climate Change*, 9, 358--369 (2019). `doi: 10.1038/s41558-019-0456-2 <https://doi.org/10.1038/s41558-019-0456-2>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. |rarr|    unicode:: U+2192 .. RIGHTWARDS ARROW
