==========
Background
==========


This documentation is intended to explain how to compute spatial and time series
estimates using GRACE/GRACE-FO time-variable gravity measurements.
``read-GRACE-harmonics`` is a Python-based geophysical software that reads
GRACE/GRACE-FO time-variable gravity solutions for estimating regional mass change.
A suite of geophysical corrections can be applied to the gravity solutions to
optimize the GRACE/GRACE-FO data for particular applications.
This software was developed with the goal of supporting science applications for
time-variable gravity.
``read-GRACE-harmonics`` provides data access utilities for ascii, netCDF4, HDF5 and gfc file formats.
``read-GRACE-harmonics`` also provides some very high-level plotting programs through the
use of `Jupyter Notebooks <./Examples.html>`_.

.. graphviz::
    :caption: Data Processing Framework
    :align: center

    digraph {
        G [label="GRACE/GRACE-FO\ntime-variable gravity" shape=box style="filled" color="darkorchid"]
        A [label="Non-tidal Ocean and\nAtmospheric Variation" shape=box style="filled" color="darkorchid"]
        I [label="Glacial Isostatic\nAdjustment" shape=box style="filled" color="darkorchid"]
        W [label="Terrestrial Water\nStorage" shape=box style="filled" color="darkorchid"]
        R [label="read-GRACE-harmonics" shape=box style="filled" color="gray"]
        S [label="Spatial Maps" shape=box style="filled" color="mediumseagreen"]
        T [label="Time Series" shape=box style="filled" color="mediumseagreen"]
        G -> R
        A -> R
        I -> R
        W -> R
        R -> S
        R -> T
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

.. __: https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/gracefo/docs/GRACE-FO_L2-UserHandbook_v1.1.pdf

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
``podaac_grace_sync.py`` for `PO.DAAC <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_grace_sync.py>`_ and
``gfz_isdc_grace_ftp.py`` for the `GFZ ISDC <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/gfz_isdc_grace_ftp.py>`_.

References
##########

.. [Tapley2019] B. D. Tapley, M. M. Watkins, F. Flechtner et al. "Contributions of GRACE to understanding climate change", *Nature Climate Change*, 9, 358--369 (2019). `doi: 10.1038/s41558-019-0456-2 <https://doi.org/10.1038/s41558-019-0456-2>`_

.. |rarr|    unicode:: U+2192 .. RIGHTWARDS ARROW