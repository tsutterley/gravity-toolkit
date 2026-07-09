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

.. |rarr|    unicode:: U+2192 .. RIGHTWARDS ARROW
