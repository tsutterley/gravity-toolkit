Background
==========

This documentation is intended to explain how to compute spatial and time series estimates using GRACE/GRACE-FO time-variable gravity measurements.

GRACE and the GRACE Follow-on (GRACE-FO) missions each consist of twin satellites in similar low Earth orbits.
The primary and secondary instrumentation onboard the GRACE/GRACE-FO satellites are the ranging instrument (GRACE has a microwave ranging instrument, GRACE-FO has both a microwave ranging instrument and a laser interferometer), the global positioning system (GPS), the accelerometers and the star cameras.
Data from these instruments are combined to estimate the distance between the two satellites, the positions of the satellites in space, the pointing vector of the satellites and any non-gravitational accelerations the satellites experience.

GRACE/GRACE-FO senses changes in gravity by measuring the change in distance between the two twin satellites:

1) As the satellites approach a mass anomaly: leading satellite "feels" a greater gravitational attraction and accelerates &rarr; **distance increases**
2) As the trailing satellite approaches: greater gravitational attraction &rarr; accelerated by the mass anomaly &rarr; **distance decreases**
3) Leading satellite passes the anomaly: gravitational attraction pulls backwards &rarr; decelerated by the mass anomaly &rarr; **distance decreases**
4) When the trailing satellite passes the anomaly and leading satellite is far from the anomaly: trailing satellite decelerated by mass anomaly &rarr; **distance increases back to standard separation**

All the onboard measurements are combined with estimates of the background gravity field, atmospheric and oceanic variability, and tides to create the [Level-2 spherical harmonic product of GRACE and GRACE-FO](https://podaac-tools.jpl.nasa.gov/drive/files/GeodeticsGravity/gracefo/docs/GRACE-FO_L2-UserHandbook_v1.1.pdf). 

There are three main processing centers that create the Level-2 spherical harmonic data as part of the GRACE/GRACE-FO Science Data System (SDS): the [University of Texas Center for Space Research (CSR)](http://www2.csr.utexas.edu/grace/), the [German Research Centre for Geosciences (GeoForschungsZentrum, GFZ)](https://www.gfz-potsdam.de/en/grace/) and the [Jet Propulsion Laboratory (JPL)](https://grace.jpl.nasa.gov/).  

GRACE/GRACE-FO data is freely available in the US from the [NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace) and internationally from the [GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/).
There are programs within this repository that can sync with both of these data archives [`podaac_grace_sync.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_grace_sync.py) and [`gfz_isdc_grace_ftp.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/gfz_isdc_grace_ftp.py).  