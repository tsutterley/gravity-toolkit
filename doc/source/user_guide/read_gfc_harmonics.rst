=====================
read_gfc_harmonics.py
=====================

- Reads gfc files and extracts spherical harmonics for SWARM and GRAZ GRACE/GRACE-FO data
- Parses date of GRACE/GRACE-FO data from filename

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
    Ylms = read_gfc_harmonics(input_file, TIDE='tide_free')

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_gfc_harmonics.py

Arguments
#########

1. ``input_file``: full path to gfc spherical harmonic data file

Keyword arguments
#################

- ``TIDE``: `tide system of output geoid <http://mitgcm.org/~mlosch/geoidcookbook/node9.html>`_ [Losch2003]_

    * ``'tide_free'``: no permanent direct and indirect tidal potentials
    * ``'mean_tide'``: permanent tidal potentials (direct and indirect)
    * ``'zero_tide'``: permanent direct tidal potential
- ``FLAG``: string denoting data lines

Returns
#######

- ``time``: mid-month date in decimal form
- ``start``: Julian dates of the start date
- ``end``: Julian dates of the start date
- ``l``: spherical harmonic degree to maximum degree of model
- ``m``: spherical harmonic order to maximum degree of model
- ``clm``: cosine spherical harmonics of input data
- ``slm``: sine spherical harmonics of input data
- ``eclm``: cosine spherical harmonic standard deviations
- ``eslm``: sine spherical harmonic standard deviations
- ``modelname``: name of the gravity model
- ``earth_gravity_constant``: GM constant of the Earth for the gravity model
- ``radius``: semi-major axis of the Earth for the gravity model
- ``max_degree``: maximum degree and order for the gravity model
- ``errors``: error type of the gravity model
- ``norm``: normalization of the spherical harmonics
- ``tide_system``: tide system of gravity model (mean_tide, zero_tide, tide_free)

References
##########

.. [Losch2003] M. Losch and V. Seufer, "How to Compute Geoid Undulations (Geoid Height Relative to a Given Reference Ellipsoid) from Spherical Harmonic Coefficients for Satellite Altimetry Applications", (2003). `eprint ID: 11802 <http://mitgcm.org/~mlosch/geoidcookbook.pdf>`_
