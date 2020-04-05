========
units.py
========

Class for converting GRACE/GRACE-FO Level-2 data to specific units

Calling Sequence
================

Calculating the degree dependent factors for converting harmonic units

.. code-block:: python

    from gravity_toolkit.units import units
    from gravity_toolkit.read_love_numbers import read_love_numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, REFERENCE='CF')
    dfactor = units(lmax=lmax).harmonic(hl,kl,ll)
    to_cmwe = dfactor.cmwe

Calculating the degree dependent factors for converting spatial units

.. code-block:: python

    from gravity_toolkit.units import units
    from gravity_toolkit.read_love_numbers import read_love_numbers
    hl,kl,ll = read_love_numbers(love_numbers_file, REFERENCE='CF')
    dfactor = units(lmax=lmax).spatial(hl,kl,ll)
    from_cmwe = dfactor.cmwe

General Attributes and Methods
==============================

.. attribute:: object.norm

    geodesy normalized spherical harmonics

.. attribute:: object.cmwe

    centimeters water equivalent

.. attribute:: object.mmwe

    millimeters water equivalent

.. attribute:: object.mmGH

    millimeters geoid height

.. attribute:: object.mmCU

    millimeters elastic crustal deformation (uplift)

.. attribute:: object.mmCH

    millimeters elastic crustal deformation (horizontal)

.. attribute:: object.cmVCU

    centimeters viscoelastic crustal uplift from `Wahr et al. (2000)`__

.. __: https://doi.org/10.1029/2000JB900113

.. attribute:: object.mVCU

    meters viscoelastic crustal uplift from `Wahr et al. (2000)`__

.. __: https://doi.org/10.1029/2000JB900113

.. attribute:: object.microGal

    microGal gravity perturbations

.. attribute:: object.mbar

    millibar equivalent surface pressure

.. attribute:: object.Pa

    pascals equivalent surface pressure

.. attribute:: object.a_axis

    semi-major axis of the WGS84 ellipsoid in cm

.. attribute:: object.flat

    flattening of the WGS84 ellipsoid

.. attribute:: object.b_axis

    semi-minor axis of the WGS84 ellipsoid in cm

.. attribute:: object.rad_e

    average radius of the Earth having the same volume as WGS84 in cm

.. attribute:: object.g_wmo

    standard gravitational acceleration in cm/s\ :sup:`2`

.. attribute:: object.rho_e

    average density of the Earth in g/cm\ :sup:`3`

.. method:: object.harmonic(hl, kl, ll)

    Calculates degree dependent factors for converting harmonic units from `Wahr et al. (1998)`__

    Inputs: `hl`, `kl`, `ll` load Love numbers to degree `lmax`

.. __: https://doi.org/10.1029/98JB02844

.. method:: object.spatial(hl, kl, ll)

    Calculates degree dependent factors for converting spatial units from `Wahr et al. (1998)`__

    Inputs: `hl`, `kl`, `ll` load Love numbers to degree `lmax`

.. __: https://doi.org/10.1029/98JB02844
