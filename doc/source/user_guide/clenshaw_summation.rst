=====================
clenshaw_summation.py
=====================

- Returns the spatial field for a series of spherical harmonics at a sequence of ungridded points
- Uses a Clenshaw summation to calculate the spherical harmonic summation [Holmes2002]_ [Tscherning1982]_

Calling Sequence
################

.. code-block:: python

   from gravity_toolkit.clenshaw_summation import clenshaw_summation
   spatial = clenshaw_summation(clm,slm,lon,lat,UNITS=1,LMAX=60,LOVE=LOVE)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/clenshaw_summation.py

Arguments
#########

1. ``clm``: cosine spherical harmonic coefficients
2. ``slm``: sine spherical harmonic coefficients
3. ``lon``: longitude of points
4. ``lat``: latitude of points

Keyword arguments
#################

- ``RAD``: Gaussian smoothing radius (km)
- ``UNITS``: output data units

   * ``1`` cm of water thickness
   * ``2`` mm of geoid height
   * ``3`` mm of elastic crustal deformation [Davis2004]_
   * ``4`` microGal gravitational perturbation
   * ``5`` millibar equivalent surface pressure
   * ``6`` cm of viscoelastic crustal uplift (GIA) [Wahr2000]_
- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``LOVE``: input load Love numbers up to degree of truncation (``hl``, ``kl``, ``ll``)
- ``ASTYPE``: floating point precision for calculating Clenshaw summation
- ``SCALE``: scaling factor to prevent underflow in Clenshaw summation

Returns
#######

- ``spatial``: spatial field

Dependencies
############

- ``gauss_weights.py``: Computes the Gaussian weights as a function of degree
- ``units.py``: Class for converting spherical harmonic data to specific units

References
##########

.. [Davis2004] J. L. Davis et al., "Climate‐driven deformation of the solid Earth from GRACE and GPS", *Geophysical Research Letters*, 31(L24605), (2004). `doi: 10.1029/2004GL021435 <https://doi.org/10.1029/2004GL021435>`_

.. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002). `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_

.. [Tscherning1982] C. C. Tscherning and K. Poder, "Some Geodetic Applications of Clenshaw Summation", *Bollettino di Geodesia e Scienze*, 4, 349--375, (1982).

.. [Wahr2000] J. Wahr, D. Wingham, and C. Bentley, "A method of combining ICESat and GRACE satellite data to constrain Antarctic mass balance", *Journal of Geophysical Research: Solid Earth*, 105(B7), 16279--16294, (2000). `doi: 10.1029/2000JB900113 <https://doi.org/10.1029/2000JB900113>`_
