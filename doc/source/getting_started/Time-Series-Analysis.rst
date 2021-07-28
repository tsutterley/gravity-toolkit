====================
Time Series Analysis
====================

Least-squares mascons are a method of extracting a regional signal from the
GRACE/GRACE-FO spherical harmonic data.
The procedure was outlined in procedure outlined in
[Tiwari2009]_ and [Jacob2012]_.
Least-squares mascons can be considered a post-processing technique for
analyzing the GRACE/GRACE-FO data.
The technique calculates the scaling factor between an input kernel and the
GRACE/GRACE-FO data for a given month.
For example, for a uniform kernel equivalent to 1 cm w.e.,
if GRACE/GRACE-FO measures 6 cm w.e. in the region, then the scaling factor would be 6.

Ideally we would want as small of kernels as possible to get as geophysically
relevant as possible (perfect fit to region shape).
However, GRACE/GRACE-FO data is limited to a specific spherical harmonic range
with noise increasing with higher degree and order (data-to-noise ratio decreases).
If the kernels are too small, there will be a higher degree of ringing when
truncated (Gibbs phenomenon) and require more information at the higher degrees
and orders to distinguish from the adjacent kernels.
Thus there is a balance between wanting kernels as geophysically relevant as
possible with wanting idealistic kernels which minimize ringing and are resolvable.
Getting the kernels "just right" in order to isolate regions of interest takes some time.

The set of least-squares mascon programs have been used in [Velicogna2014]_
and other publications for regional time series analysis.
The ``calc_mascon.py`` program additionally calculates the GRACE/GRACE-FO error
harmonics following [Wahr2006]_.

The ``calc_mascon.py`` program will output a text file of the time series for each mascon
(format: GRACE/GRACE-FO month, mid-month date in decimal-year format,
estimated monthly mass anomaly [Gt], estimated monthly error [Gt],
mascon area [km\ :sup:`2`]).

References
##########

.. [Jacob2012] T. Jacob, J. Wahr, W. T. Pfeffer, and S. Swenson, "Recent contributions of glaciers and ice caps to sea level rise", *Nature*, 482, 514--518, (2012). `doi: 10.1038/nature10847 <https://doi.org/10.1038/nature10847>`_

.. [Tiwari2009] V. M. Tiwari, J. Wahr, and S. Swenson, "Dwindling groundwater resources in northern India, from satellite gravity observations", *Geophysical Research Letters*, 36(L18401), (2009). `doi: 10.1029/2009GL039401 <https://doi.org/10.1029/2009GL039401>`_

.. [Velicogna2014] I. Velicogna, T. C. Sutterley, and M. R. van den Broeke, "Regional acceleration in ice mass loss from Greenland and Antarctica using GRACE time‐variable gravity data", *Geophysical Research Letters*, 119, 8130--8137, (2014). `doi: 10.1002/2014GL061052 <https://doi.org/10.1002/2014GL061052>`_

.. [Wahr2006] J. Wahr, S. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates", Geophysical Research Letters, 33(L06401), (2006). `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_
