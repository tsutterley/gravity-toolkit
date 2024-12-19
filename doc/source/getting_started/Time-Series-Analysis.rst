====================
Time Series Analysis
====================

Least-squares mascons are a method of extracting a regional signal from the
GRACE/GRACE-FO spherical harmonic data.
The procedure was outlined in procedure outlined in
:cite:p:`Tiwari:2009bx` and :cite:p:`Jacob:2012gv`.
Least-squares mascons can be considered a post-processing technique for
analyzing the GRACE/GRACE-FO data.
The technique calculates the scaling factor between an input kernel and the
GRACE/GRACE-FO data for a given month.
In an ideal case, the input kernel is a spatially discrete, uniformly
distributed layer of equivalent water height at colatitudes
:math:`\theta` and longitudes :math:`\phi` :cite:p:`Rowlands:2010hj`.

.. math::
    :label: 3

	k(\theta,\phi) =
	\begin{cases}
		~1~\text{cm w.e.} & \text{if}~(\theta,\phi)~\text{is in region}\\
		~0 & \text{if}~(\theta,\phi)~\text{is not in region}
	\end{cases}

GRACE/GRACE-FO measurements over the region at any given time would be a
scalar multiple of this uniform layer.
For example, for a uniform kernel equivalent to 1 cm w.e.,
if GRACE/GRACE-FO measures 6 cm w.e. in the region, then the scaling factor would be 6.
However, the GRACE/GRACE-FO harmonic solutions are truncated and typically smoothed,
meaning that the sharp 0-to-1 transitions along regional boundaries cannot be
resolved at the inherent GRACE/GRACE-FO resolution :cite:p:`Wahr:1998hy`.
The mascon kernels are instead represented as sets of truncated spherical
harmonics processed in the same manner as the GRACE/GRACE-FO data :cite:p:`Jacob:2012gv`.
In this case, each mascon kernel is a smoothed function with a total
mass equal to the idealized case.

Ideally, the final solution for the recovered mascon mass, is equal
to the true spatial average across the mascon :cite:p:`Jacob:2012gv`.
Misfits in the regression or malformed initial kernels can lead to
the leakage of GRACE/GRACE-FO signal in-between mascons or out of the system.
The least squares mascon technique assumes that the GRACE/GRACE-FO signal are
well represented as scalar multiples of the mascons at any given time.
As the initial mascon parameters are designed with uniform mass distributions,
the GRACE/GRACE-FO anomalies over each mascon must also be uniform to limit
statistical misfit :cite:p:`Jacob:2012gv`.

Hypothetically, we would want as small of kernels as possible to get as
geophysically relevant as possible (perfect fit to region shape).
However, GRACE/GRACE-FO data is limited to a specific spherical harmonic range
with noise increasing with higher degree and order (data-to-noise ratio decreases).
If the kernels are too small, there will be a higher degree of ringing when
truncated (Gibbs phenomenon) and require more information at the higher degrees
and orders to distinguish from the adjacent kernels.
Thus there is a balance between wanting kernels as geophysically relevant as
possible with wanting idealistic kernels which minimize ringing and are resolvable.
Getting the kernels "just right" in order to isolate regions of interest takes some time.

The set of least-squares mascon programs have been used in :cite:p:`Velicogna:2014km`
and other publications for regional time series analysis.
The ``calc_mascon.py`` program additionally calculates the GRACE/GRACE-FO error
harmonics following :cite:p:`Wahr:2006bx`.

The ``calc_mascon.py`` program will output a text file of the time series for each mascon
(format: GRACE/GRACE-FO month, mid-month date in decimal-year format,
estimated monthly mass anomaly [Gt], estimated monthly error [Gt],
mascon area [km\ :sup:`2`]).
