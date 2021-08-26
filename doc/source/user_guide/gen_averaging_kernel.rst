=======================
gen_averaging_kernel.py
=======================

- Generates averaging kernel coefficients which minimize the total error [Swenson2002]_
- Uses a normalized version of Christopher Jekeli's Gaussian averaging function [Jekeli1981]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_averaging_kernel import gen_averaging_kernel
    Wlms = gen_averaging_kernel(gclm,gslm,eclm,eslm,sigma,hw,UNITS=0,LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_averaging_kernel.py

Arguments
#########

1. ``gclm``: cosine spherical harmonics of exact averaging kernel
2. ``gslm``: sine spherical harmonics of exact averaging kernel
3. ``eclm``: measurement error in the cosine harmonics
4. ``eslm``: measurement error in the sine harmonics
5. ``sigma``: variance of the surface mass signal
6. ``hw``: Gaussian radius of the kernel in kilometers

Keyword arguments
#################

- ``LMAX``: Upper bound of Spherical Harmonic Degrees
- ``MMAX``: Upper bound of Spherical Harmonic Orders
- ``UNITS``: Units of input spherical harmonics

   * ``0`` fully-normalized
   * ``1`` cm of water thickness
- ``LOVE``: input load Love numbers up to degree of truncation (``hl``, ``kl``, ``ll``)

Returns
#######

- ``clm``: cosine coefficients of the averaging kernel
- ``slm``: sine coefficients of the averaging kernel

References
##########

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Swenson2002] S. Swenson and J. Wahr, "Methods for inferring regional surface‐mass anomalies from Gravity Recovery and Climate Experiment (GRACE) measurements of time‐variable gravity", *Journal of Geophysical Research: Solid Earth*, 107(B9), 2193, (2002). `doi: 10.1029/2001JB000576 <https://doi.org/10.1029/2001JB000576>`_
