================
gauss_weights.py
================

- Computes the Gaussian weights as a function of degree
- A normalized version of Christopher Jekeli's Gaussian averaging function [Jekeli1981]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gauss_weights import gauss_weights
    wl = 2.0*np.pi*gauss_weights(hw,LMAX)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gauss_weights.py

Arguments
#########

1. ``hw``: Gaussian smoothing radius in km
2. ``LMAX``: Upper bound of Spherical Harmonic Degrees

Returns
#######

- ``wl``: Gaussian weights for each degree ``l``

References
##########

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).
