=================
savitzky_golay.py
=================

- Smooth and optionally differentiate data of non-uniform sampling with a Savitzky-Golay filter [Savitzky1964]_
- A type of low-pass filter, particularly suited for smoothing noisy data
- For each point makes a least-square fit of a polynomial over an odd-sized window centered at the point

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.savitzky_golay import savitzky_golay
    sg = savitzky_golay(t_in, d_in, WINDOW=13, ORDER=2)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/savitzky_golay.py

Arguments
#########

- ``t_in``: input time array
- ``d_in``: input data array

Keyword arguments
#################

- ``WINDOW``: length of the window

    * Must be an odd integer
- ``ORDER``: order of the polynomial used in the filtering

    * Must be less than (window_size - 1)
- ``DERIV``: order of the derivative to compute
- ``RATE``: scaling factor for output data and error
- ``DATA_ERR``: estimated data error of known and equal value

Returns
#######

- ``time``: time points for window
- ``data``: smoothed time-series (or n-th derivative)
- ``error``: estimated error at time points

References
##########

.. [Savitzky1964] A. Savitzky, M. J. E. Golay, "Smoothing and Differentiation of Data by Simplified Least Squares Procedures". *Analytical Chemistry*, 36(8), 1627--1639, (1964).
