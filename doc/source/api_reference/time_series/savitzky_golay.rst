==========================
time_series.savitzky_golay
==========================

- Smooth and optionally differentiate data of non-uniform sampling with a Savitzky-Golay filter
- A type of low-pass filter, particularly suited for smoothing noisy data
- For each point makes a least-square fit of a polynomial over an odd-sized window centered at the point

Calling Sequence
################

.. code-block:: python

    import gravity_toolkit.time_series
    sg = gravity_toolkit.time_series.savitzky_golay(t_in, d_in, WINDOW=13, ORDER=2)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/time_series/savitzky_golay.py

.. autofunction:: gravity_toolkit.time_series.savitzky_golay
