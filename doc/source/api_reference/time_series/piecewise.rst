=====================
time_series.piecewise
=====================

- Fits a synthetic signal to data over a time period by ordinary or weighted least-squares for breakpoint analysis

Calling Sequence
################

.. code-block:: python

    import gravity_toolkit.time_series
    tsbeta = gravity_toolkit.time_series.piecewise(t_in, d_in, BREAKPOINT=len(t_in)//2, CYCLES=[0.5,1.0])

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/time_series/piecewise.py

.. autofunction:: gravity_toolkit.time_series.piecewise
