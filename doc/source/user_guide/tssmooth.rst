===========
tssmooth.py
===========

- Computes a moving average of a time-series using three possible routines:

    1. centered moving average
    2. 13-month Loess filter [Velicogna2009]_
    3. weighted 13-month Loess filter

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.tssmooth import tssmooth
    smth = tssmooth(t_in, d_in, HFWTH=6)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/tssmooth.py

Arguments
#########

- ``t_in``: input time array
- ``d_in``: input data array

Keyword arguments
#################

- ``MOVING``: calculates centered moving average using mean of window
- ``WEIGHT``: use smoothing algorithm that backward models dates before half-width and forward models dates after half-width

    * ``0``: use unweighted Loess filter
    * ``1``: use linear weights with Loess filter
    * ``2``: use gaussian weights with Loess filter
- ``HFWTH``: half-width of the moving average
- ``DATA_ERR``: input error for known and equal errors
- ``STDEV``: standard deviation of output error
- ``CONF``: confidence interval of output error

Returns
#######

- ``time``: time after removing start and end half-windows
- ``data``: smoothed time-series
- ``season``: seasonal component calculated by the Loess filter
- ``annual``: annual component calculated by the Loess filter
- ``semiann``: semi-annual component calculated by the Loess filter
- ``trend``: instantaneous trend calculated by the Loess filter
- ``error``: estimated error of the instantaneous trend
- ``noise``: remaining noise after removing the trend and seasonal components
- ``reduce``: original time series after removing start and end half-windows

References
##########

.. [Velicogna2009] I. Velicogna, "Increasing rates of ice mass loss from the Greenland and Antarctic ice sheets revealed by GRACE", *Geophysical Research Letters*, 36(L19503), `doi: 10.1029/2009GL040222 <https://doi.org/10.1029/2009GL040222>`_
