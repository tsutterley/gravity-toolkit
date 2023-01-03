#!/usr/bin/env python
u"""
tssmooth.py
Written by Tyler Sutterley (04/2022)

Computes a moving average of a time-series using three possible routines:
    1) centered moving average
    2) 13-month Loess filter (default)
    3) 13-month Loess filter weighted and outputs for all dates

Note: due to the missing months in the GRACE/GRACE-FO time series,
    a standard moving average will have problems if the
    missing months are not interpolated.

CALLING SEQUENCE:
    smth = tssmooth(t_in, d_in, HFWTH=6)

INPUTS:
    t_in: input time array
    d_in: input data array

OUTPUTS:
    time: time after removing start and end half-windows
    data: smoothed time-series
    season: seasonal component calculated by the Loess filter
    annual: annual component calculated by the Loess filter
    semiann: semi-annual component calculated by the Loess filter
    trend: instantaneous trend calculated by the Loess filter
    error: estimated error of the instantaneous trend
    noise: noise component after removing the Loess trend and seasonal components
    reduce: original time series after removing start and end half-windows

OPTIONS:
    MOVING: calculates centered moving average using mean of window
        mean of: (January up to December) and (February up to January)
    WEIGHT: smoothing algorithm that backward models dates before
        half-width and forward models dates after half-width
        0: use unweighted Loess filter
        1: use linear weights with Loess filter
        2: use gaussian weights with Loess filter
    HFWTH: half-width of the moving average (default = 6 for 13-month Loess)
    DATA_ERR: input error for known and equal errors (single value)
    STDEV: standard deviation of output error
    CONF: confidence interval of output error (default is for 95%)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 07/2020: added function docstrings
    Updated 10/2019: changing Y/N flags to True/False. output amplitudes
    Updated 09/2019: calculate and output annual and semi-annual phase
    Updated 08/2018: use implicit import of scipy stats and special
    Updated 09/2017: using rcond=-1 in numpy least-squares algorithms
    Updated 09/2014: added output for the reduced time-series
    Updated 06/2014: added parameter DATA_ERR for known and equal errors
    Updated 03/2014: created a new smoothing algorithm
        similar to Loess-type least-squares algorithm but has
        backward models dates before HFWTH and forward models dates after
        if all dates are available will be centrally weighted
        need to update to include error and should find a reference
            as i came up with this algorithm at 4:30am
    Updated 02/2014: minor update to if statements
    Updated 09/2013: switched MOVING flag to Y/N
        Minor change: P_cons, and P_lin changed to P_x0, and P_x1
    Updated 06/2013: added error for instantaneous trend
    Updated 04/2013: converted to python and added more outputs
    Updated 02/2013: using a centered moving average
        added seasonal option to compute the smooth seasonal variation
            calculated by the loess filter program.
        added option to compute the noise component after removing the
        smoothed trend and the seasonal component
    Updated 03/2012: added Loess smoothing following Velicogna (2009)
    Written 12/2011
"""
import warnings
import gravity_toolkit.time_series

def tssmooth(*args, **kwargs):
    """
    Computes the moving average of a time-series

        1) centered moving average
        2) 13-month Loess filter [Velicogna2009]_
        3) 13-month Loess filter weighted and outputs for all dates

    Parameters
    ----------
    t_in: float
        input time array
    d_in: float
        input data array
    HFWTH: int
        half-width of the moving average
    MOVING: bool, default False
        calculates centered moving average using mean of window
    WEIGHT: smoothing algorithm that backward models dates before
        half-width and forward models dates after half-width

            - ``0``: use unweighted Loess filter
            - ``1``: use linear weights with Loess filter
            - ``2``: use gaussian weights with Loess filter
    DATA_ERR: float or list
        input error for known and equal errors
    STDEV: float, default 0
        Standard deviation of output error
    CONF: float, default 0
        Confidence interval of output error

    Returns
    -------
    time: float
        time after removing start and end half-windows
    data: float
        smoothed time-series
    season: float
        seasonal component calculated by the Loess filter
    annual: float
        annual component calculated by the Loess filter
    semiann: float
        semi-annual component calculated by the Loess filter
    trend: float
        instantaneous trend calculated by the Loess filter
    error: float
        estimated error of the instantaneous trend
    noise: float
        noise component after removing the Loess trend and seasonal components
    reduce: float
        original time series after removing start and end half-windows

    References
    ----------
    .. [Velicogna2009] I. Velicogna, "Increasing rates of ice mass loss
        from the Greenland and Antarctic ice sheets revealed by GRACE",
        *Geophysical Research Letters*, 36(L19503),
        `doi: 10.1029/2009GL040222 <https://doi.org/10.1029/2009GL040222>`_
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.time_series instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.time_series.smooth(*args, **kwargs)
