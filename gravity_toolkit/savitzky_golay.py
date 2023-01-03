#!/usr/bin/env python
u"""
savitzky_golay.py
Written by Tyler Sutterley (04/2022)
Adapted from Numerical Recipes, Third Edition

Smooth and optionally differentiate data of non-uniform sampling
    with a Savitzky-Golay filter

Can preserves the original shape and features of the signal better
    than moving averages techniques

CALLING SEQUENCE:
    y_out = savitzky_golay(t_in, y_in, WINDOW=13, ORDER=2)

INPUTS:
    t_in: input time array
    y_in: input data array

OPTIONS:
    WINDOW: length of the window (such as 13 for annual).
        Must be an odd integer number.
    ORDER: order of the polynomial used in the filtering.
        Must be less than (window_size - 1)
    DERIV: the order of the derivative to compute
        default = 0 means only smoothing
    RATE: scaling factor for output data and error
    DATA_ERR: estimated data error of known and equal value

OUTPUTS:
    data: smoothed time-series (or n-th derivative)
    error: estimated error at time points
    time: time points for window

NOTES:
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over an odd-sized window centered at
    the point.

    A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
        Data by Simplified Least Squares Procedures. Analytical
        Chemistry, 1964, 36 (8), pp 1627-1639.
    Numerical Recipes 3rd Edition: The Art of Scientific Computing
        W.H. Press, S.A. Teukolsky, W. T. Vetterling, B.P. Flannery
        Cambridge University Press

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: added function docstrings
    Updated 08/2019: importing factorial from scipy.special.factorial
    Updated 11/2018: using future division for python3 compatibility
    Updated 08/2015: changed sys.exit to raise ValueError
    Written 06/2014
"""
import warnings
import gravity_toolkit.time_series

def savitzky_golay(*args, **kwargs):
    """
    Smooth and optionally differentiate data with a Savitzky-Golay
    filter [Savitzky1964]_ [Press2007]_

    Parameters
    ----------
    t_in: float
        time array
    y_in: float
        data magnitude array
    WINDOW: int or NoneType, default None
        Length of the window

        Must be an odd integer
    ORDER: int, defualt 2
        Order of the polynomial used in the filtering

        Must be less than (window_size - 1)
    DERIV: int, default 0
        Order of the derivative to compute
    RATE: float, default 1
        Scaling factor for output data and error
    DATA_ERR: float, default 0
        Estimated data error of known and equal value

    Returns
    -------
    data: float
        Smoothed signal (or n-th derivative)
    error: float
        Estimated error at time points
    time: float
        Time points for window

    References
    ----------
    .. [Savitzky1964] A. Savitzky, M. J. E. Golay, "Smoothing and
        Differentiation of Data by Simplified Least Squares Procedures".
        *Analytical Chemistry*, 36(8), 1627--1639, (1964).
    .. [Press2007] *Numerical Recipes 3rd Edition: The Art of Scientific
        Computing*, W.H. Press, S.A. Teukolsky, W. T. Vetterling,
        B.P. Flannery. Cambridge University Press, (2007).
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.time_series instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.time_series.savitzky_golay(*args,**kwargs)
