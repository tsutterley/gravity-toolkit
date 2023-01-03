#!/usr/bin/env python
u"""
savitzky_golay.py
Written by Tyler Sutterley (01/2023)
Adapted from Numerical Recipes, Third Edition

Smooth and optionally differentiate data of non-uniform sampling
    with a Savitzky-Golay filter

Can preserves the original shape and features of the signal better
    than moving averages techniques

CALLING SEQUENCE:
    y_out = gravity_toolkit.time_series.savitzky_golay(t_in, y_in,
        WINDOW=13, ORDER=2)

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
    Updated 01/2023: refactored time series analysis functions
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: added function docstrings
    Updated 08/2019: importing factorial from scipy.special.factorial
    Updated 11/2018: using future division for python3 compatibility
    Updated 08/2015: changed sys.exit to raise ValueError
    Written 06/2014
"""
from __future__ import print_function, division

import numpy as np
import scipy.special

def savitzky_golay(t_in, y_in, WINDOW=None, ORDER=2, DERIV=0,
    RATE=1, DATA_ERR=0):
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

    # verify that WINDOW is positive, odd and greater than ORDER+1
    if WINDOW is None:
        WINDOW = ORDER + -1*(ORDER % 2) + 3
    if WINDOW % 2 != 1 or WINDOW < 1:
        raise ValueError("WINDOW size must be a positive odd number")
    if WINDOW < ORDER + 2:
        raise ValueError("WINDOW is too small for the polynomials order")
    # remove any singleton dimensions
    t_in = np.squeeze(t_in)
    y_in = np.squeeze(y_in)
    nmax = len(t_in)

    # order range
    order_range = np.arange(ORDER+1)
    # filter half-window
    half_window = (WINDOW - 1) // 2
    # output time-series (removing half-windows on ends)
    t_out = t_in[half_window:nmax-half_window]
    # output smoothed timeseries (or derivative)
    y_out = np.zeros((nmax-2*half_window))
    y_err = np.zeros((nmax-2*half_window))
    for n in range(0, (nmax-(2*half_window))):
        yran = y_in[n + np.arange(0, 2*half_window+1)]
        # Vandermonde matrix for the time-series
        b = np.mat([[(t_in[k]-t_in[n+half_window])**i for i in order_range]
            for k in range(n, n+2*half_window+1)])
        # compute the pseudoinverse of the design matrix
        m=np.linalg.pinv(b).A[DERIV]*RATE**DERIV*scipy.special.factorial(DERIV)
        # pad the signal at the extremes with values taken from the signal
        firstvals = yran[0] - np.abs(yran[1:half_window+1][::-1] - yran[0])
        lastvals = yran[-1] + np.abs(yran[-half_window-1:-1][::-1] - yran[-1])
        yn = np.concatenate((firstvals, yran, lastvals))
        # compute the convolution and use middle value
        y_out[n] = np.convolve(m[::-1], yn, mode='valid')[half_window]
        if (DATA_ERR != 0):
            # if data error is known and of equal value
            P_err = DATA_ERR*np.ones((4*half_window+1))
            # compute the convolution and use middle value
            y_err[n] = np.sqrt(np.convolve(m[::-1]**2, P_err**2,
                mode='valid')[half_window])

    return {'data':y_out, 'error':y_err, 'time':t_out}
