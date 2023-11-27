#!/usr/bin/env python
u"""
fit.py
Written by Tyler Sutterley (06/2023)
Utilities for fitting time-series data with regression models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 06/2023: made the tidal aliasing period an option
    Written 05/2023
"""
from __future__ import annotations
import numpy as np

# PURPOSE: build a list of tidal aliasing terms for regression fit
def aliasing_terms(t_in: np.ndarray, period=161.0):
    """
    Create a list of custom design-matrix terms to account for
    tidal aliasing during the GRACE and GRACE-FO periods

    Parameters
    ----------
    t_in: float
        input time array
    period: float
        tidal aliasing period for a constituent

    Returns
    -------
    TERMS: list
        tidal aliasing terms for GRACE and GRACE-FO missions
    """
    # create output list of fit terms
    TERMS = []
    # number of time points
    nmax = len(t_in)
    # create custom terms for tidal aliasing during GRACE period
    ii, = np.nonzero(t_in[0:nmax] < 2018.0)
    SIN = np.zeros((nmax))
    COS = np.zeros((nmax))
    SIN[ii] = np.sin(np.pi*t_in[ii]*730.50/period)
    COS[ii] = np.cos(np.pi*t_in[ii]*730.50/period)
    TERMS.append(SIN)
    TERMS.append(COS)
    # create custom terms for tidal aliasing during GRACE-FO period
    ii, = np.nonzero(t_in[0:nmax] >= 2018.0)
    SIN = np.zeros((nmax))
    COS = np.zeros((nmax))
    SIN[ii] = np.sin(np.pi*t_in[ii]*730.50/period)
    COS[ii] = np.cos(np.pi*t_in[ii]*730.50/period)
    TERMS.append(SIN)
    TERMS.append(COS)
    # return the fit terms
    return TERMS
