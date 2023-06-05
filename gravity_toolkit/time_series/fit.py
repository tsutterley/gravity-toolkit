#!/usr/bin/env python
u"""
fit.py
Written by Tyler Sutterley (05/2023)
Utilities for fitting time-series data with regression models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Written 05/2023
"""
from __future__ import annotations
import numpy as np

# PURPOSE: build a list of S2 tidal aliasing terms for regression fit
def aliasing_terms(t_in: np.ndarray):
    """
    Create a list of custom design-matrix terms to account for
    S2 tidal aliasing during the GRACE and GRACE-FO periods

    Parameters
    ----------
    t_in: float
        input time array

    Returns
    -------
    TERMS: list
        S2 tidal aliasing terms for GRACE and GRACE-FO missions
    """
    # create output list of fit terms
    TERMS = []
    # number of time points
    nmax = len(t_in)
    # create custom terms for S2 tidal aliasing during GRACE period
    ii, = np.nonzero(t_in[0:nmax] < 2018.0)
    S2sin = np.zeros((nmax))
    S2cos = np.zeros((nmax))
    S2sin[ii] = np.sin(np.pi*t_in[ii]*730.50/161.0)
    S2cos[ii] = np.cos(np.pi*t_in[ii]*730.50/161.0)
    TERMS.append(S2sin)
    TERMS.append(S2cos)
    # create custom terms for S2 tidal aliasing during GRACE-FO period
    ii, = np.nonzero(t_in[0:nmax] >= 2018.0)
    S2sin = np.zeros((nmax))
    S2cos = np.zeros((nmax))
    S2sin[ii] = np.sin(np.pi*t_in[ii]*730.50/161.0)
    S2cos[ii] = np.cos(np.pi*t_in[ii]*730.50/161.0)
    TERMS.append(S2sin)
    TERMS.append(S2cos)
    # return the fit terms
    return TERMS
