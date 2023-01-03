#!/usr/bin/env python
u"""
amplitude.py
Written by Tyler Sutterley (01/2023)

Calculate the amplitude and phase of a harmonic function from calculated
    sine and cosine of a series of measurements

CALLING SEQUENCE:
    ampl,ph = gravity_toolkit.time_series.amplitude(bsin, bcos)

INPUTS:
    bsin: amplitude of the calculated sine values
    bcos: amplitude of the calculated cosine values

OUTPUTS:
    ampl: amplitude from the harmonic functions
    ph: phase from the harmonic functions in degrees

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 01/2023: refactored time series analysis functions
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 07/2020: added function docstrings
    Updated 10/2019: output both amplitude and phase
    Updated 05/2013: converted to python
    Written 07/2012:
"""
import numpy as np

def amplitude(bsin, bcos):
    """
    Calculate the amplitude and phase of a harmonic function

    Parameters
    ----------
    bsin: float
        amplitude of the calculated sine values
    bcos: float
        amplitude of the calculated cosine values

    Returns
    -------
    ampl: float
        amplitude from the harmonic functions
    ph: float
        phase from the harmonic functions in degrees
    """
    ampl = np.sqrt(bsin**2.0 + bcos**2.0)
    ph = 180.0*np.arctan2(bcos, bsin)/np.pi
    return (ampl,ph)
