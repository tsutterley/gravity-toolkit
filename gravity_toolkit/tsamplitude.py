#!/usr/bin/env python
u"""
tsamplitude.py
Written by Tyler Sutterley (07/2020)

Calculate the amplitude and phase of a harmonic function from calculated
    sine and cosine of a series of measurements

CALLING SEQUENCE:
    ampl,ph = tsamplitude(bsin, bcos)

INPUTS:
    bsin: amplitude of the calculated sine values
    bcos: amplitude of the calculated cosine values

OUTPUTS:
    ampl: amplitude from the harmonic functions
    ph: phase from the harmonic functions in degrees

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 10/2019: output both amplitude and phase
    Updated 05/2013: converted to python
    Written 07/2012:
"""
import numpy as np

def tsamplitude(bsin, bcos):
    """
    Calculate the amplitude and phase of a harmonic function

    arguments
    ---------
    bsin: amplitude of the calculated sine values
    bcos: amplitude of the calculated cosine values

    Returns
    -------
    ampl: amplitude from the harmonic functions
    ph: phase from the harmonic functions in degrees
    """
    ampl = np.sqrt(bsin**2.0 + bcos**2.0)
    ph = 180.0*np.arctan2(bcos, bsin)/np.pi
    return (ampl,ph)
