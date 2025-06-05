#!/usr/bin/env python
u"""
lomb_scargle.py
Written by Tyler Sutterley (06/2024)

Wrapper function for computing Lomb-Scargle periodograms

Probability is calculated using a calculation of the independent
    frequencies from Horne and Baliunas (1986)

INPUTS:
    t_in: input time array
    d_in: input data array

OUTPUTS:
    PowerDensity: normalized spectral power density
    probability: probability of each frequency
    frequency: considered frequencies array (omega/2pi)
    period: considered periods array (1/f)
    peak: period at peak power density
    centroid: centroid of power density and period

OPTIONS:
    NORMALIZE: compute normalized periodogram
    OMEGA: angular frequency range [min, max]
    FREQUENCY: temporal frequency range [min, max]
    PERIOD: range of periods [min, max]
    N: number of frequencies
    p: probability levels for contours

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

REFERENCES:
    N. R. Lomb, "Least-squares frequency analysis of unequally spaced data",
        Astrophysics and Space Science, 39, 447--462, 1976.
        https://doi.org/10.1007/BF00648343
    J. D. Scargle, "Studies in astronomical time series analysis. II -
        Statistical aspects of spectral analysis of unevenly spaced data",
        The Astrophysical Journal, 263, 835--853, 1982.
        https://doi.org/10.1086/160554
    J. H. Horne and S. L. Baliunas, "A Prescription for Period
        Analysis of Unevenly Sampled Time Series",
        The Astrophysical Journal, 302, 757--763, 1986.
        https://doi.org/10.1086/164037

UPDATE HISTORY:
    Updated 06/2024: add function docstrings
    Updated 01/2015: added centroid output
    Written 08/2013
"""
import numpy as np
import scipy.signal

def lomb_scargle(t_in, d_in, **kwargs):
    """
    Computes periodograms for least-squares spectral analysis following
    :cite:p:`Lomb:1976bo,Scargle:1982eu` and computes the
    frequency probabilities following :cite:t:`Horne:1986ds`

    Parameters
    ----------
    t_in: float
        input time array
    d_in: float
        input data array
    NORMALIZE: bool, default False
        Compute normalized periodogram
    OMEGA: list, default []
        Angular frequency range
    FREQUENCY: list, default []
        Temporal frequency range
    PERIOD: list, default []
        Temporal period range
    N: int, default 1000
        Number of frequencies
    p: list, default [0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6]
        Probability levels for contours

    Returns
    -------
    PowerDensity: float
        spectral power density
    probability: float
        probability of each frequency
    frequency: float
        considered frequencies array
    period: float
        periods array
    peak: float
        period at peak power density
    centroid: float
        centroid of power density and period
    """
    # default keyword arguments
    kwargs.setdefault('NORMALIZE', False)
    kwargs.setdefault('OMEGA', [])
    kwargs.setdefault('FREQUENCY', [])
    kwargs.setdefault('PERIOD', [])
    kwargs.setdefault('N', 1000)
    kwargs.setdefault('p', [0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6])

    # remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)

    # number of independent measurements
    nmax = np.count_nonzero(np.isfinite(d_in))
    nyquist = 1.0/(2.0*np.mean(t_in[1:] - t_in[0:-1]))

    # angular frequency range
    if kwargs['OMEGA']:
        OMEGA = np.atleast_1d(kwargs['OMEGA'])
    elif kwargs['FREQUENCY']:
        OMEGA = np.atleast_1d(kwargs['FREQUENCY'])/(2.0*np.pi)
    elif kwargs['PERIOD']:
        OMEGA = (2.0*np.pi)/np.atleast_1d(kwargs['PERIOD'])
    else:
        raise ValueError('Frequency range must be defined')

    # number of angular frequencies
    N = int(kwargs['N'])
    assert len(OMEGA) == 2, 'Angular frequency range must have 2 values'
    # array of angular frequencies
    angular_freq = np.linspace(OMEGA[0], OMEGA[1], N)

    # Estimate of number of independent frequencies in Lomb-Scargle
    # analysis based on sample size.  From Horne and Baliunas,
    # "A Prescription for Period Analysis of Unevenly Sampled Time Series",
    # The Astrophysical Journal, 392: 757-763, 1986.
    independent_freq = np.round(-6.362 + 1.193*nmax + 0.00098*nmax**2)
    # if less than 1 independent frequency: set equal to 1
    independent_freq = np.maximum(independent_freq, 1)

    # scaling the date (t[0] = 0)
    t = t_in - t_in[0]
    # periods and frequencies considered
    frequency = angular_freq/(2.0*np.pi)
    period = 2.0*np.pi/angular_freq
    # scaling the data to be mean 0 with variance 1
    data_norm = (d_in - d_in.mean())/d_in.std()
    # computing the lomb-scargle periodogram
    # "normalized" spectral density refers to variance term in denominator
    # PowerDensity has exponential probability distribution with unit mean
    # can calculate normalized as described in Scipy reference
    PowerDensity = scipy.signal.lombscargle(t, data_norm, angular_freq,
        normalize=kwargs['NORMALIZE'])
    # probability of frequencies (NULL test, significance of peak)
    probability = 1.0 - (1.0-np.exp(-PowerDensity))**independent_freq
    # probability contours
    p = np.atleast_1d(kwargs['p'])
    contour = -np.log(1.0 - (1 - p)**(1.0/independent_freq))
    # period at peak (maximum probability)
    ipeak = np.argmax(PowerDensity)
    peak = period[ipeak]
    # period at signal centroid
    centroid = np.sum(period*PowerDensity)/np.sum(PowerDensity)

    return {'PowerDensity':PowerDensity, 'Probability':probability, 
        'frequency':frequency, 'period':period, 'contour':contour,
        'Nyquist':nyquist, 'peak':peak, 'centroid':centroid}
