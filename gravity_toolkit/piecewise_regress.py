#!/usr/bin/env python
u"""
piecewise_regress.py
Written by Tyler Sutterley (04/2022)

Fits a synthetic signal to data over a time period by ordinary or weighted
    least-squares for breakpoint analysis

Derivation of Sharp Breakpoint Piecewise Regression:
    https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/02-0472
        y = beta_0 + beta_1*t + e (for x <= alpha)
        y = beta_0 + beta_1*t + beta_2*(t-alpha) + e (for x > alpha)

Fit significance derivations are based on Burnham and Anderson (2002)
    Model Selection and Multimodel Inference

CALLING SEQUENCE:
    pcwbeta = piecewise_regress(tdec,data,CYCLES=[0.5,1.0],BREAKPOINT=ind)

INPUTS:
    t_in: input time array
    d_in: input data array

OUTPUTS:
    beta: regressed coefficients array
    error: regression fit error for each coefficient for an input deviation
        STDEV: standard deviation of output error
        CONF: confidence interval of output error
    std_err: standard error for each coefficient
    R2: coefficient of determination (r**2).
        Proportion of variability accounted by the model
    R2Adj: adjusted r**2. adjusts the r**2 for the number of terms in the model
    MSE: mean square error
    WSSE: Weighted sum of squares error
    NRMSE: normalized root mean square error
    AIC: Akaike information criterion (Second-Order, AICc)
    BIC: Bayesian information criterion (Schwarz criterion)
    model: modeled timeseries
    simple: modeled timeseries without oscillating components
    residual: model residual
    DOF: degrees of freedom
    N: number of terms used in fit
    cov_mat: covariance matrix

OPTIONS:
    BREAK_TIME: breakpoint time for piecewise regression
    BREAKPOINT: breakpoint indice of piecewise regression
    DATA_ERR: data precision
        single value if equal
        array if unequal for weighted least squares
    WEIGHT: Set if measurement errors for use in weighted least squares
    CYCLES: list of cyclical terms (0.5=semi-annual, 1=annual)
    STDEV: standard deviation of output error
    CONF: confidence interval of output error
    AICc: use second order AIC

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: added function docstrings
    Updated 10/2019: changing Y/N flags to True/False
    Updated 01/2019: added option S2 to include 161-day tidal aliasing terms
    Updated 12/2018: put transpose of design matrix within FIT_TYPE if statement
    Updated 08/2018: import packages before function definition
    Updated 09/2017: use rcond=-1 in numpy least-squares algorithms
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 11/2014: added simple output for model without climate oscillations
    Updated 07/2014: import scipy.stats and scipy.special
    Updated 06/2014: changed message to sys.exit
    Updated 02/2014: minor update to if statements
    Updated 10/2013: updated Log-likelihood (converted from Least-Squares (LS)
        log-likelihood to maximum likelihood (ML) log-likelihood)
        Added calculation for AICc (corrected for small sample size
    Updated 09/2013: updated weighted least-squares and added AIC, BIC and LOGLIK
        options for parameter evaluation
        Added cases with known standard deviations and weighted least-squares
        Minor change: P_cons, P_lin1 and P_lin2 changed to P_x0, P_x1a and P_x1b
    Updated 07/2013: updated errors for beta2
    Updated 05/2013: converted to python
    Updated 06/2012: added options for MSE and NRMSE.  adjusted rsq_adj
    Updated 06/2012: changed design matrix creation to 'FIT_TYPE'
        added r_square to calcuating goodness of fit compared to others
    Updated 03/2012: combined polynomial and harmonic regression functions
    Updated 01/2012: added std weighting for a error weighted least-squares
    Written 10/2011
"""
import warnings
import gravity_toolkit.time_series

def piecewise_regress(*args, **kwargs):
    """
    Fits a synthetic signal to data over a time period by ordinary or
    weighted least-squares for breakpoint analysis [Toms2003]_

    Parameters
    ----------
    t_in: float
        input time array
    d_in: float
        input data array
    BREAK_TIME: float or NoneType, default None
        breakpoint time for piecewise regression
    BREAKPOINT: int or NoneType, default None
        breakpoint indice of piecewise regression
    CYCLES: list, default [0.5, 1.0]
        list of cyclical terms in fractions of year
    DATA_ERR: float or list
        data precision

            - single value if equal
            - array if unequal for weighted least squares
    WEIGHT: bool, default False
        Use weighted least squares with measurement errors
    STDEV: float, default 0
        Standard deviation of output error
    CONF: float, default 0
        Confidence interval of output error
    AICc: bool, default False
        Use second order AIC for small sample sizes [Burnham2002]_

    Returns
    -------
    beta: float
        regressed coefficients array
    error: float
        regression fit error for each coefficient for an input deviation

            - ``STDEV``: standard deviation of output error
            - ``CONF``: confidence interval of output error
    std_err: float
        standard error for each coefficient
    R2: float
        coefficient of determination (r\ :sup:`2`)
    R2Adj: float
        r\ :sup:`2` adjusted for the number of terms in the model
    MSE: float
        mean square error
    WSSE: float
        Weighted sum of squares error
    NRMSE: float
        normalized root mean square error
    AIC: float
        Akaike information criterion
    BIC: float
        Bayesian information criterion (Schwarz criterion)
    model: float
        modeled timeseries
    simple: float
        modeled timeseries without oscillating components
    residual: float
        model residual
    DOF: int
        degrees of freedom
    N: int
        number of terms used in fit
    cov_mat: float
        covariance matrix

    References
    ----------
    .. [Toms2003] J. D. Toms and M. L. Lesperance,
        "Piecewise Regression: A Tool For Identifying Ecological
        Thresholds", *Ecology*, 84, 2034-2041, (2003).
        `doi: 10.1890/02-0472 <https://doi.org/10.1890/02-0472>`_
    .. [Burnham2002] K. P. Burnham and D. R. Anderson,
        *Model Selection and Multimodel Inference*,
        2nd Edition, 488 pp., (2002).
        `doi: 10.1007/b97636 <https://doi.org/10.1007/b97636>`_
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.time_series instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.time_series.piecewise(*args,**kwargs)
