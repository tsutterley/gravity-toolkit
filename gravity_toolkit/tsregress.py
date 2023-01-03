#!/usr/bin/env python
u"""
tsregress.py
Written by Tyler Sutterley (04/2022)

Fits a synthetic signal to data over a time period by ordinary or weighted
    least-squares

Fit significance derivations are based on Burnham and Anderson (2002)
    Model Selection and Multimodel Inference

CALLING SEQUENCE:
    tsbeta = tsregress(t_in, d_in, ORDER=1, CYCLES=[0.5,1.0], CONF=0.95)

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
    DATA_ERR: data precision
        single value if equal
        array if unequal for weighted least squares
    WEIGHT: Set if measurement errors for use in weighted least squares
    RELATIVE: relative period
    ORDER: maximum polynomial order in fit (0=constant, 1=linear, 2=quadratic)
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
    Updated 07/2020: added function docstrings
    Updated 10/2019: changing Y/N flags to True/False
    Updated 12/2018: put transpose of design matrix within FIT_TYPE if statement
    Updated 08/2018: import packages before function definition
    Updated 10/2017: output a seasonal model (will be 0 if no oscillating terms)
    Updated 09/2017: using rcond=-1 in numpy least-squares algorithms
    Updated 03/2017: added a catch for zero error in weighted least-squares
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 09/2014: made AICc option for second order AIC
        previously was default with no option for standard AIC
    Updated 07/2014: output the covariance matrix Hinv
        import scipy.stats and scipy.special
    Updated 06/2014: changed message to sys.exit
        new output for number of terms
    Updated 04/2014: added parameter RELATIVE for the relative time
    Updated 02/2014: minor update to if statements.  output simple regression
    Updated 10/2013:
        Added calculation for AICc (corrected for small sample size)
        Added output DOF (degrees of freedom, nu)
    Updated 09/2013: updated weighted least-squares and added AIC, BIC and
        LOGLIK options for parameter evaluation
        Fixed case with known standard deviations
        Changed Weight flag to Y/N
        Minor change: P_cons, P_lin and P_quad changed to P_x0, P_x1 and P_x2
    Updated 07/2013: added output for the modelled time-series
    Updated 05/2013: converted to Python
    Updated 02/2013: added in case for equal data error
        and made associated edits.. added option WEIGHT
    Updated 10/2012: added in additional FIT_TYPES that do not have a trend
        added option for the Schwarz Criterion for model selection
    Updated 08/2012: added option for changing the confidence interval or
        adding a standard deviation of the error
        changed 'std' to data_err for measurement errors
    Updated 06/2012: added options for MSE and NRMSE.  adjusted rsq_adj
    Updated 06/2012: changed design matrix creation to 'FIT_TYPE'
        added r_square to calcuating goodness of fit compared to others
    Updated 03/2012: combined polynomial and harmonic regression functions
    Updated 01/2012: added std weighting for a error weighted least-squares
    Written 10/2011
"""
import warnings
import gravity_toolkit.time_series

def tsregress(*args, **kwargs):
    """
    Fits a synthetic signal to data over a time period by
        ordinary or weighted least-squares

    Parameters
    ----------
    t_in: float
        input time array
    d_in: float
        input data array
    ORDER: int, default 1
        maximum polynomial order in fit

            * ``0``: constant
            * ``1``: linear
            * ``2``: quadratic
    CYCLES: list, default [0.5, 1.0]
        list of cyclical terms
    DATA_ERR: float or list
        Data precision

            - single value if equal
            - array if unequal for weighted least squares
    WEIGHT: bool, default False
        Use weighted least squares with measurement errors
    RELATIVE: float or List, default Ellipsis
        Epoch for calculating relative dates

            - float: use exact value as epoch
            - list: use mean from indices of available times
            - ``Ellipsis``: use mean of all available times
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

    Reference
    ---------
    .. [Burnham2002] K. P. Burnham and D. R. Anderson,
        *Model Selection and Multimodel Inference*,
        2nd Edition, 488 pp., (2002).
        `doi: 10.1007/b97636 <https://doi.org/10.1007/b97636>`_
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use gravity_toolkit.time_series instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return gravity_toolkit.time_series.regress(*args, **kwargs)
