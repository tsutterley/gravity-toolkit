#!/usr/bin/env python
u"""
piecewise.py
Written by Tyler Sutterley (01/2023)

Fits a synthetic signal to data over a time period by ordinary or weighted
    least-squares for breakpoint analysis

Derivation of Sharp Breakpoint Piecewise Regression:
    https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/02-0472
        y = beta_0 + beta_1*t + e (for x <= alpha)
        y = beta_0 + beta_1*t + beta_2*(t-alpha) + e (for x > alpha)

Fit significance derivations are based on Burnham and Anderson (2002)
    Model Selection and Multimodel Inference

CALLING SEQUENCE:
    pcwbeta = gravity_toolkit.time_series.piecewise(tdec, data,
        CYCLES=[0.5,1.0], BREAKPOINT=ind)

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
    Updated 01/2023: refactored time series analysis functions
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
import numpy as np
import scipy.stats
import scipy.special

def piecewise(t_in, d_in, BREAK_TIME=None, BREAKPOINT=None,
    CYCLES=[0.5,1.0], DATA_ERR=0, WEIGHT=False, STDEV=0, CONF=0,
    AICc=False):
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

    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    nmax = len(t_in)

    # If indice of cutoff time entered: will calculate cutoff time
    # If cutoff time entered: will find the cutoff indice
    if BREAKPOINT is not None:
        tco = t_in[BREAKPOINT]
        nco = np.squeeze(BREAKPOINT)
    elif BREAK_TIME is not None:
        nco = np.argmin(np.abs(t_in - BREAK_TIME))
        tco = np.copy(BREAK_TIME)

    # create design matrix for sharp breakpoint piecewise regression
    # y = beta_0 + beta_1*t + e (for x <= alpha)
    # y = beta_0 + beta_1*t + beta_2*(t-alpha) + e (for x > alpha)
    DMAT = []
    # add polynomial orders (0=constant, 1=linear)
    for o in range(2):
        DMAT.append(t_in**o)
    # Linear Term 2 (change from linear term1: trend2 = beta1+beta2)
    P_x1 = np.zeros((nmax))
    P_x1[nco:] = t_in[nco:] - tco
    DMAT.append(P_x1)
    # add cyclical terms (0.5=semi-annual, 1=annual)
    for c in CYCLES:
        DMAT.append(np.sin(2.0*np.pi*t_in/np.float64(c)))
        DMAT.append(np.cos(2.0*np.pi*t_in/np.float64(c)))
    # take the transpose of the design matrix
    DMAT = np.transpose(DMAT)

    # Calculating Least-Squares Coefficients
    if WEIGHT:
        # Weighted Least-Squares fitting
        if (np.ndim(DATA_ERR) == 0):
            raise ValueError('Input DATA_ERR for Weighted Least-Squares')
        # check if any error values are 0 (prevent infinite weights)
        if np.count_nonzero(DATA_ERR == 0.0):
            # change to minimum floating point value
            DATA_ERR[DATA_ERR == 0.0] = np.finfo(np.float64).eps
        # Weight Precision
        wi = np.squeeze(DATA_ERR**(-2))
        # If uncorrelated weights are the diagonal
        W = np.diag(wi)
        # Least-Squares fitting
        # Temporary Matrix: Inv(X'.W.X)
        TM1 = np.linalg.inv(np.dot(np.transpose(DMAT),np.dot(W,DMAT)))
        # Temporary Matrix: (X'.W.Y)
        TM2 = np.dot(np.transpose(DMAT),np.dot(W,d_in))
        # Least Squares Solutions: Inv(X'.W.X).(X'.W.Y)
        beta_mat = np.dot(TM1,TM2)
    else:# Standard Least-Squares fitting (the [0] denotes coefficients output)
        beta_mat = np.linalg.lstsq(DMAT,d_in,rcond=-1)[0]
        # Weights are equal
        wi = 1.0

    # Calculating trend2 = beta1 + beta2
    # beta2 = change in linear term from beta1
    beta_out = np.copy(beta_mat)# output beta
    beta_out[2] = beta_mat[1] + beta_mat[2]

    # number of terms in least-squares solution
    n_terms = len(beta_mat)
    # modelled time-series
    mod = np.dot(DMAT,beta_mat)
    # time-series residuals
    res = d_in[0:nmax] - np.dot(DMAT,beta_mat)
    # Fitted Values without climate oscillations
    simple = np.dot(DMAT[:,0:3],beta_mat[0:3])

    # Error Analysis
    # nu = Degrees of Freedom = number of measurements-number of parameters
    nu = nmax - n_terms

    # calculating R^2 values
    # SStotal = sum((Y-mean(Y))**2)
    SStotal = np.dot(np.transpose(d_in[0:nmax] - np.mean(d_in[0:nmax])),
        (d_in[0:nmax] - np.mean(d_in[0:nmax])))
    # SSerror = sum((Y-X*B)**2)
    SSerror = np.dot(np.transpose(d_in[0:nmax] - np.dot(DMAT,beta_mat)),
        (d_in[0:nmax] - np.dot(DMAT,beta_mat)))
    # R**2 term = 1- SSerror/SStotal
    rsquare = 1.0 - (SSerror/SStotal)
    # Adjusted R**2 term: weighted by degrees of freedom
    rsq_adj = 1.0 - (SSerror/SStotal)*np.float64((nmax-1.0)/nu)
    # Fit Criterion
    # number of parameters including the intercept and the variance
    K = np.float64(n_terms + 1)
    # Log-Likelihood with weights (if unweighted, weight portions == 0)
    # log(L) = -0.5*n*log(sigma^2) - 0.5*n*log(2*pi) - 0.5*n
    #log_lik = -0.5*nmax*(np.log(2.0 * np.pi) + 1.0 + np.log(np.sum((res**2)/nmax)))
    log_lik = 0.5*(np.sum(np.log(wi)) - nmax*(np.log(2.0 * np.pi) + 1.0 -
        np.log(nmax) + np.log(np.sum(wi * (res**2)))))

    # Aikaike's Information Criterion
    AIC = -2.0*log_lik + 2.0*K
    if AICc:
        # Second-Order AIC correcting for small sample sizes (restricted)
        # Burnham and Anderson (2002) advocate use of AICc where
        # ratio num/K is small
        # A small ratio is defined in the definition at approximately < 40
        AIC += (2.0*K*(K+1.0))/(nmax - K - 1.0)

    # Bayesian Information Criterion (Schwarz Criterion)
    BIC = -2.0*log_lik + np.log(nmax)*K

    # Error Analysis
    if WEIGHT:
        # WEIGHTED LEAST-SQUARES CASE (unequal error)
        # Covariance Matrix
        Hinv = np.linalg.inv(np.dot(np.transpose(DMAT),np.dot(W,DMAT)))
        # Normal Equations
        NORMEQ = np.dot(Hinv,np.transpose(np.dot(W,DMAT)))
        temp_err = np.zeros((n_terms))
        # Propagating RMS errors
        for i in range(0,n_terms):
            temp_err[i] = np.sqrt(np.sum((NORMEQ[i,:]*DATA_ERR)**2))

        # Recalculating beta2 error
        beta_err = np.copy(temp_err)
        beta_err[2] = np.sqrt(temp_err[1]**2 + temp_err[2]**2)
        # Weighted sum of squares Error
        WSSE = np.dot(np.transpose(wi*(d_in[0:nmax] - np.dot(DMAT,beta_mat))),
            wi*(d_in[0:nmax] - np.dot(DMAT,beta_mat)))/np.float64(nu)

        return {'beta':beta_out, 'error':beta_err, 'R2':rsquare,
            'R2Adj':rsq_adj, 'WSSE':WSSE, 'AIC':AIC, 'BIC':BIC,
            'LOGLIK':log_lik, 'model':mod, 'residual':res,
            'N':n_terms, 'DOF':nu, 'cov_mat':Hinv}

    elif ((not WEIGHT) and (DATA_ERR != 0)):
        # LEAST-SQUARES CASE WITH KNOWN AND EQUAL ERROR
        P_err = DATA_ERR*np.ones((nmax))
        Hinv = np.linalg.inv(np.dot(np.transpose(DMAT),DMAT))
        # Normal Equations
        NORMEQ = np.dot(Hinv,np.transpose(DMAT))
        temp_err = np.zeros((n_terms))
        for i in range(0,n_terms):
            temp_err[i] = np.sum((NORMEQ[i,:]*P_err)**2)
        # Recalculating beta2 error
        beta_err = np.copy(temp_err)
        beta_err[2] = np.sqrt(temp_err[1]**2 + temp_err[2]**2)
        # Mean square error
        MSE = np.dot(np.transpose(d_in[0:nmax] - np.dot(DMAT,beta_mat)),
            (d_in[0:nmax] - np.dot(DMAT,beta_mat)))/np.float64(nu)

        return {'beta':beta_out, 'error':beta_err, 'R2':rsquare,
            'R2Adj':rsq_adj, 'MSE':MSE, 'AIC':AIC, 'BIC':BIC,
            'LOGLIK':log_lik, 'model':mod, 'residual':res,
            'N':n_terms, 'DOF':nu, 'cov_mat':Hinv}
    else:
        # STANDARD LEAST-SQUARES CASE
        # Regression with Errors with Unknown Standard Deviations
        # MSE = (1/nu)*sum((Y-X*B)**2)
        # Mean square error
        MSE = np.dot(np.transpose(d_in[0:nmax] - np.dot(DMAT,beta_mat)),
            (d_in[0:nmax] - np.dot(DMAT,beta_mat)))/np.float64(nu)
        # Root mean square error
        RMSE = np.sqrt(MSE)
        # Normalized root mean square error
        NRMSE = RMSE/(np.max(d_in[0:nmax])-np.min(d_in[0:nmax]))
        # Covariance Matrix
        # Multiplying the design matrix by itself
        Hinv = np.linalg.inv(np.dot(np.transpose(DMAT),DMAT))
        # Taking the diagonal components of the cov matrix
        hdiag = np.diag(Hinv)
        # set either the standard deviation or the confidence interval
        if (STDEV != 0):
            # Setting the standard deviation of the output error
            alpha = 1.0 - scipy.special.erf(STDEV/np.sqrt(2.0))
        elif (CONF != 0):
            # Setting the confidence interval of the output error
            alpha = 1.0 - CONF
        else:
            # Default is 95% confidence interval
            alpha = 1.0 - (0.95)
        # Student T-Distribution with D.O.F. nu
        # t.ppf parallels tinv in matlab
        tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)
        # beta_err is the error for each coefficient
        # beta_err = t(nu,1-alpha/2)*standard error
        temp_std = np.sqrt(MSE*hdiag)
        temp_err = tstar*temp_std

        # Recalculating standard error for beta2
        st_err = np.copy(temp_std)
        st_err[2] = np.sqrt(temp_std[1]**2 + temp_std[2]**2)
        # Recalculating beta2 error
        beta_err = np.copy(temp_err)
        beta_err[2] = np.sqrt(temp_err[1]**2 + temp_err[2]**2)

        return {'beta':beta_out, 'error':beta_err, 'std_err':st_err, 'R2':rsquare,
            'R2Adj':rsq_adj, 'MSE':MSE, 'NRMSE':NRMSE, 'AIC':AIC, 'BIC':BIC,
            'LOGLIK':log_lik, 'model':mod, 'simple': simple, 'residual':res,
            'N':n_terms, 'DOF': nu, 'cov_mat':Hinv}
