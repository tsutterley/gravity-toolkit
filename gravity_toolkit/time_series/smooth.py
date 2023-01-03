#!/usr/bin/env python
u"""
smooth.py
Written by Tyler Sutterley (01/2023)

Computes a moving average of a time-series using three possible routines:
    1) centered moving average
    2) 13-month Loess filter (default)
    3) 13-month Loess filter weighted and outputs for all dates

Note: due to the missing months in the GRACE/GRACE-FO time series,
    a standard moving average will have problems if the
    missing months are not interpolated.

CALLING SEQUENCE:
    smth = gravity_toolkit.time_series.smooth(t_in, d_in, HFWTH=6)

INPUTS:
    t_in: input time array
    d_in: input data array

OUTPUTS:
    time: time after removing start and end half-windows
    data: smoothed time-series
    season: seasonal component calculated by the Loess filter
    annual: annual component calculated by the Loess filter
    semiann: semi-annual component calculated by the Loess filter
    trend: instantaneous trend calculated by the Loess filter
    error: estimated error of the instantaneous trend
    noise: noise component after removing the Loess trend and seasonal components
    reduce: original time series after removing start and end half-windows

OPTIONS:
    MOVING: calculates centered moving average using mean of window
        mean of: (January up to December) and (February up to January)
    WEIGHT: smoothing algorithm that backward models dates before
        half-width and forward models dates after half-width
        0: use unweighted Loess filter
        1: use linear weights with Loess filter
        2: use gaussian weights with Loess filter
    HFWTH: half-width of the moving average (default = 6 for 13-month Loess)
    DATA_ERR: input error for known and equal errors (single value)
    STDEV: standard deviation of output error
    CONF: confidence interval of output error (default is for 95%)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

UPDATE HISTORY:
    Updated 01/2023: refactored time series analysis functions
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 07/2020: added function docstrings
    Updated 10/2019: changing Y/N flags to True/False. output amplitudes
    Updated 09/2019: calculate and output annual and semi-annual phase
    Updated 08/2018: use implicit import of scipy stats and special
    Updated 09/2017: using rcond=-1 in numpy least-squares algorithms
    Updated 09/2014: added output for the reduced time-series
    Updated 06/2014: added parameter DATA_ERR for known and equal errors
    Updated 03/2014: created a new smoothing algorithm
        similar to Loess-type least-squares algorithm but has
        backward models dates before HFWTH and forward models dates after
        if all dates are available will be centrally weighted
        need to update to include error and should find a reference
            as i came up with this algorithm at 4:30am
    Updated 02/2014: minor update to if statements
    Updated 09/2013: switched MOVING flag to Y/N
        Minor change: P_cons, and P_lin changed to P_x0, and P_x1
    Updated 06/2013: added error for instantaneous trend
    Updated 04/2013: converted to python and added more outputs
    Updated 02/2013: using a centered moving average
        added seasonal option to compute the smooth seasonal variation
            calculated by the loess filter program.
        added option to compute the noise component after removing the
        smoothed trend and the seasonal component
    Updated 03/2012: added Loess smoothing following Velicogna (2009)
    Written 12/2011
"""
import numpy as np
import scipy.stats
import scipy.special

def smooth(t_in, d_in, HFWTH=6, MOVING=False, DATA_ERR=0, WEIGHT=0,
    STDEV=0, CONF=0):
    """
    Computes the moving average of a time-series

        1) centered moving average
        2) 13-month Loess filter [Velicogna2009]_
        3) 13-month Loess filter weighted and outputs for all dates

    Parameters
    ----------
    t_in: float
        input time array
    d_in: float
        input data array
    HFWTH: int
        half-width of the moving average
    MOVING: bool, default False
        calculates centered moving average using mean of window
    WEIGHT: smoothing algorithm that backward models dates before
        half-width and forward models dates after half-width

            - ``0``: use unweighted Loess filter
            - ``1``: use linear weights with Loess filter
            - ``2``: use gaussian weights with Loess filter
    DATA_ERR: float or list
        input error for known and equal errors
    STDEV: float, default 0
        Standard deviation of output error
    CONF: float, default 0
        Confidence interval of output error

    Returns
    -------
    time: float
        time after removing start and end half-windows
    data: float
        smoothed time-series
    season: float
        seasonal component calculated by the Loess filter
    annual: float
        annual component calculated by the Loess filter
    semiann: float
        semi-annual component calculated by the Loess filter
    trend: float
        instantaneous trend calculated by the Loess filter
    error: float
        estimated error of the instantaneous trend
    noise: float
        noise component after removing the Loess trend and seasonal components
    reduce: float
        original time series after removing start and end half-windows

    References
    ----------
    .. [Velicogna2009] I. Velicogna, "Increasing rates of ice mass loss
        from the Greenland and Antarctic ice sheets revealed by GRACE",
        *Geophysical Research Letters*, 36(L19503),
        `doi: 10.1029/2009GL040222 <https://doi.org/10.1029/2009GL040222>`_
    """

    # remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    nmax = len(t_in)

    # Indice with start of seasonal terms:
    SEAS = 2

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

    # moving average algorithm
    if MOVING:
        # Centered moving average using the mean of each window
        # equal to mean of Jan:Dec and Feb:Jan+1 for HFWTH 6
        # problematic with GRACE due to missing months within time-series
        # output time
        tout = t_in[HFWTH:nmax-HFWTH]
        smth = np.zeros((nmax-2*HFWTH))
        for k in range(0, (nmax-(2*HFWTH))):
            # centered moving average sum[2:i-1] + 0.5[1] + 0.5[i]
            smth[k] = np.sum(d_in[k+1:k+2*HFWTH]) + 0.5*(d_in[k]+d_in[k+2*HFWTH])
        dsmth = smth/(2*HFWTH)
        return {'data':dsmth, 'time':tout}
    elif WEIGHT in (1,2):
        # weighted moving average calculated from the least-squares of window
        # and removing An/SAn signal.  models entire range of dates
        # for a HFWTH of 6 (remove annual)
        # will fit linear model to data for 13 months
        # creates a weight array ranging from 1:HFWTH+1:-1 for linear
        # or a gaussian function centered on HFWTH
        # which favors the regression with the date centered
        # smoothed time-series = sum(smth*weights)/sum(weights)the weight array
        # output time = input time
        tout = np.copy(t_in)
        if (WEIGHT == 1):
            # linear weights (range from 1:HFWTH+1:-1)
            wi = np.concatenate((np.arange(1,HFWTH+2,dtype=np.float64),
                np.arange(HFWTH,0,-1,dtype=np.float64)),axis=0)
        elif (WEIGHT == 2):
            # gaussian weights
            # default standard deviation of 2
            stdev = 2.0
            # gaussian function over range 2*HFWTH
            # centered on HFWTH
            xi=np.arange(0, 2*HFWTH+1)
            wi=np.exp(-(xi-HFWTH)**2/(2.0*stdev**2))/(stdev*np.sqrt(2.0*np.pi))

        dsmth = np.zeros((nmax))
        dseason = np.zeros((nmax))
        dannual = np.zeros((nmax))
        annamp = np.zeros((nmax))
        annphase = np.zeros((nmax))
        dsemian = np.zeros((nmax))
        semiamp = np.zeros((nmax))
        semiphase = np.zeros((nmax))
        weight = np.zeros((nmax))
        for i in range(0, (nmax-(2*HFWTH))):
            ran = i + np.arange(0, 2*HFWTH+1)
            P_x0 = np.ones((2*HFWTH+1))# Constant Term
            P_x1 = t_in[ran]# Linear Term
            # Annual term = 2*pi*t*harmonic
            P_asin = np.sin(2*np.pi*t_in[ran])
            P_acos = np.cos(2*np.pi*t_in[ran])
            #Semi-Annual = 4*pi*t*harmonic
            P_ssin = np.sin(4*np.pi*t_in[ran])
            P_scos = np.cos(4*np.pi*t_in[ran])
            # x0,x1,AS,AC,SS,SC
            TMAT = np.array([P_x0, P_x1, P_asin, P_acos, P_ssin, P_scos])
            TMAT = np.transpose(TMAT)
            # Least-Squares fitting
            # (the [0] denotes coefficients output)standard
            beta_mat = np.linalg.lstsq(TMAT,d_in[ran],rcond=-1)[0]
            # Calculating the output components
            # add weighted smoothed time series
            dsmth[ran] += wi*np.dot(TMAT[:,0:SEAS],beta_mat[0:SEAS])
            # seasonal component
            dseason[ran] += wi*np.dot(TMAT[:,SEAS:],beta_mat[SEAS:])
            # annual component
            AS,AC = beta_mat[SEAS:SEAS+2]
            dannual[ran] += wi*np.dot(TMAT[:,SEAS:SEAS+2],[AS,AC])
            annamp[ran] += wi*np.sqrt(AS**2 + AC**2)
            annphase[ran] += wi*np.arctan2(AC,AS)*180.0/np.pi
            # semi-annual component
            SS,SC = beta_mat[SEAS+2:SEAS+4]
            dsemian[ran] += wi*np.dot(TMAT[:,SEAS+2:SEAS+4],[SS,SC])
            semiamp[ran] += wi*np.sqrt(SS**2 + SC**2)
            semiphase[ran] += wi*np.arctan2(SC,SS)*180.0/np.pi
            # add weights
            weight[ran] += wi
        # divide weighted smoothed time-series by weights
        # to get output smoothed time-series
        dsmth /= weight
        dseason /= weight
        dannual /= weight
        annamp /= weight
        annphase /= weight
        dsemian /= weight
        semiamp /= weight
        semiphase /= weight
        # noise = data - smoothed - seasonal
        dnoise = d_in - dsmth - dseason
        return {'data':dsmth, 'seasonal':dseason, 'annual':dannual,
            'annamp':annamp, 'annphase':annphase, 'semiann':dsemian,
            'semiamp':semiamp, 'semiphase':semiphase, 'noise':dnoise,
            'time':tout, 'weight':weight}
    else:
        # Moving average calculated from least-squares of window
        # and removing An/SAn signal
        # output time
        tout = t_in[HFWTH:nmax-HFWTH]
        dsmth = np.zeros((nmax-2*HFWTH))
        dtrend = np.zeros((nmax-2*HFWTH))
        derror = np.zeros((nmax-2*HFWTH))
        dseason = np.zeros((nmax-2*HFWTH))
        dannual = np.zeros((nmax-2*HFWTH))
        annamp = np.zeros((nmax-2*HFWTH))
        annphase = np.zeros((nmax-2*HFWTH))
        dsemian = np.zeros((nmax-2*HFWTH))
        semiamp = np.zeros((nmax-2*HFWTH))
        semiphase = np.zeros((nmax-2*HFWTH))
        dnoise = np.zeros((nmax-2*HFWTH))
        dreduce = np.zeros((nmax-2*HFWTH))
        for i in range(0, (nmax-(2*HFWTH))):
            ran = i + np.arange(0, 2*HFWTH+1)
            P_x0 = np.ones((2*HFWTH+1))# Constant Term
            P_x1 = t_in[ran]# Linear Term
            # Annual term = 2*pi*t*harmonic
            P_asin = np.sin(2*np.pi*t_in[ran])
            P_acos = np.cos(2*np.pi*t_in[ran])
            #Semi-Annual = 4*pi*t*harmonic
            P_ssin = np.sin(4*np.pi*t_in[ran])
            P_scos = np.cos(4*np.pi*t_in[ran])
            # x0,x1,AS,AC,SS,SC
            TMAT = np.array([P_x0, P_x1, P_asin, P_acos, P_ssin, P_scos])
            TMAT = np.transpose(TMAT)
            # Least-Squares fitting
            # (the [0] denotes coefficients output)
            beta_mat = np.linalg.lstsq(TMAT,d_in[ran],rcond=-1)[0]
            n_terms = len(beta_mat)

            if (DATA_ERR != 0):
                # LEAST-SQUARES CASE WITH KNOWN AND EQUAL ERROR
                P_err = DATA_ERR*np.ones((2*HFWTH+1))
                Hinv = np.linalg.inv(np.dot(np.transpose(TMAT),TMAT))
                #    Normal Equations
                NORMEQ = np.dot(Hinv,np.transpose(TMAT))
                beta_err = np.zeros((n_terms))
                for n in range(0,n_terms):
                    beta_err[n] = np.sqrt(np.sum((NORMEQ[n,:]*P_err)**2))
            else:
                # Error Analysis
                # Degrees of Freedom
                nu = (2*HFWTH+1) - n_terms
                # Mean square error
                MSE = np.dot(np.transpose(d_in[ran] - np.dot(TMAT,beta_mat)),
                    (d_in[ran] - np.dot(TMAT,beta_mat)))/nu
                # Covariance Matrix
                # Multiplying the design matrix by itself
                Hinv = np.linalg.inv(np.dot(np.transpose(TMAT),TMAT))
                # Taking the diagonal components of the cov matrix
                hdiag = np.diag(Hinv)

                # STANDARD LEAST-SQUARES CASE
                # Regression with Errors with Unknown Standard Deviations
                # Student T-Distribution with D.O.F. nu
                #    t.ppf parallels tinv in matlab
                tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)
                # beta_err is the error for each coefficient
                # beta_err = t(nu,1-alpha/2)*standard error
                st_err = np.sqrt(MSE*hdiag)
                beta_err = tstar*st_err

            # Calculating the output components
            # smoothed time series
            dsmth[i] = np.dot(TMAT[HFWTH,0:SEAS],beta_mat[0:SEAS])
            dtrend[i] = np.copy(beta_mat[1])# Instantaneous data trend
            derror[i] = np.copy(beta_err[1])# Error in trend
            # seasonal component
            dseason[i] = np.dot(TMAT[HFWTH,SEAS:],beta_mat[SEAS:])
            # annual component
            AS,AC = beta_mat[SEAS:SEAS+2]
            dannual[i] = np.dot(TMAT[HFWTH,SEAS:SEAS+2],[AS,AC])
            annphase[i] = np.arctan2(AC,AS)*180.0/np.pi
            annamp[i] = np.sqrt(AS**2 + AC**2)
            # semi-annual component
            SS,SC = beta_mat[SEAS+2:SEAS+4]
            dsemian[i] = np.dot(TMAT[HFWTH,SEAS+2:SEAS+4],[SS,SC])
            semiamp[i] = np.sqrt(SS**2 + SC**2)
            semiphase[i] = np.arctan2(SC,SS)*180.0/np.pi
            # noise component
            dnoise[i] = d_in[i+HFWTH] - dsmth[i] - dseason[i]
            # reduced time-series
            dreduce[i] = d_in[i+HFWTH]

        return {'data':dsmth, 'trend':dtrend, 'error':derror,
            'seasonal':dseason, 'annual':dannual, 'annphase':annphase,
            'annamp':annamp, 'semiann':dsemian, 'semiamp':semiamp,
            'semiphase':semiphase, 'noise':dnoise, 'time':tout, 'reduce':dreduce}
