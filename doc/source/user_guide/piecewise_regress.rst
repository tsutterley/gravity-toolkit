====================
piecewise_regress.py
====================

- Fits a synthetic signal to data over a time period by ordinary or weighted least-squares for breakpoint analysis [Toms2003]_
- Fit significance derivations are based on [Burnham2002]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.piecewise_regress import piecewise_regress
    tsbeta = piecewise_regress(t_in, d_in, BREAKPOINT=len(t_in)//2, CYCLES=[0.5,1.0])

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/piecewise_regress.py

Arguments
#########

- ``t_in``: input time array
- ``d_in``: input data array

Keyword arguments
#################

- ``BREAK_TIME``: breakpoint time for piecewise regression
- ``BREAKPOINT``: breakpoint indice of piecewise regression
- ``DATA_ERR``: data precision (single value or array)
- ``WEIGHT``: use weighted least squares
- ``CYCLES``: list of cyclical terms to include in fit
- ``STDEV``: standard deviation of output error
- ``CONF``: confidence interval of output error (default is for 95%)
- ``AICc``: use second order AIC for small sample sizes

Returns
#######

- ``beta``: regressed coefficients array
- ``error``: regression fit error for each coefficient for an input deviation
- ``std_err``: standard error for each coefficient
- ``R2``: coefficient of determination (r\ :sup:`2`)
- ``R2Adj``: coefficient of determination adjusted for the number of terms in the model
- ``MSE``: mean square error
- ``WSSE``: Weighted sum of squares error
- ``NRMSE``: normalized root mean square error
- ``AIC``: Akaike information criterion
- ``BIC``: Bayesian information criterion (Schwarz criterion)
- ``model``: modeled timeseries
- ``simple``: modeled timeseries without oscillating components
- ``residual``: model residual
- ``DOF``: degrees of freedom
- ``N``: number of terms used in fit
- ``cov_mat``: covariance matrix

References
##########

.. [Toms2003] J. D. Toms and M. L. Lesperance, "Piecewise Regression: A Tool For Identifying Ecological Thresholds", *Ecology*, 84, 2034-2041, (2003). `doi: 10.1890/02-0472 <https://doi.org/10.1890/02-0472>`_

.. [Burnham2002] K. P. Burnham and D. R. Anderson, *Model Selection and Multimodel Inference*, 2nd Edition, 488 pp., (2002). `doi: 10.1007/b97636 <https://doi.org/10.1007/b97636>`_
