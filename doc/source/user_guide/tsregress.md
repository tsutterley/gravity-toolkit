tsregress.py
============

 - Fits a synthetic signal to the data over the time period by least-squares or weighted least-squares
 - Fit significance derivations are based on Burnham and Anderson (2002) Model Selection and Multimodel Inference

#### Calling Sequence
```python
from gravity_toolkit.tsregress import tsregress
tsbeta = tsregress(t_in, d_in, ORDER=1, CYCLES=[0.5,1.0], CONF=0.95)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/tsregress.py)

#### Inputs
 - `t_in`: input time array
 - `d_in`: input data array

#### Options
 - `DATA_ERR`: data precision (single value or array)
 - `WEIGHT`: use weighted least squares
 - `RELATIVE`: relative time period
 - `ORDER`: maximum polynomial order in fit
    0) constant
    1) linear
    2) quadratic
 - `CYCLES`: list of cyclical terms to include in fit
 - `STDEV`: standard deviation of output error
 - `CONF`: confidence interval of output error (default is for 95%)
 - `AICc`: use second order AIC for small sample sizes

#### Outputs
 - `beta`: regressed coefficients array
 - `error`: regression fit error for each coefficient for an input deviation
 - `std_err`: standard error for each coefficient
 - `R2`: coefficient of determination (r<sup>2</sup>)
 - `R2Adj`: coefficient of determination adjusted for the number of terms in the model
 - `MSE`: mean square error
 - `WSSE`: Weighted sum of squares error
 - `NRMSE`: normalized root mean square error
 - `AIC`: Akaike information criterion
 - `BIC`: Bayesian information criterion (Schwarz criterion)
 - `model`: modeled timeseries
 - `simple`: modeled timeseries without oscillating components
 - `residual`: model residual
 - `DOF`: degrees of freedom
 - `N`: number of terms used in fit
 - `cov_mat`: covariance matrix
