plm_columbo.py
==============

 - Computes fully-normalized associated Legendre Polynomials and their first derivative for a vector of x values using a standard forward column method
 - Uses the Colombo (1981) recursion relation listed in the [Geoid Cookbook](http://mitgcm.org/~mlosch/geoidcookbook.pdf) and [Holmes-Featherstone (2002)](https://doi.org/10.1007/s00190-002-0216-2)

#### Calling Sequence
```python
from gravity_toolkit.plm_colombo import plm_colombo
plm,dplm = plm_colombo(LMAX, x)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/gravity_toolkit/plm_columbo.py)

#### Inputs
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `x`: elements ranging from -1 to 1. Typically cos(theta), where theta is the colatitude in radians
        
#### Options
 - `ASTYPE`: output variable type. Default is 64-bit floating point

#### Outputs
 - `plms`: Legendre polynomials of x (geodesy normalization)
 - `dplms`: first differentials of Legendre polynomials of x
