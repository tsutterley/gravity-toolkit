plm_columbo.py
==============

 - Computes fully-normalized associated Legendre Polynomials for a vector of x values using a standard forward column method
 - Uses the Colombo (1981) recursion relation listed in the [Geoid Cookbook](http://mitgcm.org/~mlosch/geoidcookbook.pdf) and [Holmes-Featherstone (2002)](https://doi.org/10.1007/s00190-002-0216-2)

#### Calling Sequence
```python
from gravity_toolkit.plm_colombo import plm_colombo
plm,dplm = plm_colombo(LMAX, x)
```

#### Inputs
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `x`: typically cos(theta), where theta is the colatitude in radians

#### Options
 - `ASTYPE`: output variable type. Default is 64-bit floating point

#### Outputs
 - `plms`: Legendre polynomials of x (geodesy normalization)
 - `dplms`: first differentials of Legendre polynomials of x
