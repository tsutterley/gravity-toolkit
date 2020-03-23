plm_holmes.py
=============

 - Computes fully-normalized associated Legendre Polynomials for a vector of x values using the [Holmes and Featherstone (2002)](https://doi.org/10.1007/s00190-002-0216-2) recursion relation  
 - Recursion relation is stable up to very high degree and order  

#### Calling Sequence
```
from gravity_toolkit.plm_holmes import plm_holmes
plm,dplm = plm_holmes(LMAX, x)
```

#### Inputs
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `x`: typically cos(theta), where theta is the colatitude in radians

#### Options
 - `ASTYPE`: output variable type. Default is 64-bit floating point

#### Outputs
 - `plms`: Legendre polynomials of x (geodesy normalization)
 - `dplms`: first differentials of Legendre polynomials of x
