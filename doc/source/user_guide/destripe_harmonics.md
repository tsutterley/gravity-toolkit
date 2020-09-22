destripe_harmonics.py
=====================

 - Filters spherical harmonic coefficients for correlated "striping" errors following [Swenson and Wahr (2006)](http://dx.doi.org/10.1029/2005GL025285)  

#### Calling Sequence
```python
from gravity_toolkit.destripe_harmonics import destripe_harmonics
Ylms = destripe_harmonics(clm,slm,LMAX=60)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/destripe_harmonics.py)

#### Inputs
 1. `clm`: cosine spherical harmonic coefficients  
 2. `slm`: sine spherical harmonic coefficients  

#### Options
 - `LMIN`: Lower bound of Spherical Harmonic Degrees
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `MMAX`: Upper bound of Spherical Harmonic Orders
 - `ROUND`: use round to find nearest even (True) or use floor (False)
 - `NARROW`: Clm=Slm=0 if number of points is less than window size (False)

#### Outputs
 - `Wclm`: filtered cosine spherical harmonic coefficients
 - `Wslm`: filtered sine spherical harmonic coefficients
