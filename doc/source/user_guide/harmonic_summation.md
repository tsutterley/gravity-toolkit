harmonic_summation.py
====================

- Returns the spatial field for a series of spherical harmonics

#### Calling Sequence
```python
from gravity_toolkit.harmonic_summation import harmonic_summation
spatial = harmonic_summation(clm,slm,lon,lat,LMAX=60)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/harmonic_summation.py)

#### Inputs:
1. `clm`: cosine spherical harmonic coefficients
2. `slm`: sine spherical harmonic coefficients
3. `lon`: longitude
4. `lat`: latitude

#### Options:
- `LMIN`: Lower bound of Spherical Harmonic Degrees
- `LMAX`: Upper bound of Spherical Harmonic Degrees
- `MMAX`: Upper bound of Spherical Harmonic Orders
- `PLM`: Fully-normalized associated Legendre polynomials

#### Outputs:
- `spatial`: spatial field [lon,lat]

#### Dependencies
- `plm_holmes.py`: Computes fully-normalized associated Legendre polynomials
