gen_pressure_stokes.py
======================

 - Converts pressure fields from the spatial domain to spherical harmonic coefficients

#### Calling Sequence
```python
from gravity_toolkit.gen_pressure_stokes import gen_pressure_stokes
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_stokes(PG, R, lon, lat, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_pressure_stokes.py)

#### Inputs
 - `PG`: pressure/gravity ratio
 - `R`: Earth's radius at each data point
 - `lon`: longitude array
 - `lat`: latitude array

#### Options
 - `LMAX`:  maximum spherical harmonic degree of the output harmonics
 - `MMAX`: maximum spherical harmonic order of the output harmonics
 - `PLM`: input Legendre polynomials (for improving computational time)
 - `LOVE`: input load Love numbers up to degree `LMAX` (hl,kl,ll)

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
 - `slm`: Sine spherical harmonic coefficients (geodesy normalization)
 - `l`: spherical harmonic degree to LMAX
 - `m`: spherical harmonic order to MMAX
