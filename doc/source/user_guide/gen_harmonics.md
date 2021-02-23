gen_harmonics.py
================

- Converts data from the spatial domain to spherical harmonic coefficients
- Does not compute the solid Earth elastic response

#### Calling Sequence
```python
from gravity_toolkit.gen_harmonics import gen_harmonics
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_harmonics(data, lon, lat, LMAX=LMAX, PLM=PLM)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_harmonics.py)

#### Inputs
- `data`: data magnitude
- `lon`: longitude array
- `lat`: latitude array

#### Options
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `PLM`: input Legendre polynomials (for improving computational time)

#### Outputs
- `clm`: Cosine spherical harmonic coefficients (4-pi normalized)
- `slm`: Sine spherical harmonic coefficients (4-pi normalized)
- `l`: spherical harmonic degree to LMAX
- `m`: spherical harmonic order to MMAX
