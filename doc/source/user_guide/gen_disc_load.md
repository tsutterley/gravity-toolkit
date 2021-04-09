gen_disc_load.py
================

- Calculates gravitational spherical harmonic coefficients for a uniform disc load

#### Calling Sequence
```python
from gravity_toolkit.gen_disc_load import gen_disc_load
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_disc_load(data, lon, lat, area, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_disc_load.py)

#### Arguments
- `data`: data magnitude (Gt)
- `lon`: longitude of disc center
- `lat`: latitude of disc center
- `area`: area of disc (km<sup>2</sup>)

#### Keyword arguments
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `PLM`: input Legendre polynomials for `cos(theta)` (disc center)
- `LOVE`: input load Love numbers up to degree of truncation

#### Returns
- `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
- `slm`: Sine spherical harmonic coefficients (geodesy normalization)
- `l`: spherical harmonic degree to `LMAX`
- `m`: spherical harmonic order to `MMAX`
