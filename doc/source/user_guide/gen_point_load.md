gen_point_load.py
=================

- Calculates gravitational spherical harmonic coefficients for point masses

#### Calling Sequence
```python
from gravity_toolkit.gen_point_load import gen_point_load
Ylms = gen_point_load(data, lon, lat, UNITS=1, LMAX=LMAX, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_point_load.py)

#### Inputs
- `data`: data magnitude
- `lon`: longitude of points
- `lat`: latitude of points

#### Options
- `UNITS`: input data units
   * `1`: grams of mass (default)
   * `2`: gigatonnes of mass (Gt)
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `LOVE`: input load Love numbers up to degree of truncation

#### Outputs
- `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
- `slm`: Sine spherical harmonic coefficients (geodesy normalization)
- `l`: spherical harmonic degree to LMAX
- `m`: spherical harmonic order to MMAX
