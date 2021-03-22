gen_stokes.py
=============

- Converts data from the spatial domain to spherical harmonic coefficients

#### Calling Sequence
```python
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_stokes(data, lon, lat, UNITS=1, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_stokes.py)

#### Arguments
- `data`: data matrix
- `lon`: longitude array
- `lat`: latitude array

#### Keyword arguments
- `UNITS`: input data units
   * `1` cm water equivalent thickness (cm w.e., g/cm<sup>2</sup>)
   * `2` gigatonnes of mass (Gt)
   * `3` mm water equivalent thickness (mm w.e., kg/m<sup>2</sup>)
- `LMIN`: minimum spherical harmonic degree of the output harmonics
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `PLM`: input Legendre polynomials (for improving computational time)
- `LOVE`: input load Love numbers up to degree of truncation

#### Returns
- `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
- `slm`: Sine spherical harmonic coefficients (geodesy normalization)
- `l`: spherical harmonic degree to LMAX
- `m`: spherical harmonic order to MMAX
