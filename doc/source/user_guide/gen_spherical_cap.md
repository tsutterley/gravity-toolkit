gen_spherical_cap.py
====================

- Calculates gravitational spherical harmonic coefficients for a spherical cap

#### Calling Sequence
```python
from gravity_toolkit.gen_spherical_cap import gen_spherical_cap
from gravity_toolkit.plm_holmes import plm_holmes
PLM,dPLM = plm_holmes(LMAX, np.cos(th))
Ylms = gen_spherical_cap(data, lon, lat, UNITS=1, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_spherical_cap.py)

#### Inputs
- `data`: data magnitude
- `lon`: longitude of spherical cap center
- `lat`: latitude of spherical cap center

#### Options
- `LMAX`:  maximum spherical harmonic degree of the output harmonics
- `MMAX`: maximum spherical harmonic order of the output harmonics
- `RAD_CAP`: spherical cap radius in degrees
- `RAD_KM`: spherical cap radius in kilometers
- `AREA`: spherical cap area (cm<sup>2</sup>)
- `UNITS`: input data units
   * `1`: cm water equivalent thickness (cm w.e., g/cm<sup>2</sup>)
   * `2` gigatonnes of mass (Gt)
   * `3` mm water equivalent thickness (mm w.e., kg/m<sup>2</sup>)
- `PLM`: input Legendre polynomials for spherical cap centers
- `LOVE`: input load Love numbers up to degree of truncation

#### Outputs
- `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
- `slm`: Sine spherical harmonic coefficients (geodesy normalization)
- `l`: spherical harmonic degree to LMAX
- `m`: spherical harmonic order to MMAX
