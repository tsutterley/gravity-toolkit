clenshaw_summation.py
=====================

- Returns the spatial field for a series of spherical harmonics at a sequence of ungridded points

#### Calling Sequence
```python
from gravity_toolkit.clenshaw_summation import clenshaw_summation
spatial = clenshaw_summation(clm,slm,lon,lat,UNITS=1,LMAX=60,LOVE=LOVE)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/clenshaw_summation.py)

#### Inputs:
1. `clm`: cosine spherical harmonic coefficients
2. `slm`: sine spherical harmonic coefficients
3. `lon`: longitude of points
4. `lat`: latitude of points

#### Options:
- `RAD`: Gaussian smoothing radius (km)
- `UNITS`: output data units
   * `1` cm of water thickness
   * `2` mm of geoid height
   * `3` mm of elastic crustal deformation ([Davis et al., 2004](https://doi.org/10.1029/2004GL021435))
   * `4` microGal gravitational perturbation
   * `5` Pa, equivalent surface pressure in Pascals
   * `6` cm of viscoelastic rustal uplift (GIA) ([Wahr et al., 2000](https://doi.org/10.1029/2000JB900113))
- `LMAX`: Upper bound of Spherical Harmonic Degrees
- `LOVE`: input load Love numbers up to degree of truncation (`hl`,`kl`,`ll`)
- `ASTYPE`: floating point precision for calculating Clenshaw summation
- `SCALE`: scaling factor to prevent underflow in Clenshaw summation

#### Outputs:
- `spatial`: spatial field

#### Dependencies
- `gauss_weights.py`: Computes the Gaussian weights as a function of degree
- `units.py`: Class for converting spherical harmonic data to specific units

#### References
- [Holmes and Featherstone, Journal of Geodesy (2002)](https://doi.org/10.1007/s00190-002-0216-2)
- Tscherning and Poder, Bollettino di Geodesia e Scienze (1982)
