legendre_polynomials.py
=======================

- Computes fully normalized Legendre polynomials for an array of x values and their first derivative following [Hofmann-Wellenhof and Moritz, 2006](http://www.springerlink.com/content/978-3-211-33544-4)
- Calculates Legendre polynomials for zonal harmonics (order 0)

#### Calling Sequence
```python
from gravity_toolkit.legendre_polynomials import legendre_polynomials
pl,dpl = legendre_polynomials(LMAX, x)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/legendre_polynomials.py)

#### Inputs
- `LMAX`: Upper bound of Spherical Harmonic Degrees
- `x`: elements ranging from -1 to 1. Typically cos(theta), where theta is the colatitude in radians

#### Options
- `ASTYPE`: output variable type. Default is 64-bit floating point

#### Outputs
- `pl`: Legendre polynomials of x (geodesy normalization)
- `dpl`: first differentials of Legendre polynomials of x
