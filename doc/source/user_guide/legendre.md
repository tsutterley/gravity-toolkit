legendre.py
===========

- Computes associated Legendre functions of degree l evaluated for elements x
- l must be a scalar integer and x must contain real values ranging -1 <= x <= 1

#### Calling Sequence
```python
from gravity_toolkit.legendre import legendre
Pl = legendre(l, x)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/legendre.py)

#### Arguments
- `l`: degree of Legrendre polynomials
- `x`: elements ranging from -1 to 1. Typically cos(theta), where theta is the colatitude in radians

#### Returns
- `Pl`: Legendre polynomials of degree l for orders 0 to l

#### Keyword arguments
- `NORMALIZE`: output Fully Normalized Associated Legendre Functions