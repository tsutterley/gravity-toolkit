geocenter.py
============

- Calculates the geocenter variation (in mm) from degree 1 spherical harmonic coefficients
- Calculates the Degree 1 Stokes Coefficients of a geocenter variation (in mm)

#### Calling Sequence
```python
from gravity_toolkit.geocenter import geocenter
xyz = geocenter(C10=C10, C11=C11, S11=S11)
Ylms = geocenter(X=x, Y=y, Z=z, INVERSE=True)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/geocenter.py)

#### Options
- `C10`: Cosine spherical harmonic of degree 1 and order 0
- `C11`: Cosine spherical harmonic of degree 1 and order 1
- `S11`: Sine spherical harmonic of degree 1 and order 1
- `X`: X-component of geocenter variation
- `Y`: Y-component of geocenter variation
- `Z`: Z-component of geocenter variation
- `RADIUS`:  Earth's radius for calculating spherical harmonics
- `INVERSE`: calculate the spherical harmonics coefficients from geocenter
