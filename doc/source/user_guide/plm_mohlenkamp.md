plm_mohlenkamp.py
=================

 - Computes fully-normalized associated Legendre Polynomials for a vector of x values using Martin Mohlenkamp's recursion relation as listed in his [Guide to Spherical Harmonics](http://www.ohiouniversityfaculty.com/mohlenka/research/uguide.pdf)  
 - Derived from [Gabor Szeg&ouml; (1939)](https://people.math.osu.edu/nevai.1/AT/SZEGO/szego=szego1975=ops=OCR.pdf) recurrence formula for Jacobi Polynomials (Pg 71)


#### Calling Sequence
```python
from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp
plm = plm_mohlenkamp(LMAX, x)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/gravity_toolkit/plm_mohlenkamp.py)

#### Inputs
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `x`: elements ranging from -1 to 1. Typically cos(theta), where theta is the colatitude in radians

#### Options
 - `MMAX`: Upper bound of Spherical Harmonic Orders (default = LMAX)

#### Outputs
 - `plms`: Legendre polynomials of x (geodesy normalization)
