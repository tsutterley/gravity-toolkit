degree_amplitude.py
====================

 - Calculates the amplitude of each spherical harmonic degree

#### Calling Sequence
```python
from gravity_toolkit.degree_amplitude import degree_amplitude
amp = degree_amplitude(clm,slm,LMAX=60)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/degree_amplitude.py)

#### Inputs:
 1. `clm`: cosine spherical harmonic coefficients
 2. `slm`: sine spherical harmonic coefficients

#### Options:
 - `LMAX`: Upper bound of Spherical Harmonic Degrees
 - `MMAX`: Upper bound of Spherical Harmonic Orders

#### Outputs:
 - `amp`: degree amplitude
