read_ICGEM_harmonics.py
=======================

 - Reads gravity model files and extracts spherical harmonic data from the [GFZ International Centre for Global Earth Models (ICGEM)](http://icgem.gfz-potsdam.de/)

#### Calling Sequence
```python
from gravity_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
Ylms = read_ICGEM_harmonics(model_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/gravity_toolkit/read_ICGEM_harmonics.py)

#### Inputs
 1. `model_file`: full path to GFZ ICGEM gfc spherical harmonic data file

#### Options
 - `FLAG`: string denoting data lines  

#### Outputs
 - `clm`: cosine spherical harmonics of input data
 - `slm`: sine spherical harmonics of input data
 - `eclm`: cosine spherical harmonic standard deviations of type errors
 - `eslm`: sine spherical harmonic standard deviations of type errors
 - `modelname`: name of the gravity model
 - `earth_gravity_constant`: GM constant of the Earth for the gravity model
 - `radius`: semi-major axis of the Earth for the gravity model
 - `max_degree`: maximum degree and order for the gravity model
 - `errors`: error type of the gravity model
 - `norm`: normalization of the spherical harmonics
 - `tide_system`: tide system of gravity model (mean_tide, zero_tide, tide_free)
