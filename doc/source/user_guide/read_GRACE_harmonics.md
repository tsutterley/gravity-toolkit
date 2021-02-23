read_GRACE_harmonics.py
=======================

- Reads GRACE/GRACE-FO files and extracts spherical harmonic data and drift rates (RL04)
- Adds drift rates to clm and slm for release 4 harmonics
- Correct GSM data for drift in pole tide following Wahr et al. (2015)
- Extracts start and end date of GRACE/GRACE-FO files and calculates mean of range

#### Calling Sequence
```python
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
CSR_L2_input = read_GRACE_harmonics('GSM-2_2002095-2002120_0021_UTCSR_0060_0005.gz',60)
GFZ_L2_input = read_GRACE_harmonics('GSM-2_2002094-2002120_0024_EIGEN_G---_005a.gz',90)
JPL_L2_input = read_GRACE_harmonics('GSM-2_2002091-2002120_0018_JPLEM_0001_0005.gz',60)
JPLMSC_input = read_GRACE_harmonics('GSM-2_2003001-2003031_0029_JPLMSC_0719_0005',719)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_GRACE_harmonics.py)

#### Inputs
1. full path to input GRACE file
2. spherical harmonic degree of truncation (`LMAX`)

#### Options
- `MMAX`: spherical harmonic order of truncation (default is `LMAX`)
- `POLE_TIDE`: correct GSM data for pole tide drift following [Wahr et al. (2015)](https://doi.org/10.1002/2015JB011986)

#### Outputs
- `time`: mid-month date of GRACE file in year-decimal
- `start`: start date of range as Julian day
- `end`: end date of range as Julian day
- `clm`: cosine spherical harmonics of input data
- `slm`: sine spherical harmonics of input data
- `eclm`: cosine spherical harmonic uncalibrated standard deviations
- `eslm`: sine spherical harmonic uncalibrated standard deviations
- `header`: text header of the GRACE file (will parse new YAML headers)
