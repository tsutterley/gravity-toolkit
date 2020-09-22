read_CSR_monthly_6x1.py
=======================

 - Reads in monthly 5x5 spherical harmonic coefficients with 1 coefficient from degree 6 all calculated from satellite laser ranging (SLR) measurements calculated by the [University of Texas Center for Space Research (CSR)](https://doi.org/10.1029/2010JB000850)

#### Calling Sequence
```python
from gravity_toolkit.read_CSR_monthly_6x1 import read_CSR_monthly_6x1
Ylms = read_CSR_monthly_6x1(input_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_CSR_monthly_6x1.py)

#### Inputs
 - `input_file`: input satellite laser ranging file

#### Options
 - `HEADER`: file contains header text to be skipped (default: True)

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients
 - `slm`: Sine spherical harmonic coefficients
 - `error/clm`: Cosine spherical harmonic coefficient uncertainty
 - `error/slm`: Sine spherical harmonic coefficients uncertainty
 - `MJD`: output date as Modified Julian Day
 - `date`: output date in year-decimal
