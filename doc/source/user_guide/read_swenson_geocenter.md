read_swenson_geocenter.py
========================

- Reads [monthly geocenter coefficients](https://github.com/swensosc/GRACE_Tiles/blob/master/ancillary_data/gad_gsm.rl05.txt) from GRACE measurements and Ocean Models of Degree 1 provided by [Sean Swenson](https://doi.org/10.1029/2007JB005338) in mm w.e.


#### Calling Sequence
```python
from gravity_toolkit.read_swenson_geocenter import read_swenson_geocenter
deg1_input = read_swenson_geocenter(geocenter_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_swenson_geocenter.py)

#### Inputs
- `geocenter_file`: degree 1 file

#### Options
- `HEADER`: file contains header text to be skipped (default: True)

#### Outputs
- `C10`: Cosine degree 1, order 0 spherical harmonic coefficients
- `C11`: Cosine degree 1, order 1 spherical harmonic coefficients
- `S11`: Sine degree 1, order 1 spherical harmonic coefficients
- `month`: GRACE/GRACE-FO month (April 2002 = 004)
- `time`: date of GRACE/GRACE-FO month in decimal format
