read_tellus_geocenter.py
========================

- Reads monthly geocenter spherical harmonic data files from [GRACE Tellus Technical Notes (TN-13)](https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1) calculated following [Swenson et al. (2008)](https://doi.org/10.1029/2007JB005338)

#### Calling Sequence
```python
from gravity_toolkit.read_tellus_geocenter import read_tellus_geocenter
deg1_input = read_tellus_geocenter(geocenter_file, JPL=True)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_tellus_geocenter.py)

#### Inputs
- `geocenter_file`: degree 1 file
    * CSR: `TN-13_GEOC_CSR_RL06.txt`
    * GFZ: `TN-13_GEOC_GFZ_RL06.txt`
    * JPL: `TN-13_GEOC_JPL_RL06.txt`

#### Options
- `HEADER`: file contains header text to be skipped (default: True)
- `JPL`: use JPL TN-13 geocenter files calculated following [Sun et al., (2016)](https://doi.org/10.1007/s00190-015-0852-y)

#### Outputs
- `C10`: Cosine degree 1, order 0 spherical harmonic coefficients
- `C11`: Cosine degree 1, order 1 spherical harmonic coefficients
- `S11`: Sine degree 1, order 1 spherical harmonic coefficients
- `eC10`: Cosine degree 1, order 0 spherical harmonic coefficients Error
- `eC11`: Cosine degree 1, order 1 spherical harmonic coefficients Error
- `eS11`: Sine degree 1, order 1 spherical harmonic coefficients Error
- `month`: GRACE/GRACE-FO month (Apr 2002 = 004)
- `time`: date of GRACE/GRACE-FO month in decimal format
