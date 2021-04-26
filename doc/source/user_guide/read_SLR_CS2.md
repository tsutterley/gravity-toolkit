read_SLR_CS2.py
===============

- Reads monthly degree 2,m (figure axis and azimuthal dependence) spherical harmonic data files from satellite laser ranging (SLR)

#### Calling Sequence
```python
from gravity_toolkit.read_SLR_CS2 import read_SLR_CS2
SLR_CS2 = read_SLR_CS2(SLR_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_CS2.py)

#### Arguments
- `SLR_file`: oblateness file from satellite laser ranging
    * CSR 2,1: `C21_S21_RL06.txt`
    * CSR 2,2: `C22_S22_RL06.txt`

#### Returns
- `C2m`: cosine degree 2 order m spherical harmonic coefficients
- `S2m`: sine degree 2 order m spherical harmonic coefficients
- `eC2m`: cosine degree 2 order m cosine spherical harmonic coefficient errors
- `eS2m`: sine degree 2 order m sine spherical harmonic coefficient errors
- `month`: GRACE/GRACE-FO month of measurement
- `time`: date of SLR measurement
