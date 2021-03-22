read_SLR_C30.py
===============

- Reads monthly degree 3 zonal spherical harmonic data files from satellite laser ranging (SLR)

#### Calling Sequence
```python
from gravity_toolkit.read_SLR_C30 import read_SLR_C30
SLR_C30 = read_SLR_C30(SLR_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C30.py)

#### Arguments
- `SLR_file`: low degree zonal file from satellite laser ranging
    * CSR: `CSR_Monthly_5x5_Gravity_Harmonics.txt`
    * GSFC: `TN-14_C30_C30_GSFC_SLR.txt`
    * LARES: `C30_LARES_filtered.txt`

#### Keyword arguments
- `HEADER`: file contains header text to be skipped (default: True)
- `C30_MEAN`: mean C30 to add to LARES C30 anomalies

#### Returns
- `data`: cosine degree 3 order 0 spherical harmonic coefficients (C30)
- `error`: cosine degree 3 order 0 spherical harmonic coefficient errors (eC30)
- `month`: GRACE/GRACE-FO month of measurement
- `time`: date of SLR measurement
