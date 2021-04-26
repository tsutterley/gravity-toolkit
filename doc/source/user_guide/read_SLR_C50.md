read_SLR_C50.py
===============

- Reads monthly degree 5 zonal spherical harmonic data files from satellite laser ranging (SLR)

#### Calling Sequence
```python
from gravity_toolkit.read_SLR_C50 import read_SLR_C50
SLR_C50 = read_SLR_C50(SLR_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_C50.py)

#### Arguments
- `SLR_file`: low degree zonal file from satellite laser ranging
    * CSR: `CSR_Monthly_5x5_Gravity_Harmonics.txt`
    * GSFC: `GSFC_SLR_C20_C30_C50_GSM_replacement.txt`
    * LARES: `C50_LARES_filtered.txt`

#### Keyword arguments
- `HEADER`: file contains header text to be skipped (default: True)
- `C50_MEAN`: mean C50 to add to LARES C50 anomalies

#### Returns
- `data`: cosine degree 5 order 0 spherical harmonic coefficients (C50)
- `error`: cosine degree 5 order 0 spherical harmonic coefficient errors (eC50)
- `month`: GRACE/GRACE-FO month of measurement
- `time`: date of SLR measurement
