read_SLR_C20.py
===============

 - Reads monthly oblateness (degree 2 zonal) spherical harmonic data files from satellite laser ranging (SLR)

#### Calling Sequence
```python
from gravity_toolkit.read_SLR_C20 import read_SLR_C20
SLR_C20 = read_SLR_C20(SLR_file)
```

#### Inputs
 - `SLR_file`: oblateness file from satellite laser ranging
    - RL04: TN-05_C20_SLR.txt
    - RL05: TN-07_C20_SLR.txt
    - RL06: TN-11_C20_SLR.txt
    - CSR: C20_RL05.txt
    - GSFC: TN-14_C30_C30_GSFC_SLR.txt

#### Options
 - `HEADER`: file contains header text to be skipped (default: True)
 - `AOD`: remove background De-aliasing product from the SLR solution (for CSR)

#### Outputs
 - `data`: cosine degree 2 order 0 spherical harmonic coefficients (C20)
 - `error`: cosine degree 2 order 0 spherical harmonic coefficient errors (eC20)
 - `month`: GRACE/GRACE-FO month of measurement
 - `time`: date of SLR measurement
