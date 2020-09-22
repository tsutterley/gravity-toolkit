read_SLR_geocenter.py
=====================

 - Reads monthly geocenter spherical harmonic data files from [satellite laser ranging (SLR)](ftp://ftp.csr.utexas.edu/pub/slr/geocenter/)

#### Calling Sequence
```python
from gravity_toolkit.read_SLR_geocenter import read_SLR_geocenter
deg1_input = read_SLR_geocenter(geocenter_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_SLR_geocenter.py)

#### Inputs
 - `geocenter_file`: degree 1 file
     - RL04: GCN_RL04.txt
     - RL05: GCN_RL05.txt
     - RL06: GCN_RL06.txt
     - CF-CM: GCN_L1_L2_30d_CF-CM.txt

#### Options
 - `RADIUS`: Earth's radius for calculating spherical harmonics from SLR data
 - `skiprows`: Rows of data to skip when importing data

#### Outputs
 - `C10`: Cosine degree 1, order 0 spherical harmonic coefficients
 - `C11`: Cosine degree 1, order 1 spherical harmonic coefficients
 - `S11`: Sine degree 1, order 1 spherical harmonic coefficients
 - `eC10`: Cosine degree 1, order 0 spherical harmonic coefficients Error
 - `eC11`: Cosine degree 1, order 1 spherical harmonic coefficients Error
 - `eS11`: Sine degree 1, order 1 spherical harmonic coefficients Error
 - `month`: GRACE/GRACE-FO month (Apr 2002 = 004)
 - `time`: date of GRACE/GRACE-FO month in decimal format
