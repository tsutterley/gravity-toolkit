grace_date.py
=============

- Reads GRACE/GRACE-FO index file from `podaac_grace_sync.py` or `gfz_isdc_grace_ftp.py`
- Parses dates of each GRACE/GRACE-FO file and assigns the month number
- Creates an index of dates for GRACE/GRACE-FO files

#### Calling Sequence
```python
from gravity_toolkit.grace_date import grace_date
grace_files = grace_date(base_dir, PROC=PROC, DREL=DREL, DSET=DSET)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/grace_date.py)

#### Inputs
1. Working data directory for GRACE/GRACE-FO data

#### Options
- `PROC`: GRACE/GRACE-FO data processing center (CSR/CNES/JPL/GFZ)
   * `'CSR'`: University of Texas Center for Space Research
   * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
   * `'JPL'`: Jet Propulsion Laboratory
   * `'CNES'`: French Centre National D'Etudes Spatiales
- `DREL`: GRACE/GRACE-FO data release (RL04/RL05/RL06)
- `DSET`: GRACE/GRACE-FO dataset (GAA/GAB/GAC/GAD/GSM)
   * `'GAA'`: non-tidal atmospheric correction
   * `'GAB'`: non-tidal oceanic correction
   * `'GAC'`: combined non-tidal atmospheric and oceanic correction
   * `'GAD'`: ocean bottom pressure product
   * `'GSM'`: corrected monthly static gravity field product
- `OUTPUT`: create index file of dates for GRACE/GRACE-FO data
- `MODE`: permissions mode of output file

#### Outputs
- dictionary of files mapped by GRACE/GRACE-FO month
