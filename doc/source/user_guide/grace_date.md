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

#### Arguments
1. Working data directory for GRACE/GRACE-FO data

#### Keyword arguments
- `PROC`: Data processing center or satellite mission
   * `'CSR'`: University of Texas Center for Space Research
   * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
   * `'JPL'`: Jet Propulsion Laboratory
   * `'CNES'`: French Centre National D'Etudes Spatiales
   * `'GRAZ'`: Institute of Geodesy from GRAZ University of Technology
   * `'COSTG'`: Combination Service for Time-variable Gravity Fields
   * `'Swarm'`: Time-variable gravity data from Swarm satellites
- `DREL`: GRACE/GRACE-FO/Swarm data release
- `DSET`: GRACE/GRACE-FO/Swarm dataset
   * `'GAA'`: non-tidal atmospheric correction
   * `'GAB'`: non-tidal oceanic correction
   * `'GAC'`: combined non-tidal atmospheric and oceanic correction
   * `'GAD'`: ocean bottom pressure product
   * `'GSM'`: corrected monthly static gravity field product
- `OUTPUT`: create index file of dates for GRACE/GRACE-FO data
- `MODE`: permissions mode of output file

#### Returns
- dictionary of files mapped by GRACE/GRACE-FO month
