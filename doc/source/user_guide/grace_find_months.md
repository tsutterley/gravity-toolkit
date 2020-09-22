grace_find_months.py
====================

 - Parses date index file from grace_date.py
 - Finds the months available for a GRACE/GRACE-FO product
 - Finds the all months missing from the product

#### Calling Sequence
```python
from gravity_toolkit.grace_find_months import grace_find_months
grace_months = grace_find_months(base_dir, PROC, DREL, DSET=DSET)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/grace_find_months.py)

#### Inputs
 1. `base_dir`: Working data directory for GRACE/GRACE-FO data
 2. `PROC`: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)
    * `'CSR'`: University of Texas Center for Space Research
    * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
    * `'JPL'`: Jet Propulsion Laboratory
    * `'CNES'`: French Centre National D'Etudes Spatiales
 3. `DREL`: GRACE data release (RL04, RL05, RL06)

#### Options
 - `DSET`: GRACE dataset (GSM, GAC, GAD, GAB, GAA)
    * `'GAA'`: non-tidal atmospheric correction
    * `'GAB'`: non-tidal oceanic correction
    * `'GAC'`: combined non-tidal atmospheric and oceanic correction
    * `'GAD'`: GRACE/GRACE-FO ocean bottom pressure product
    * `'GSM'`: corrected monthly GRACE/GRACE-FO static field product

#### Outputs
 - `start`: First month in a GRACE/GRACE-FO dataset
 - `end`: Last month in a GRACE/GRACE-FO dataset
 - `missing`: missing months in a GRACE/GRACE-FO dataset
 - `months`: all available months in a GRACE/GRACE-FO dataset
 - `time`: center dates of all available months in a GRACE/GRACE-FO dataset
