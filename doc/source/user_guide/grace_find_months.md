grace_find_months.py
====================

- Parses date index file from `grace_date` program
- Finds the months available for a GRACE/GRACE-FO/Swarm product
- Finds the all months missing from the product

#### Calling Sequence
```python
from gravity_toolkit.grace_find_months import grace_find_months
grace_months = grace_find_months(base_dir, PROC, DREL, DSET=DSET)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/grace_find_months.py)

#### Arguments
1. `base_dir`: Working data directory for GRACE/GRACE-FO data
2. `PROC`: Data processing center or satellite mission
    * `'CSR'`: University of Texas Center for Space Research
    * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
    * `'JPL'`: Jet Propulsion Laboratory
    * `'CNES'`: French Centre National D'Etudes Spatiales
    * `'GRAZ'`: Institute of Geodesy from GRAZ University of Technology
    * `'COSTG'`: Combination Service for Time-variable Gravity Fields
    * `'Swarm'`: Time-variable gravity data from Swarm satellites
3. `DREL`: GRACE/GRACE-FO/Swarm data release (RL04, RL05, RL06)

#### Keyword arguments
- `DSET`: GRACE/GRACE-FO/Swarm dataset
    * `'GAA'`: non-tidal atmospheric correction
    * `'GAB'`: non-tidal oceanic correction
    * `'GAC'`: combined non-tidal atmospheric and oceanic correction
    * `'GAD'`: ocean bottom pressure product
    * `'GSM'`: corrected monthly static field product

#### Returns
- `start`: First month in a GRACE/GRACE-FO dataset
- `end`: Last month in a GRACE/GRACE-FO dataset
- `missing`: missing months in a GRACE/GRACE-FO dataset
- `months`: all available months in a GRACE/GRACE-FO dataset
- `time`: center dates of all available months in a GRACE/GRACE-FO dataset
