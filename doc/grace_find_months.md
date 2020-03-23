grace_find_months.py
====================

 - Finds the months available for a GRACE/GRACE-FO product
 - Finds the all months missing from the product

#### Calling Sequence
```
from gravity_toolkit.grace_find_months import grace_find_months
grace_months = grace_find_months(base_dir, PROC, DREL, DSET=DSET)
```

#### Inputs
 - `base_dir`: Working data directory for GRACE/GRACE-FO data
 - `PROC`: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)
 - `DREL`: GRACE/GRACE-FO data release (RL04, RL05, RL06)

#### Options
 - `DSET`: GRACE dataset (GSM, GAC, GAD, GAB, GAA)

#### Outputs
 - `start`: First month in a GRACE/GRACE-FO dataset
 - `end`: Last month in a GRACE/GRACE-FO dataset
 - `missing`: missing months in a GRACE/GRACE-FO dataset
 - `months`: all available months in a GRACE/GRACE-FO dataset
 - `time`: center dates of all available months in a GRACE/GRACE-FO dataset
