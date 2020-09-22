grace_months_index.py
=====================

 - Reads GRACE/GRACE-FO date files from `grace_date.py`  
 - Creates an index of dates for all GRACE/GRACE-FO processing centers  

#### Calling Sequence
```python
from gravity_toolkit.grace_months_index import grace_months_index
grace_months_index(base_dir, DREL=DREL, MODE=)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/grace_months_index.py)

#### Inputs  
 1. Working data directory for GRACE/GRACE-FO data

#### Options  
 - `DREL`: GRACE/GRACE-FO data release (RL04/RL05/RL06)
 - `MODE`: Permissions mode of output file

#### Outputs
 - dictionary of files mapped by GRACE/GRACE-FO month  
