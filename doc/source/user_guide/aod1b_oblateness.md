aod1b_oblateness.py
===================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product
    - `atm`: atmospheric loading from ECMWF
    - `ocn`: oceanic loading from OMCT/MPIOM
    - `glo`: global atmospheric and oceanic loading
    - `oba`: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of oblateness variations at 6-hour intervals

#### Calling Sequence
```python
from gravity_toolkit.aod1b_oblateness import aod1b_oblateness
aod1b_oblateness(base_dir, DREL='RL06', DSET='glo', CLOBBER=True)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/aod1b_oblateness.py)

#### Inputs
 1. `base_dir`: working data directory  

#### Options
 - `DREL`: GRACE/GRACE-FO data release (RL05 or RL06)  
 - `DSET`: GRACE/GRACE-FO dataset (atm, ocn, glo, oba)  
 - `CLOBBER`: overwrite existing data  
 - `MODE`: Permission mode of directories and files  
 - `VERBOSE`: Output information for each output file  
