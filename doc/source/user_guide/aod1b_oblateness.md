aod1b_oblateness.py
===================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product
    - `atm`: atmospheric loading from ECMWF
    - `ocn`: oceanic loading from OMCT/MPIOM
    - `glo`: global atmospheric and oceanic loading
    - `oba`: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of oblateness variations at 6-hour intervals

#### Calling Sequence
```bash
python aod1b_oblateness.py --release RL06 --product glo
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/aod1b_oblateness.py)

#### Command Line Options
 - `-D X`, `--directory X`: Working Data Directory
 - `-r X`, `--release X`: GRACE/GRACE-FO Data Release (RL05 or RL06)
 - `-p X`, `--product X`: GRACE/GRACE-FO dealiasing product (atm, ocn, glo, oba)
 - `-C`, `--clobber`: Overwrite existing data
 - `-M X`, `--mode X`: Permission mode of directories and files
 - `-V`, `--verbose`: Output information for each output file
