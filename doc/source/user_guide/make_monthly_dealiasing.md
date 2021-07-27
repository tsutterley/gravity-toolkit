make_monthly_dealiasing.py
==========================

- Reads GRACE/GRACE-FO level-1b dealiasing data files for a specific product and outputs monthly the mean for a specific GRACE/GRACE-FO processing center and data release
    * `'GAA'`: atmospheric loading from ECMWF
    * `'GAB'`: oceanic loading from OMCT/MPIOM
    * `'GAC'`: global atmospheric and oceanic loading
    * `'GAD'`: ocean bottom pressure from OMCT/MPIOM
- Creates monthly files of oblateness variations at 3 or 6-hour intervals

#### Calling Sequence
```bash
python make_monthly_dealiasing.py --release RL06 --product GAD
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/make_monthly_dealiasing.py)

#### Command Line Options
- `-D X`, `--directory X`: Working Data Directory
- `-c X`, `--center X`: GRACE/GRACE-FO Processing Center for dates
- `-r X`, `--release X`: GRACE/GRACE-FO Data Release (RL05 or RL06)
- `-p X`, `--product X`: GRACE/GRACE-FO dealiasing product (GAA, GAB, GAC, GAD)
- `-L X`, `--lmax X`: Maximum spherical harmonic degree and order for output
- `-F X`, `--format X`: Output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
    * `'SHM'`
- `-C`, `--clobber`: Overwrite existing data
- `-M X`, `--mode X`: Permission mode of directories and files
- `-V`, `--verbose`: Output information for each output file
