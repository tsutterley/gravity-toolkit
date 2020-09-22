regress_grace_maps.py
=====================

 - Reads in GRACE/GRACE-FO spatial files from `grace_spatial_maps.py` and fits a regression model at each grid point

#### Calling Sequence
```bash
python regress_grace_maps.py --order=1 --cycles=0.5,1.0 --mode=0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/regress_grace_maps.py)

#### Inputs
   parameter file containing specific variables for the analysis

#### Command Line Options
 - `-P X`, `--np=X`: Run in parallel with X number of processes
 - `-S X`, `--start=X`: starting GRACE/GRACE-FO month for time series regression
 - `-E X`, `--end=X`: ending GRACE/GRACE-FO month for time series regression
 - `--order=X`: regression fit polynomial order
 - `--cycles=X`: list of regression fit cyclical terms as wavelength in decimal years
 - `-V`, `--verbose`: verbose output of processing run
 - `-M X`, `--mode=X`: permissions mode of output files
 - `-l`, `--log`: output log file for each job
