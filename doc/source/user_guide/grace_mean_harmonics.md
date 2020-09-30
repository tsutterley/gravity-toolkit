grace_mean_harmonics.py
=======================

 - Calculates the temporal mean of the GRACE/GRACE-FO spherical harmonics for a specified date range
 - Used to estimate the static gravitational field over a given date rage

#### Calling Sequence
```bash
python grace_mean_harmonics.py --mode 0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_mean_harmonics.py)

#### Inputs
   parameter file containing specific variables for the analysis

#### Command Line Options
 - `-V`, `--verbose`: verbose output of processing run
 - `-M X`, `--mode X`: permissions mode of output files
 - `-l`, `--log`: output log file for each job
