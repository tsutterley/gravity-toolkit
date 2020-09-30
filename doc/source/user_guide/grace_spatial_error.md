grace_spatial_error.py
======================

 - Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports spatial error field following [Wahr et al. (2006)](https://doi.org/10.1029/2005GL025305)
 - Filters and smooths data with specified processing algorithms
 - Converts data to specified units and performs a spherical harmonic summation to convert error field to the spatial domain

#### Calling Sequence
```bash
python grace_spatial_error.py --mode 0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_error.py)

#### Inputs
   parameter file containing specific variables for the analysis

#### Command Line Options
 - `-P X`, `--np X`: Run in parallel with X number of processes
 - `-V`, `--verbose`: verbose output of processing run
 - `-M X`, `--mode X`: permissions mode of output files
 - `-l`, `--log`: output log file for each job
