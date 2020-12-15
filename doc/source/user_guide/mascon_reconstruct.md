mascon_reconstruct.py
=====================

 - Calculates the equivalent spherical harmonics from a mascon time series

#### Calling Sequence
```bash
python mascon_reconstruct.py --mode 0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/mascon_reconstruct.py)

#### Inputs
   parameter file containing specific variables for the analysis

#### Command Line Options
 - `-P X`, `--np X`: Run in parallel with X number of processes
 - `-n X`, `--love X`: Load Love numbers dataset
      * `0`: Han and Wahr (1995) values from PREM
      * `1`: Gegout (2005) values from PREM
      * `2`: Wang et al. (2012) values from PREM
 - `-r X`, `--reference X`: Reference frame for load love numbers
      * `'CF'`: Center of Surface Figure (default)
      * `'CM'`: Center of Mass of Earth System
      * `'CE'`: Center of Mass of Solid Earth
 - `-l`, `--log`: output log file for each job
 - `-M X`, `--mode X`: permissions mode of output files
