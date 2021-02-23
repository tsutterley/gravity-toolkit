grace_spatial_maps.py
=====================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports monthly spatial fields
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms
- Converts data to specified units and performs a spherical harmonic summation to convert to the spatial domain

#### Calling Sequence
```bash
python grace_spatial_maps.py --mode 0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py)

#### Inputs
parameter file containing specific variables for the analysis

#### Command Line Options
- `-P X`, `--np X`: Run in parallel with X number of processes
- `-D X`, `--directory X`: Working data directory
- `-n X`, `--love X`: Load Love numbers dataset
     * `0`: Han and Wahr (1995) values from PREM
     * `1`: Gegout (2005) values from PREM
     * `2`: Wang et al. (2012) values from PREM
- `-r X`, `--reference X`: Reference frame for load love numbers
     * `'CF'`: Center of Surface Figure (default)
     * `'CM'`: Center of Mass of Earth System
     * `'CE'`: Center of Mass of Solid Earth
- `-V`, `--verbose`: verbose output of processing run
- `-M X`, `--mode X`: permissions mode of output files
- `-l`, `--log`: output log file for each job
