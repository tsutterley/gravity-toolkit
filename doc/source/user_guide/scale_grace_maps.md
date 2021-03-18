scale_grace_maps.py
===================

- Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports scaled spatial fields
- Correct spherical harmonics with the specified GIA model group
- Filters and smooths data with specified processing algorithms
- Converts data to centimeters water equivalent, performs a spherical harmonic summation to convert to the spatial domain
- Scales the spatial fields following [Landerer and Swenson (2012)](https://doi.org/10.1029/2011WR011453)
- Calculates the scaled spatial error field following [Wahr et al. (2006)](https://doi.org/10.1029/2005GL025305)

#### Calling Sequence
```bash
python scale_grace_maps.py --mode 0o775 parameter_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/scale_grace_maps.py)

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
