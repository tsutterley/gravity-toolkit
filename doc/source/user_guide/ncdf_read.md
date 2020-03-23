ncdf_read.py
============

 - Reads spatial data from COARDS-compliant netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_read import ncdf_read
file_inp = ncdf_read(filename, DATE=True, VERBOSE=False)
```

#### Inputs
 - `filename`: netCDF4 file to be opened and read

#### Options
 - `DATE`: netCDF4 file has date information
 - `MISSING`: netCDF4 variables have missing values
 - `VERBOSE`: will print to screen the netCDF4 structure parameters
 - `VARNAME`: z variable name in netCDF4 file
 - `LONNAME`: longitude variable name in netCDF4 file
 - `LATNAME`: latitude variable name in netCDF4 file
 - `TIMENAME`: time variable name in netCDF4 file
 - `ATTRIBUTES`: netCDF4 variables contain attribute parameters
 - `TITLE`: netCDF4 file contains description attribute parameter

#### Outputs
 - `data`: z value of dataset
 - `lon`: longitudinal array
 - `lat`: latitudinal array
 - `time`: time value of dataset (if specified by DATE)
 - `attributes`: netCDF4 attributes (for variables and title)
