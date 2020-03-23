hdf5_read.py
============

 - Reads spatial data from HDF5 files

#### Calling Sequence
```
from gravity_toolkit.hdf5_read import hdf5_read
file_inp = hdf5_read(filename, DATE=True, VERBOSE=False)
```

#### Inputs
 - `filename`: HDF5 file to be opened and read

#### Options
 - `DATE`: HDF5 file has date information
 - `MISSING`: HDF5 variables have missing values
 - `VERBOSE`: will print to screen the HDF5 structure parameters
 - `VARNAME`: z variable name in HDF5 file
 - `LONNAME`: longitude variable name in HDF5 file
 - `LATNAME`: latitude variable name in HDF5 file
 - `TIMENAME`: time variable name in HDF5 file
 - `ATTRIBUTES`: HDF5 variables contain attribute parameters
 - `TITLE`: HDF5 file contains description attribute parameter

#### Outputs
 - `data`: z value of dataset
 - `lon`: longitudinal array
 - `lat`: latitudinal array
 - `time`: time value of dataset (if specified by DATE)
 - `attributes`: HDF5 attributes (for variables and title)
