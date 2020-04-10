hdf5_write.py
=============

 - Writes spatial data to HDF5 files    

#### Calling Sequence
```python
from gravity_toolkit.hdf5_write import hdf5_write
hdf5_write(data, lon, lat, tim, FILENAME=output_netcdf4_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/gravity_toolkit/hdf5_write.py)

#### Inputs
 - `data`: z data
 - `lon`: longitude array
 - `lat`: latitude array
 - `tim`: time array

#### Options
 - `FILENAME`: output HDF5 filename
 - `VARNAME`: z variable name in HDF5 file
 - `LONNAME`: longitude variable name in HDF5 file
 - `LATNAME`: latitude variable name in HDF5 file
 - `UNITS`: z variable units
 - `LONGNAME`: z variable description
 - `FILL_VALUE`: missing value for z variable
 - `TIME_UNITS`: time variable units
 - `TIME_LONGNAME`: time variable description
 - `TITLE`: title attribute of dataset
 - `CLOBBER`: will overwrite an existing HDF5 file
 - `VERBOSE`: will print to screen the HDF5 structure parameters
