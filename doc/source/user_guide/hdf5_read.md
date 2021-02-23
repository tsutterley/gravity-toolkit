hdf5_read.py
============

- Reads spatial data from HDF5 files

#### Calling Sequence
```python
from gravity_toolkit.hdf5_read import hdf5_read
file_inp = hdf5_read(filename, DATE=True, VERBOSE=False)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/hdf5_read.py)

#### Inputs
- `filename`: HDF5 file to be opened and read

#### Options
- `DATE`: HDF5 file has date information
- `VERBOSE`: will print to screen the HDF5 structure parameters
- `VARNAME`: z variable name in HDF5 file
- `LONNAME`: longitude variable name in HDF5 file
- `LATNAME`: latitude variable name in HDF5 file
- `TIMENAME`: time variable name in HDF5 file
- `COMPRESSION`: netCDF4 file is compressed or streaming as bytes
    * `'gzip'`
    * `'zip'`
    * `'bytes'`

#### Outputs
- `data`: z value of dataset
- `lon`: longitudinal array
- `lat`: latitudinal array
- `time`: time value of dataset (if specified by DATE)
- `attributes`: HDF5 attributes (for variables and title)
