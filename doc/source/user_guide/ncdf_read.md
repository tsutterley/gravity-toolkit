ncdf_read.py
============

- Reads spatial data from COARDS-compliant netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_read import ncdf_read
file_inp = ncdf_read(filename, DATE=True, VERBOSE=False)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ncdf_read.py)

#### Arguments
- `filename`: netCDF4 file to be opened and read

#### Keyword arguments
- `DATE`: netCDF4 file has date information
- `VARNAME`: z variable name in netCDF4 file
- `LONNAME`: longitude variable name in netCDF4 file
- `LATNAME`: latitude variable name in netCDF4 file
- `TIMENAME`: time variable name in netCDF4 file
- `COMPRESSION`: netCDF4 file is compressed or streaming as bytes
    * `'gzip'`
    * `'zip'`
    * `'bytes'`

#### Returns
- `data`: z value of dataset
- `lon`: longitudinal array
- `lat`: latitudinal array
- `time`: time value of dataset (if specified by `DATE`)
- `attributes`: netCDF4 attributes (for variables and title)
