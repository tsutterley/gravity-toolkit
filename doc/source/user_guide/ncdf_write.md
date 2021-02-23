ncdf_write.py
=============

- Writes spatial data to COARDS-compliant netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_write import ncdf_write
ncdf_write(data, lon, lat, tim, FILENAME=output_netcdf4_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ncdf_write.py)

#### Inputs
- `data`: z data
- `lon`: longitude array
- `lat`: latitude array
- `tim`: time array

#### Options
- `FILENAME`: output netCDF4 filename
- `VARNAME`: z variable name in netCDF4 file
- `LONNAME`: longitude variable name in netCDF4 file
- `LATNAME`: latitude variable name in netCDF4 file
- `UNITS`: z variable units
- `LONGNAME`: z variable description
- `FILL_VALUE`: missing value for z variable
- `TIME_UNITS`: time variable units
- `TIME_LONGNAME`: time variable description
- `TITLE`: title attribute of dataset
- `REFERENCE`: reference attribute of dataset
- `DATE`: data has date information
- `CLOBBER`: will overwrite an existing netCDF4 file
- `VERBOSE`: will print to screen the netCDF4 structure parameters
