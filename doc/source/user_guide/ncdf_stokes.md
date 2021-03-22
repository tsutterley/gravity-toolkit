ncdf_stokes.py
==============

- Writes spherical harmonic coefficients to netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_stokes import ncdf_stokes
ncdf_stokes(clm1, slm1, linp, minp, tinp, month, FILENAME=output_netcdf4_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ncdf_stokes.py)

#### Arguments
- `clm1`: cosine spherical harmonic coefficients
- `slm1`: sine spherical harmonic coefficients
- `linp`: spherical harmonic degree (l)
- `minp`: spherical harmonic order (m)
- `tinp`: date of measurement
- `month`: GRACE/GRACE-FO month

#### Keyword arguments
- `FILENAME`: output filename netCDF4
- `UNITS`: spherical harmonic units
- `TIME_UNITS`: time variable units
- `TIME_LONGNAME`: time variable description
- `MONTHS_NAME`: name of months variable within netCDF4 file
- `MONTHS_UNITS`: months variable units
- `MONTHS_LONGNAME`: months variable description
- `TITLE`: title attribute of dataset
- `REFERENCE`: reference attribute of dataset
- `CLOBBER`: will overwrite an existing netCDF4 file
- `VERBOSE`: will print to screen the netCDF4 structure parameters
- `DATE`: harmonics have date information
