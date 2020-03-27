ncdf_read_stokes.py
===================

- Reads spherical harmonic data from netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
file_inp = ncdf_read_stokes(filename, DATE=True, VERBOSE=False)
```

#### Inputs
 - `filename`: netCDF4 file to be opened and read

#### Options
 - `DATE`: netCDF4 file has date information
 - `VERBOSE`: will print to screen the netCDF4 structure parameters

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients
 - `slm`: Sine spherical harmonic coefficients
 - `l`: spherical harmonic degree
 - `m`: spherical harmonic order
 - `time`: time of measurement (if specified by DATE)
 - `month`: GRACE/GRACE-FO month (if specified by DATE)
 - `attributes`: netCDF4 attributes for:
    * spherical harmonics (`clm`,`slm`)
    * variables (`l`,`m`,`time`,`month`)
    * title
