ncdf_read_stokes.py
===================

- Reads spherical harmonic data from netCDF4 files

#### Calling Sequence
```python
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
file_inp = ncdf_read_stokes(filename, DATE=True, VERBOSE=False)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ncdf_read_stokes.py)

#### Arguments
- `filename`: netCDF4 file to be opened and read

#### Keyword arguments
- `DATE`: netCDF4 file has date information
- `VERBOSE`: will print to screen the netCDF4 structure parameters
- `COMPRESSION`: netCDF4 file is compressed or streaming as bytes
    * `'gzip'`
    * `'zip'`
    * `'bytes'`

#### Returns
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
