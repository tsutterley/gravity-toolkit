hdf5_stokes.py
==============

- Writes spherical harmonic coefficients to HDF5 files

#### Calling Sequence
```python
from gravity_toolkit.hdf5_stokes import hdf5_stokes
hdf5_stokes(clm, slm, linp, minp, tinp, month, FILENAME=output_HDF5_file)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/hdf5_stokes.py)

#### Arguments
- `clm`: Cosine spherical harmonic coefficients
- `slm`: Sine spherical harmonic coefficients
- `linp`: spherical harmonic degree (l)
- `minp`: spherical harmonic order (m)
- `tinp`: date of measurement
- `month`: GRACE/GRACE-FO month

#### Keyword arguments
- `FILENAME`: output filename HDF5
- `UNITS`: spherical harmonic units
- `TIME_UNITS`: time variable units
- `TIME_LONGNAME`: time variable description
- `MONTHS_NAME`: name of months variable within HDF5 file
- `MONTHS_UNITS`: months variable units
- `MONTHS_LONGNAME`: months variable description
- `TITLE`: description attribute of dataset
- `REFERENCE`: reference attribute of dataset
- `CLOBBER`: will overwrite an existing HDF5 file
- `VERBOSE`: will print to screen the HDF5 structure parameters
- `DATE`: harmonics have date information
