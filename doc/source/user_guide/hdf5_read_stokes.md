hdf5_read_stokes.py
===================

 - Reads spherical harmonic data from HDF5 files

#### Calling Sequence
```
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes
file_inp = hdf5_read_stokes(filename, DATE=True, VERBOSE=False)
```

#### Inputs
 - `filename`: HDF5 file to be opened and read

#### Options
 - `DATE`: HDF5 file has date information
 - `VERBOSE`: will print to screen the HDF5 structure parameters

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients
 - `slm`: Sine spherical harmonic coefficients
 - `l`: spherical harmonic degree
 - `m`: spherical harmonic order
 - `time`: time of measurement (if specified by DATE)
 - `month`: GRACE/GRACE-FO month (if specified by DATE)
 - `attributes`: HDF5 attributes for:
     - spherical harmonics (`clm`,`slm`)
     - variables (`l`,`m`,`time`,`month`)
     - title
