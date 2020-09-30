convert_harmonics.py
====================

 - Converts a file from the spatial domain into the spherical harmonic domain

#### Calling Sequence
```bash
python convert_harmonics.py --format netCDF4 --lmax 60 --units 1 input_file output_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/convert_harmonics.py)

#### Command Line Options
 - `--help`: list the command line options
 - `--lmax X`: maximum spherical harmonic degree
 - `--mmax X`: maximum spherical harmonic order
 - `--reference X`: Reference frame for load love numbers
 - `-R X`, `--radius X`: Gaussian smoothing radius (km)
 - `-D`, `--destripe`: use a decorrelation filter (destriping filter)
 - `-U X`, `--units X`: output units
    1. cm of water thickness (cm.w.e.)
    2. Gigatonnes (Gt)
    3. mm of water thickness (kg/m^2)
 - `-S X`, `--spacing X`: spatial resolution of output data (dlon,dlat)
 - `-I X`, `--interval X`: output grid interval
    1. (0:360, 90:-90)
    2. (degree spacing/2)
 - `--missing`: input spatial fields have missing values
 - `--fill-value X`: set fill_value for input spatial fields
 - `--header X`: number of header rows to skip in input ascii files
 - `--delimiter X`: delimiter in input ascii files
 - `-F X`, `--format X`: input and output data format
    * `'ascii'`
    * `'netCDF4'`
    * `'HDF5'`
 - `-V`, `--verbose`: verbose output of processing run
 - `-M X`, `--mode X`: Permissions mode of the files created
