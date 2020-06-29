combine_harmonics.py
====================

 - Converts a file from the spherical harmonic domain into the spatial domain

#### Calling Sequence
```bash
python combine_harmonics --format=2 --lmax=60 --units=1 input_file output_file
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/scripts/combine_harmonics.py)

#### Command Line Options
 - `--help`: list the command line options
 - `--lmax=X`: maximum spherical harmonic degree
 - `--mmax=X`: maximum spherical harmonic order
 - `--love=X`: path to load love numbers file
 - `--reference=X`: Reference frame for load love numbers
 - `-R X`, `--radius=X`: Gaussian smoothing radius (km)
 - `-D`, `--destripe`: use a decorrelation filter (destriping filter)
 - `-U X`, `--units=X`: output units
    1. cm of water thickness
    2. mm of geoid height
    3. mm of elastic crustal deformation
    4. microGal gravitational perturbation
    5. Pa, equivalent surface pressure in Pascals
 - `-S X`, `--spacing=X`: spatial resolution of output data (dlon,dlat)
 - `-I X`, `--interval=X`: output grid interval
    1. (0:360, 90:-90)
    2. (degree spacing/2)
    3. non-global grid (set bounds with --bounds)
 - `-B X`, `--bounds=X`: bounding box for interval 3 (minlon,maxlon,minlat,maxlat)
 - `-O`, `--ocean`: redistribute total mass over the ocean
 - `--mask=X`: input land-sea function (netCDF4) with variable LSMASK as mask
 - `--mean=X`: mean file to remove from the harmonic data
 - `-F X`, `--format=X`: input and output data format
    1. ascii format
    2. netCDF4 format
    3. HDF5 format
 - `-V`, `--verbose`: verbose output of processing run
 - `-M X`, `--mode=X`: Permissions mode of the files created