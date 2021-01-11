ocean_stokes.py
===============

 - Reads a land-sea mask and converts to a series of spherical harmonics

#### Calling Sequence
```python
from gravity_toolkit.ocean_stokes import ocean_stokes
ocean_Ylms = ocean_stokes(LANDMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ocean_stokes.py)

#### Inputs
 - `LANDMASK`: mask file from [Sutterley et al. (2020)](https://doi.org/10.6084/m9.figshare.9702338)
     * updated 1,0.5 and 0.25 degree masks from [ORNL as part of ISLSCP](https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html)
 - `LMAX`:  maximum spherical harmonic degree of the output harmonics

#### Options
 - `MMAX`: maximum spherical harmonic order of the output harmonics
 - `LOVE`: input load Love numbers up to degree LMAX (hl,kl,ll)
 - `VARNAME`: variable name for mask in netCDF4 file
 - `SIMPLIFY`: simplify land mask by removing isolated points

#### Outputs
 - `clm`: Cosine spherical harmonic coefficients (geodesy normalization)
 - `slm`: Sine spherical harmonic coefficients (geodesy normalization)
 - `l`: spherical harmonic degree to LMAX
 - `m`: spherical harmonic order to MMAX
