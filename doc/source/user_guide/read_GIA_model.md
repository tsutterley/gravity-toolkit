read_GIA_model.py
=================

- Reads Glacial Isostatic Adjustment (GIA) files that can come in various formats depending on the group
- Outputs spherical harmonics for the GIA rates and the GIA model parameters
- Can also output fully normalized harmonics to netCDF4 or HDF5 formats

#### Calling Sequence
```python
from gravity_toolkit.read_GIA_model import read_GIA_model
GIA_Ylms = read_GIA_model('Stokes.R2_65_.2_1.5_L120',GIA='IJ05-R2',LMAX=60)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_GIA_model.py)


#### Inputs:
1. `input_file`: GIA file to read (ascii, netCDF4 or HDF5)

#### Outputs:
- `GIA`: GIA model type to read and output
   * `'IJ05-R2'`: [Ivins R2 GIA Models](https://doi.org/10.1002/jgrb.50208)
   * `'W12a'`: [Whitehouse GIA Models](https://doi.org/10.1111/j.1365-246X.2012.05557.x)
   * `'SM09'`: [Simpson/Milne GIA Models](https://doi.org/10.1029/2010JB007776)
   * `'ICE6G'`: [ICE-6G GIA Models](https://doi.org/10.1002/2014JB011176)
   * `'Wu10'`: [Wu (2010) GIA Correction](https://doi.org/10.1038/ngeo938)
   * `'AW13-ICE6G'`: [Geruo A ICE-6G GIA Models](https://doi.org/10.1093/gji/ggs030)
   * `'Caron'`: [Caron JPL GIA Assimilation](https://doi.org/10.1002/2017GL076644)
   * `'ICE6G-D'`: [ICE-6G Version-D GIA Models](https://doi.org/10.1002/2016JB013844)
   * `'netCDF4'`: reformatted GIA in netCDF4 format
   * `'HDF5'`: reformatted GIA in HDF5 format
- `LMAX`: maximum degree of spherical harmonics
- `DATAFORM`: Spherical harmonic data output format
   * `None`: output only as variables
   * `'netCDF4'`: output to netCDF4 format (.nc)
   * `'HDF5'`: output to HDF5 format (.H5)
- `MODE`: permissions mode of output spherical harmonic files

#### Outputs:
- `clm`: cosine spherical harmonic of GIA rate
- `slm`: sine spherical harmonic of GIA rate
- `l`: spherical harmonic degree
- `m`: spherical harmonic order
- `title`: parameters of GIA model

#### References:
1. E. R. Ivins, T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer, and K. M. Simon, "Antarctic contribution to sea level rise observed by GRACE with improved GIA correction." Journal of Geophysical Research: Solid Earth, 118(6), 3126-3141 (2013). [https://doi.org/10.1002/jgrb.50208](https://doi.org/10.1002/jgrb.50208)

2. P. L. Whitehouse, M. J. Bentley, G. A. Milne, M. A. King, and I. D. Thomas, "A new glacial isostatic adjustment model for Antarctica: calibrated and tested using observations of relative sea-level change and present-day uplift rates." Geophysical Journal International, 190(3), 1464-1482 (2012). [https://doi.org/10.1111/j.1365-246X.2012.05557.x](https://doi.org/10.1111/j.1365-246X.2012.05557.x)

3. M. J. R. Simpson, L. Wake, G. A. Milne, and P. Huybrechts, "The influence of decadal- to millennial-scale ice mass changes on present-day vertical land motion in Greenland: Implications for the interpretation of GPS observations." Journal of Geophysical Research: Solid Earth, 116(B2), B02406 (2011). [https://doi.org/10.1029/2010JB007776](https://doi.org/10.1029/2010JB007776)

4. W. R. Peltier, D. F. Argus, and R. Drummond, "Space geodesy constrains ice age terminal deglaciation: The global ICE‐6G_C (VM5a) model." Journal of Geophysical Research: Solid Earth, 120(1), 450-487 (2015). [https://doi.org/10.1002/2014JB011176](https://doi.org/10.1002/2014JB011176)

5. X. Wu, M. B. Heflin, H. Schotman, B. L. A. Vermeersen, D. Dong, R. S. Gross, E. R. Ivins, A. W. Moore, S. E. Owen, "Simultaneous estimation of global present-day water transport and glacial isostatic adjustment." Nature Geoscience, 3(9), 642-646 (2010). [https://doi.org/10.1038/ngeo938](https://doi.org/10.1038/ngeo938)

6. G. A, J. Wahr, S. Zhong, "Computations of the viscoelastic response of a 3-D compressible Earth to surface loading: an application to Glacial Isostatic Adjustment in Antarctica and Canada." Geophysical Journal International, 192(2), 557-572 (2013). [https://doi.org/10.1093/gji/ggs030](https://doi.org/10.1093/gji/ggs030)

7. L. Caron, E. R. Ivins, E. Larour, S. Adhikari, J. Nilsson, and G. Blewitt, "GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science." Geophysical Research Letters, 45(5), 2203-2212 (2018). [https://doi.org/10.1002/2017GL076644](https://doi.org/10.1002/2017GL076644)

8. W. R. Peltier, D. F. Argus, and R. Drummond, "Comment on 'An Assessment of the ICE-6G_C (VM5a) Glacial Isostatic Adjustment Model' by Purcell et al." Journal of Geophysical Research: Solid Earth, 123(2), 2019-2028 (2018). [https://doi.org/10.1002/2016JB013844](https://doi.org/10.1002/2016JB013844)
