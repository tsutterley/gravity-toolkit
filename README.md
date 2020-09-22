read-GRACE-harmonics
====================

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/read-grace-harmonics/badge/?version=latest)](https://read-grace-harmonics.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/read-GRACE-harmonics/main)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/read-GRACE-harmonics/main)

Python tools for obtaining and working with Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

#### Resources
- [NASA GRACE mission site](https://www.nasa.gov/mission_pages/Grace/index.html)
- [NASA GRACE-FO mission site](https://www.nasa.gov/missions/grace-fo)
- [JPL GRACE Tellus site](https://grace.jpl.nasa.gov/)
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)

#### Programs
- [`aod1b_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/aod1b_geocenter.md) - Creates monthly files of geocenter variations due to non-tidal atmospheric or oceanic variation at 6-hour intervals
- [`aod1b_oblateness`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/aod1b_oblateness.md) - Creates monthly files of oblateness variations due to non-tidal atmospheric or oceanic variation at 6-hour intervals
- [`clenshaw_summation`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/clenshaw_summation.md) - Returns the spatial field for a series of spherical harmonics at a sequence of ungridded points
- [`combine_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/combine_harmonics.md) - Converts a file from the spherical harmonic domain into the spatial domain
- [`convert_calendar_decimal`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/convert_calendar_decimal.md) - Converts from calendar date into decimal years taking into account leap years
- [`convert_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/convert_harmonics.md) - Converts a file from the spatial domain into the spherical harmonic domain
- [`convert_julian`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/convert_julian.md) - Return the calendar date and time given Julian date
- [`degree_amplitude`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/degree_amplitude.md) - Calculates the amplitude of each spherical harmonic degree
- [`destripe_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/destripe_harmonics.md) - Filters spherical harmonic coefficients for correlated "striping" errors
- [`gauss_weights`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gauss_weights.md) - Computes the Gaussian weights as a function of degree
- [`gen_disc_load`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gen_disc_load.md) - Calculates gravitational spherical harmonic coefficients for a uniform disc load
- [`gen_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gen_harmonics.md) - Calculates the spherical harmonic coefficients of a spatial field
- [`gen_point_load`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gen_point_load.md) - Calculates gravitational spherical harmonic coefficients for point masses
- [`gen_spherical_cap`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gen_spherical_cap.md) - Calculates gravitational spherical harmonic coefficients for a spherical cap
- [`gen_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gen_stokes.md) - Returns a series of spherical harmonics for an input spatial field
- [`geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/geocenter.md) - Converts degree 1 spherical harmonic coefficients to geocenter variations
- [`gfz_isdc_dealiasing_ftp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gfz_isdc_dealiasing_ftp.md) - Syncs GRACE Level-1b dealiasing products from the GFZ Information System and Data Center (ISDC)
- [`gfz_isdc_grace_ftp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/gfz_isdc_grace_ftp.md) - Syncs GRACE/GRACE-FO and auxiliary data from the GFZ Information System and Data Center (ISDC)
- [`grace_date`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_date.md) - Calculates dates of each GRACE/GRACE-FO file and assigns the month number
- [`grace_find_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_find_months.md) - Finds the months available for a GRACE/GRACE-FO product
- [`grace_input_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_input_months.md) - Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range
- [`grace_mean_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_mean_harmonics.md) - Calculates the temporal mean of the GRACE/GRACE-FO spherical harmonics for a specified date range
- [`grace_months_index`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_months_index.md) - Creates an index of dates for all GRACE/GRACE-FO processing centers
- [`grace_spatial_error`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_spatial_error.md) - Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports spatial error field
- [`grace_spatial_maps`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/grace_spatial_maps.md) - Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports monthly spatial fields
- [`harmonic_summation`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/harmonic_summation.md) - Returns the spatial field for a series of spherical harmonics
- [`harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/harmonics.rst) - Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data
- [`hdf5_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/hdf5_read_stokes.md) - Reads spherical harmonic data from HDF5 files
- [`hdf5_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/hdf5_read.md) - Reads spatial data from HDF5 files
- [`hdf5_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/hdf5_stokes.md) - Writes spherical harmonic data to HDF5 files
- [`hdf5_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/hdf5_write.md) - Writes spatial data to HDF5 files
- [`legendre`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/legendre.md) - Computes associated Legendre functions for a specific spherical harmonic degree
- [`legendre_polynomials`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/legendre_polynomials.md) - Computes fully normalized Legendre polynomials and their first derivative
- [`ncdf_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/ncdf_read_stokes.md) - Reads spherical harmonic data from netCDF4 files
- [`ncdf_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/ncdf_read.md) - Reads spatial data from netCDF4 files
- [`ncdf_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/ncdf_stokes.md) - Writes spherical harmonic data to netCDF4 files
- [`ncdf_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/ncdf_write.md) - Writes spatial data to netCDF4 files
- [`ocean_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/ocean_stokes.md) - Reads a land-sea mask and converts to a series of spherical harmonics
- [`plm_columbo`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/plm_columbo.md) - Computes fully-normalized associated Legendre Polynomials using the Colombo (1981) recursion relation
- [`plm_holmes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/plm_holmes.md) - Computes fully-normalized associated Legendre Polynomials using the Holmes and Featherstone (2002) recursion relation
- [`plm_mohlenkamp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/plm_mohlenkamp.md) - Computes fully-normalized associated Legendre Polynomials using Martin Mohlenkamp's recursion relation
- [`podaac_grace_sync`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/podaac_grace_sync.md) - Syncs GRACE/GRACE-FO and auxiliary data from the NASA JPL PO.DAAC Drive Server
- [`podaac_webdav`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/podaac_webdav.md) - Retrieves and prints a user's PO.DAAC WebDAV credentials
- [`read_CSR_monthly_6x1`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_CSR_monthly_6x1.md) - Reads the monthly low-degree spherical harmonic data files from satellite laser ranging (SLR)
- [`read_GIA_model`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_GIA_model.md) - Reads Glacial Isostatic Adjustment (GIA) files and extracts spherical harmonic data
- [`read_GRACE_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_GRACE_harmonics.md) - Reads GRACE/GRACE-FO files and extracts spherical harmonic data
- [`read_ICGEM_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_ICGEM_harmonics.md) - Reads gfc files and extracts gravity model spherical harmonics from the GFZ ICGEM
- [`read_love_numbers`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_love_numbers.md) - Reads sets of load Love numbers output from the Preliminary Reference Earth Model (PREM)
- [`read_SLR_C20`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_SLR_C20.md) - Reads monthly oblateness spherical harmonic data files from satellite laser ranging (SLR)
- [`read_SLR_C30`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_SLR_C30.md) - Reads monthly degree 3 zonal spherical harmonic data files from satellite laser ranging (SLR)
- [`read_SLR_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_SLR_geocenter.md) - Reads monthly geocenter spherical harmonic data files from satellite laser ranging (SLR)
- [`read_tellus_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/read_tellus_geocenter.md) - Reads monthly geocenter spherical harmonic data files from GRACE Tellus Technical Notes
- [`regress_grace_maps`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/regress_grace_maps.md) - Reads in GRACE/GRACE-FO spatial files and fits a regression model at each grid point
- [`run_grace_date`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/run_grace_date.md) - Wrapper program for running GRACE date and months programs
- [`spatial`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/spatial.rst) - Spatial data class for reading, writing and processing spatial data
- [`tsamplitude`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/tsamplitude.md) - Calculate the amplitude and phase of a harmonic function from calculated sine and cosine of a series of measurements
- [`tsregress`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/tsregress.md) - Fits a synthetic signal to data over a time period by least-squares or weighted least-squares
- [`tssmooth`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/tssmooth.md) - Computes a moving average of a time-series
- [`units`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/units.rst) - Class for converting GRACE/GRACE-FO Level-2 data to specific units
- [`utilities.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/doc/source/user_guide/utilities.rst) - Download and management utilities for syncing time and auxiliary files  

#### Dependencies
- [numpy: Scientific Computing Tools For Python](https://www.numpy.org)
- [scipy: Scientific Tools for Python](https://docs.scipy.org/doc/)
- [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)
- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)
- [future: Compatibility layer between Python 2 and Python 3](https://python-future.org/)
- [matplotlib: Python 2D plotting library](https://matplotlib.org/)
- [cartopy: Python package designed for geospatial data processing](https://scitools.org.uk/cartopy/docs/latest/)
- [netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/)
- [h5py: Python interface for Hierarchal Data Format 5 (HDF5)](https://www.h5py.org/)
- [read-GRACE-geocenter: Python reader for GRACE/GRACE-FO geocenter data](https://github.com/tsutterley/read-GRACE-geocenter/)

#### References
T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self‐Consistent Ice Mass Balance
and Regional Sea Level From Time‐Variable Gravity", *Earth and Space Science*, 7,
(2020). [doi:10.1029/2019EA000860](https://doi.org/10.1029/2019EA000860)

T. C. Sutterley and I. Velicogna, "Improved estimates of geocenter variability
from time-variable gravity and ocean model outputs", *Remote Sensing*, 11(18),
2108, (2019). [doi:10.3390/rs11182108](https://doi.org/10.3390/rs11182108)

J. Wahr, S. C. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates",
*Geophysical Research Letters*, 33(6), L06401, (2006).
[doi: 10.1029/2005GL025305](https://doi.org/10.1029/2005GL025305)

J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity
field: Hydrological and oceanic effects and their possible detection using
GRACE", *Journal of Geophysical Research: Solid Earth*, 103(B12), (1998).
[doi: 10.1029/98JB02844](https://doi.org/10.1029/98JB02844)

D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified
earth, and a further analysis of postglacial rebound", *Geophysical Journal
International*, 120(2), (1995).
[doi: 10.1111/j.1365-246X.1995.tb01819.x](https://doi.org/10.1111/j.1365-246X.1995.tb01819.x)

#### Data Repositories
T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Ice Mass and Regional Sea Level
Estimates from Time-Variable Gravity", (2020).
[doi:10.6084/m9.figshare.9702338](https://doi.org/10.6084/m9.figshare.9702338)

T. C. Sutterley and I. Velicogna, "Geocenter Estimates from Time-Variable
Gravity and Ocean Model Outputs", (2019).
[doi:10.6084/m9.figshare.7388540](https://doi.org/10.6084/m9.figshare.7388540)

#### Download
The program homepage is:  
https://github.com/tsutterley/read-GRACE-harmonics  
A zip archive of the latest version is available directly at:  
https://github.com/tsutterley/read-GRACE-harmonics/archive/main.zip  

#### Disclaimer
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL), the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.

#### License
The content of this project is licensed under the [Creative Commons Attribution 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/) and the source code is licensed under the [MIT license](LICENSE).
