read-GRACE-harmonics
====================

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/read-grace-harmonics/badge/?version=latest)](https://read-grace-harmonics.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/read-GRACE-harmonics/master)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/read-GRACE-harmonics/master)

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
- [`aod1b_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/aod1b_geocenter.md) - Creates monthly files of geocenter variations due to non-tidal atmospheric or oceanic variation at 6-hour intervals  
- [`combine_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/combine_harmonics.md) - Returns the spatial field for a series of spherical harmonics  
- [`convert_calendar_decimal`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/convert_calendar_decimal.md) - Converts from calendar date into decimal years taking into account leap years    
- [`convert_julian`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/convert_julian.md) - Return the calendar date and time given Julian date  
- [`destripe_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/destripe_harmonics.md) - Filters spherical harmonic coefficients for correlated "striping" errors  
- [`gauss_weights`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/gauss_weights.md) - Computes the Gaussian weights as a function of degree  
- [`gen_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/gen_stokes.md) - Returns a series of spherical harmonics for an input spatial field  
- [`geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/geocenter.md) - Converts degree 1 spherical harmonic coefficients to geocenter variations  
- [`gfz_isdc_grace_ftp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/gfz_isdc_grace_ftp.md) - Syncs GRACE/GRACE-FO and auxiliary data from the GFZ Information System and Data Center (ISDC)  
- [`grace_date`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/grace_date.md) - Calculates dates of each GRACE/GRACE-FO file and assigns the month number  
- [`grace_find_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/grace_find_months.md) - Finds the months available for a GRACE/GRACE-FO product  
- [`grace_input_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/grace_input_months.md) - Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range  
- [`harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/harmonics.rst) - Spherical harmonic data class for processing GRACE/GRACE-FO Level-2 data
- [`hdf5_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/hdf5_read_stokes.md) - Reads spherical harmonic data from HDF5 files  
- [`hdf5_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/hdf5_read.md) - Reads spatial data from HDF5 files  
- [`hdf5_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/hdf5_stokes.md) - Writes spherical harmonic data to HDF5 files  
- [`hdf5_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/hdf5_write.md) - Writes spatial data to HDF5 files  
- [`ncdf_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/ncdf_read_stokes.md) - Reads spherical harmonic data from netCDF4 files  
- [`ncdf_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/ncdf_read.md) - Reads spatial data from netCDF4 files  
- [`ncdf_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/ncdf_stokes.md) - Writes spherical harmonic data to netCDF4 files  
- [`ncdf_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/ncdf_write.md) - Writes spatial data to netCDF4 files  
- [`ocean_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/ocean_stokes.md) - Reads a land-sea mask and converts to a series of spherical harmonics  
- [`plm_columbo`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/plm_columbo.md) - Computes fully-normalized associated Legendre Polynomials using the Colombo (1981) recursion relation  
- [`plm_holmes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/plm_holmes.md) - Computes fully-normalized associated Legendre Polynomials using the Holmes and Featherstone (2002) recursion relation  
- [`plm_mohlenkamp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/plm_mohlenkamp.md) - Computes fully-normalized associated Legendre Polynomials using Martin Mohlenkamp's recursion relation  
- [`podaac_grace_sync`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/podaac_grace_sync.md) - Syncs GRACE/GRACE-FO and auxiliary data from the NASA JPL PO.DAAC Drive Server  
- [`read_CSR_monthly_6x1`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_CSR_monthly_6x1.md) - Reads the monthly low-degree spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_GIA_model`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_GIA_model.md) - Reads Glacial Isostatic Adjustment (GIA) files and extracts spherical harmonic data  
- [`read_GRACE_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_GRACE_harmonics.md) - Reads GRACE/GRACE-FO files and extracts spherical harmonic data  
- [`read_ICGEM_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_ICGEM_harmonics.md) - Reads gfc files and extracts gravity model spherical harmonics from the GFZ ICGEM  
- [`read_love_numbers`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_love_numbers.md) - Reads sets of load Love numbers output from the Preliminary Reference Earth Model (PREM)  
- [`read_SLR_C20`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_SLR_C20.md) - Reads monthly oblateness spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_SLR_C30`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_SLR_C30.md) - Reads monthly degree 3 zonal spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_SLR_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_SLR_geocenter.md) - Reads monthly geocenter spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_tellus_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/read_tellus_geocenter.md) - Reads monthly geocenter spherical harmonic data files from GRACE Tellus Technical Notes  
- [`units`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/source/user_guide/units.rst) - Class for converting GRACE/GRACE-FO Level-2 data to specific units

#### Dependencies
- [numpy: Scientific Computing Tools For Python](http://www.numpy.org)  
- [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)  
- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)  
- [future: Compatibility layer between Python 2 and Python 3](http://python-future.org/)  
- [matplotlib: Python 2D plotting library](http://matplotlib.org/)  
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

#### Download
The program homepage is:   
https://github.com/tsutterley/read-GRACE-harmonics   
A zip archive of the latest version is available directly at:    
https://github.com/tsutterley/read-GRACE-harmonics/archive/master.zip  

#### Disclaimer  
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL), the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.  
