read-GRACE-harmonics
====================

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/LICENSE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/read-GRACE-harmonics/master)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/read-GRACE-harmonics/master)

Python tools for obtaining and working with Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions  

#### Resources  
- [NASA GRACE mission site](http://www.nasa.gov/mission_pages/Grace/index.html)  
- [JPL GRACE Tellus site](http://grace.jpl.nasa.gov/)  
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)  
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)  
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)  

#### Programs
- [`aod1b_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/aod1b_geocenter.md) - Creates monthly files of geocenter variations due to non-tidal atmospheric or oceanic variation at 6-hour intervals  
- [`combine_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/combine_harmonics.md) - Returns the spatial field for a series of spherical harmonics  
- [`convert_calendar_decimal`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/convert_calendar_decimal.md) - Converts from calendar date into decimal years taking into account leap years    
- [`convert_julian`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/convert_julian.md) - Return the calendar date and time given Julian date  
- [`destripe_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/destripe_harmonics.md) - Filters spherical harmonic coefficients for correlated "striping" errors  
- [`gauss_weights`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/gauss_weights.md) - Computes the Gaussian weights as a function of degree  
- [`geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/geocenter.md) - Converts degree 1 spherical harmonic coefficients to geocenter variations  
- [`grace_date`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/grace_date.md) - Calculates dates of each GRACE/GRACE-FO file and assigns the month number  
- [`grace_find_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/grace_find_months.md) - Finds the months available for a GRACE/GRACE-FO product  
- [`grace_input_months`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/grace_input_months.md) - Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range  
- [`hdf5_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/hdf5_read_stokes.md) - Reads spherical harmonic data from HDF5 files  
- [`hdf5_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/hdf5_read.md) - Reads spatial data from HDF5 files  
- [`hdf5_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/hdf5_stokes.md) - Writes spherical harmonic data to HDF5 files  
- [`hdf5_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/hdf5_write.md) - Writes spatial data to HDF5 files  
- [`ncdf_read_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/ncdf_read_stokes.md) - Reads spherical harmonic data from netCDF4 files  
- [`ncdf_read`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/ncdf_read.md) - Reads spatial data from netCDF4 files  
- [`ncdf_stokes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/ncdf_stokes.md) - Writes spherical harmonic data to netCDF4 files  
- [`ncdf_write`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/ncdf_write.md) - Writes spatial data to netCDF4 files  
- [`plm_columbo`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/plm_columbo.md) - Computes fully-normalized associated Legendre Polynomials using the Colombo (1981) recursion relation  
- [`plm_holmes`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/plm_holmes.md) - Computes fully-normalized associated Legendre Polynomials using the Holmes and Featherstone (2002) recursion relation  
- [`plm_mohlenkamp`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/plm_mohlenkamp.md) - Computes fully-normalized associated Legendre Polynomials using Martin Mohlenkamp's recursion relation  
- [`read_CSR_monthly_6x1`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_CSR_monthly_6x1.md) - Reads the monthly low-degree spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_GRACE_harmonics`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_GRACE_harmonics.md) - Reads GRACE/GRACE-FO files and extracts spherical harmonic data  
- [`read_love_numbers`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_love_numbers.md) - Reads sets of load Love numbers output from the Preliminary Reference Earth Model (PREM)  
- [`read_SLR_C20`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_SLR_C20.md) - Reads monthly oblateness spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_SLR_C30`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_SLR_C30.md) - Reads monthly degree 3 zonal spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_SLR_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_SLR_geocenter.md) - Reads monthly geocenter spherical harmonic data files from satellite laser ranging (SLR)  
- [`read_tellus_geocenter`](https://github.com/tsutterley/read-GRACE-harmonics/blob/master/doc/read_tellus_geocenter.md) - Reads monthly geocenter spherical harmonic data files from GRACE Tellus Technical Notes  

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
I. Velicogna, Y. Mohajerani, G. A, F. Landerer, J. Mouginot, B. No&euml;l,
E. Rignot, T. C. Sutterley, M. van den Broeke, J. M. van Wessem, and D. Wiese,
"Continuity of ice sheet mass loss in Greenland and Antarctica from the GRACE
and GRACE Follow‐On missions", *Geophysical Research Letters*, 47,
(2020). [doi:10.1029/2020GL087291]( https://doi.org/10.1029/2020GL087291)  

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
