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

#### Dependencies
- [numpy: Scientific Computing Tools For Python](https://www.numpy.org)
- [scipy: Scientific Tools for Python](https://docs.scipy.org/doc/)
- [dateutil: powerful extensions to datetime](https://dateutil.readthedocs.io/en/stable/)
- [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)
- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)
- [future: Compatibility layer between Python 2 and Python 3](https://python-future.org/)
- [matplotlib: Python 2D plotting library](https://matplotlib.org/)
- [cartopy: Python package designed for geospatial data processing](https://scitools.org.uk/cartopy/docs/latest/)
- [netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/)
- [h5py: Python interface for Hierarchal Data Format 5 (HDF5)](https://www.h5py.org/)
- [read-GRACE-geocenter: Python reader for GRACE/GRACE-FO geocenter data](https://github.com/tsutterley/read-GRACE-geocenter/)
- [geoid-toolkit: Python utilities for calculating geoid heights from static gravity field coefficients](https://github.com/tsutterley/geoid-toolkit/)

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
