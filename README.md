read-GRACE-harmonics
====================

Reads Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions  

- [NASA GRACE mission site](http://www.nasa.gov/mission_pages/Grace/index.html)  
- [JPL GRACE Tellus site](http://grace.jpl.nasa.gov/)  
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)  
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)  
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)  

#### Calling Sequence
```
from read_GRACE_harmonics import read_GRACE_harmonics
CSR_L2_input = read_GRACE_harmonics('GSM-2_2002095-2002120_0021_UTCSR_0060_0005.gz',60)
GFZ_L2_input = read_GRACE_harmonics('GSM-2_2002094-2002120_0024_EIGEN_G---_005a.gz',90)
JPL_L2_input = read_GRACE_harmonics('GSM-2_2002091-2002120_0018_JPLEM_0001_0005.gz',60)
JPLMSC_input = read_GRACE_harmonics('GSM-2_2003001-2003031_0029_JPLMSC_0719_0005',719)
```

#### Inputs
 1. full path to input GRACE file  
 2. spherical harmonic degree of truncation (`LMAX`)  

#### Options
 - `MMAX`: spherical harmonic order of truncation (default is `LMAX`)  
 - `POLE_TIDE`: correct GSM data for pole tide drift following [Wahr et al. (2015)](https://doi.org/10.1002/2015JB011986)  

#### Outputs
 - `time`: mid-month date of GRACE file in year-decimal  
 - `start`: start date of range as Julian day  
 - `end`: end date of range as Julian day  
 - `clm`: cosine spherical harmonics of input data  
 - `slm`: sine spherical harmonics of input data  
 - `eclm`: cosine spherical harmonic uncalibrated standard deviations  
 - `eslm`: sine spherical harmonic uncalibrated standard deviations  
 - `header`: text header of the GRACE file (will parse new YAML headers)  

#### Dependencies
 - [numpy: Scientific Computing Tools For Python](http://www.numpy.org)  
 - [PyYAML: YAML parser and emitter for Python](https://github.com/yaml/pyyaml)  

#### Download
The program homepage is:   
https://github.com/tsutterley/read-GRACE-harmonics   
A zip archive of the latest version is available directly at:    
https://github.com/tsutterley/read-GRACE-harmonics/archive/master.zip  

#### Disclaimer  
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL), the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.  
