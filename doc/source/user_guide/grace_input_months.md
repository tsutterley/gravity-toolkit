grace_input_months.py
=====================

 - Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range  
 - Replaces Degree 1 with with input values (if specified)  
 - Replaces C20 with SLR values (if specified)  
 - Replaces C30 with SLR values for months 179+ (if specified)  
 - Corrects for ECMWF atmospheric "jumps" using the GAE, GAF and GAG files following [Fagiolini et al. (2015)](https://doi.org/10.1093/gji/ggv276)  
 - Corrects for Pole Tide drift following [Wahr et al. (2015)](https://doi.org/10.1002/2015JB011986)  
 - Removes a temporal average gravity field to get geopotential anomalies  

#### Calling Sequence
```
from gravity_toolkit.grace_input_months import grace_input_months
GRACE_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
    start_mon, end_mon, missing, SLR_C20, DEG1, SLR_C30=SLR_C30)
```

#### Inputs
 - `base_dir`: Working data directory for GRACE/GRACE-FO data
 - `PROC`: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)  
 - `DREL`: GRACE/GRACE-FO data release (RL04, RL05, RL06)  
 - `DSET`: GRACE/GRACE-FO data product (GAA, GAB, GAC, GAD, GSM)  
 - `LMAX`: Upper bound of Spherical Harmonic Degrees  
 - `start_mon`: starting month to consider in analysis  
 - `end_mon`: ending month to consider in analysis  
 - `missing`: missing months to not consider in analysis  
 - `SLR_C20`: Replaces C20 with values from Satellite Laser Ranging (SLR)  
   - `None`: use original values  
   - `CSR`: use values from CSR (TN-07, TN-09, TN-11)  
   - `GSFC`: use values from GSFC (TN-14)  
 - `DEG1`: Use Degree 1 coefficients  
   - `None`: No degree 1  
   - `Tellus`: [GRACE/GRACE-FO TN-13 coefficients from PO.DAAC](https://grace.jpl.nasa.gov/data/get-data/geocenter/)  
   - `SLR`: [Satellite laser ranging coefficients from CSR](ftp://ftp.csr.utexas.edu/pub/slr/geocenter/)  
   - `SLF`: [Sutterley and Velicogna coefficients, Remote Sensing (2019)](https://doi.org/10.6084/m9.figshare.7388540)  

#### Options
 - `MMAX`: Upper bound of Spherical Harmonic Orders  
 - `SLR_C30`: Replaces C30 with values from Satellite Laser Ranging (SLR)  
    - `None`: use original values  
    - `CSR`: use values from CSR (5x5 with 6,1)  
    - `GSFC`: use values from GSFC (TN-14)  
 - `POLE_TIDE`: correct GSM data for pole tide drift  
 - `ATM`: correct data with ECMWF "jump" corrections GAE, GAF and GAG  
 - `MODEL_DEG1`: least-squares model missing degree 1 coefficients  
 - `DEG1_GIA`: GIA-correction used when calculating degree 1 coefficients  
 - `MEAN`: remove mean of harmonics  

#### Outputs
 - `clm`: GRACE/GRACE-FO cosine spherical harmonics  
 - `slm`: GRACE/GRACE-FO sine spherical harmonics  
 - `time`: time of each GRACE measurement (mid-month)  
 - `month`: GRACE/GRACE-FO months of input datasets  
 - `title`: string denoting low degree zonals replacement, geocenter usage and corrections  
 - `mean`: mean spherical harmonic fields as a dictionary with fields clm/slm  
 - `directory`: directory of exact GRACE/GRACE-FO product  
