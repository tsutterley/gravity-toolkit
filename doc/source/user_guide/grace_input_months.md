grace_input_months.py
=====================

 - Reads GRACE/GRACE-FO files for a specified spherical harmonic degree and order and for a specified date range  
 - Replaces Degree 1 with with input values (if specified)  
 - Replaces C20 with SLR values (if specified)  
 - Replaces C30 with SLR values for months 179+ (if specified)  
 - Corrects for ECMWF atmospheric "jumps" using the GAE, GAF and GAG files following [Fagiolini et al. (2015)](https://doi.org/10.1093/gji/ggv276)  
 - Corrects for Pole Tide drift following [Wahr et al. (2015)](https://doi.org/10.1002/2015JB011986)  

#### Calling Sequence
```python
from gravity_toolkit.grace_input_months import grace_input_months
GRACE_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
    start_mon, end_mon, missing, SLR_C20, DEG1, SLR_C30=SLR_C30)
```

#### Inputs
 1. `base_dir`: Working data directory for GRACE/GRACE-FO data
 2. `PROC`: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)  
    * `'CSR'`: University of Texas Center for Space Research  
    * `'GFZ'`: German Research Centre for Geosciences (GeoForschungsZentrum)
    * `'JPL'`: Jet Propulsion Laboratory    
    * `'CNES'`: French Centre National D'Etudes Spatiales
 3. `DREL`: GRACE/GRACE-FO data release (RL04, RL05, RL06)  
 4. `DSET`: GRACE/GRACE-FO data product (GAA, GAB, GAC, GAD, GSM)  
    * `'GAA'`: non-tidal atmospheric correction  
    * `'GAB'`: non-tidal oceanic correction  
    * `'GAC'`: combined non-tidal atmospheric and oceanic correction  
    * `'GAD'`: GRACE/GRACE-FO ocean bottom pressure product  
    * `'GSM'`: corrected monthly GRACE/GRACE-FO static field product
 5. `LMAX`: Upper bound of Spherical Harmonic Degrees  
 6. `start_mon`: starting month to consider in analysis  
 7. `end_mon`: ending month to consider in analysis  
 8. `missing`: missing months to not consider in analysis  
 9. `SLR_C20`: Replaces C20 with values from Satellite Laser Ranging (SLR)  
    * `None`: use original values  
    * `'CSR'`: use values from CSR (TN-07, TN-09, TN-11)  
    * `'GSFC'`: use values from GSFC (TN-14)  
 10. `DEG1`: Use Degree 1 coefficients  
    * `None`: No degree 1  
    * `'Tellus'`: [GRACE/GRACE-FO TN-13 coefficients from PO.DAAC](https://grace.jpl.nasa.gov/data/get-data/geocenter/)  
    * `'SLR'`: [Satellite laser ranging coefficients from CSR](ftp://ftp.csr.utexas.edu/pub/slr/geocenter/)  
    * `'SLF'`: [Sutterley and Velicogna coefficients, Remote Sensing (2019)](https://doi.org/10.6084/m9.figshare.7388540)  

#### Options
 - `MMAX`: Upper bound of Spherical Harmonic Orders  
 - `SLR_C30`: Replaces C30 with values from Satellite Laser Ranging (SLR)  
    * `None`: use original values  
    * `'CSR'`: use values from CSR (5x5 with 6,1)  
    * `'GSFC'`: use values from GSFC (TN-14)  
 - `POLE_TIDE`: correct GSM data for pole tide drift  
 - `ATM`: correct data with ECMWF "jump" corrections GAE, GAF and GAG  
 - `MODEL_DEG1`: least-squares model missing degree 1 coefficients  
 - `DEG1_GIA`: GIA-correction used when calculating degree 1 coefficients  

#### Outputs
 - `clm`: GRACE/GRACE-FO cosine spherical harmonics to degree/order LMAX and MMAX  
 - `slm`: GRACE/GRACE-FO sine spherical harmonics to degree/order LMAX and MMAX  
 - `time`: time of each GRACE measurement (mid-month)  
 - `month`: GRACE/GRACE-FO months of input datasets  
 - `l`: spherical harmonic degree to LMAX
 - `m`: spherical harmonic order to MMAX
 - `title`: string denoting low degree zonals replacement, geocenter usage and corrections  
 - `directory`: directory of exact GRACE/GRACE-FO product  
