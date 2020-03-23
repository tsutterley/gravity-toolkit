GRACE Data File Formats
=======================

#### Product Identifier  
GRACE Level-2 products consist of spherical harmonic coefficients of the Earth's gravitational field.  The data files are typically gzipped ascii files with names formatted as the following: `PID-2_YYYYDOY-yyyydoy_ndays_center_flag_rrrr`    
 - `PID` is a product identification string (for standard products: GSM, GAD, GAC, GAA, GAB)   
 - `-2` denotes that the data is a GRACE Level-2 product  
 - `YYYYDOY` denotes the start date (year and day-of-year) of the measurement range  
 - `yyyydoy` denotes the end date (year and day-of-year) of the measurement range  
 - `ndays` is the number of calendar days used to produce the monthly estimate  
 - `center` is an institution specific string (UTCSR for CSR, JPLEM for JPL Spherical Harmonics, JPLMSC for JPL Mascons, EIGEN for GFZ)  
 - `flag` is a 4 character processing center dependent string (CSR denotes the maximum degree and possibly maximum order of the solutions, JPL denotes if the data is an intermediate release, GFZ denotes constrained versus unconstrained solutions)  

 - `rrrr` is a 4 character release string which is typically a 4 digit number (GFZ datasets can denote intermediate releases in the 4<sup>th</sup> character)  

#### Character Description

1st | Description   
--- | -----------    
*G* | Geopotential coefficients

2nd | Description   
--- | -----------   
*S* | Estimate made from only GRACE data
*C* | Combination estimate from GRACE and terrestrial gravity information
*E* | Any background model specified as a time-series
*A* | Average of any background model over a time period  

3rd | Description   
--- | -----------   
*M* | Estimate of the Static field <sup>1</sup>
*U* | Geopotential estimate relative to the background gravity model
*T* | Total background gravity model except for background static model
*A* | Non-tidal atmosphere <sup>2</sup>
*B* | Non-tidal Oceans <sup>2</sup>
*C* | Combination of non-tidal atmosphere and ocean <sup>2</sup>
*D* | Ocean bottom pressure product <sup>3</sup>     

[1]: Data files for this product also contain records with the epochs and rates used to model secular changes in the background gravity model.   
[2]: see [AOD1B Description Document](https://podaac-tools.jpl.nasa.gov/drive/files/allData/gracefo/docs/AOD1B_PDD_RL06_v6.1.pdf).    
[3]: Summation of non-tidal ocean and atmosphere over the ocean.  Atmospheric pressure values are equal to zero over the land.    

#### Summary

- **GSM:** Geopotential coefficients of the static gravity field estimated from GRACE satellite data (produced by all centers).  
- **GAA:** Non-tidal atmosphere geopotential coefficients averaged over certain time period (produced by GFZ & JPL).  
- **GAB:** Non-tidal ocean geopotential coefficients averaged over certai n time period (produced by GFZ & JPL).  
- **GAC:** Combination of non-tidal atmosphere and ocean averaged over certain time period (produced by all centers).  
- **GAD:** Ocean bottom pressure product.  Combination of surface pressure and ocean pressure over the oceans, zero over the land (produced by all centers).  

#### Reference  
[Dah-Ning Yuan, _Level-2 Gravity Field Product User Handbook_, Technical Document GRACE-FO D-103922, Revision 1.1, NASA Jet Propulsion Laboratory (2019).](https://podaac-tools.jpl.nasa.gov/drive/files/allData/gracefo/docs/GRACE-FO_L2-UserHandbook_v1.1.pdf)  
