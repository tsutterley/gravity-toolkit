Time Series Analysis
====================

Least-squares mascons are a method of extracting a regional signal from the GRACE/GRACE-FO spherical harmonic data.  The procedure was outlined in procedure outlined in [Tiwari et al. (2009)](https://doi.org/10.1029/2009GL039401) and [Jacob et al. (2012)](https://doi.org/10.1038/nature10847).  Least-squares mascons can be considered a post-processing technique for analyzing the GRACE/GRACE-FO data. The technique calculates the scaling factor between an input kernel and the GRACE/GRACE-FO data for a given month. For example, for a uniform kernel equivalent to 1 cm w.e., if GRACE/GRACE-FO measures 6 cm w.e. in the region, then the scaling factor would be 6.  

Ideally we would want as small of kernels as possible to get as geophysically relevant as possible (perfect fit to region shape).  However, GRACE/GRACE-FO data is limited to a specific spherical harmonic range with noise increasing with higher degree and order (data-to-noise ratio decreases).  If the kernels are too small, there will be a higher degree of ringing when truncated (Gibbs phenomenon) and require more information at the higher degrees and orders to distinguish from the adjacent kernels.  Thus there is a balance between wanting kernels as geophysically relevant as possible with wanting idealistic kernels which minimize ringing and are resolvable. Getting the kernels "just right" in order to isolate regions of interest takes some time.  

The set of least-squares mascon programs have been used in [Velicogna et al. (2014)](https://doi.org/10.1002/2014GL061052) and other publications for regional time series analysis. The [`calc_mascon.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/calc_mascon.py) program additionally calculates the GRACE/GRACE-FO error harmonics following [Wahr et al. (2006)](https://doi.org/10.1029/2005GL025305).

The [`calc_mascon.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/calc_mascon.py) program will output a text file of the time series for each mascon (format: GRACE/GRACE-FO month, mid-month date in decimal-year format, estimated monthly mass anomaly [Gt], estimated monthly error [Gt], mascon area [km<sup>2</sup>]).

#### Parameter Files

The [`calc_mascon.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/calc_mascon.py) program works by accepting parameter file system arguments (`sys.argv`) listed after the program. These parameter files can be ran individually or in a series. 
```bash
python calc_mascon.py inp1 inp2 inp3
python calc_mascon.py --np=2 inp1 inp2 inp3
python calc_mascon.py -P 2 inp1 inp2 inp3
```
The program can be run in series (default) or in parallel with the python multiprocessing package using options `--np` and `-P`.  The number after the `--np` and `-P` options declares the number of processes to run in parallel (here 2).  

The parameter files are read line by line to fill a python dictionary variable mapping specific parameter names with their values.  The parameter names are the first column in the file and the parameter values are the second column.  For parameters consisting of lists (e.g. GRACE/GRACE-FO missing months), the parameter values are separated by commas.

- Column 1: parameter name (such as `LMAX`)  
- Column 2: parameter value (e.g. `60`)  
- Column 3: comments (which are discarded)

#### Dataset Parameters

- `PROC`: GRACE Processing Center (CSR, GFZ, JPL, CNES)
- `DREL`: GRACE data release for given processing center
- `DSET`: GRACE data product (see [GRACE Data File Formats](./GRACE-Data-File-Formats.md))
- `LMIN`: minimum spherical harmonic degree (lower bound of truncation)
- `LMAX`: maximum spherical harmonic degree (upper bound of truncation)
- `MMAX`: maximum spherical harmonic order (None if LMAX)
- `START`: first month to be analyzed
- `END`: last month to be analyzed
- `MISSING`: GRACE/GRACE-FO months that are not be analyzed (see available GRACE/GRACE-FO months)
- `RAD`: Gaussian smoothing radius in km ([Jekeli, 1981](http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf))  
- `DESTRIPE`: filter coefficients using destriping procedure ([Swenson et al., 2006](https://doi.org/10.1029/2005GL025285))
- `SLR_C20`: replace C<sub>20</sub> coefficients with values from Satellite Laser Ranging (SLR)  
	* None: use original values
	* CSR: use values from CSR (TN-07,TN-09,TN-11)
	* GSFC: use values from GSFC (TN-14)
- `SLR_C30`: replace C<sub>30</sub> coefficients with values from Satellite Laser Ranging (SLR)  
	* None: use original values
	* CSR: use values from CSR (5x5 with 6,1)
	* GSFC: use values from GSFC (TN-14)
	* LARES: use filtered values from CSR (John Ries)
- `DEG1`: account for variations in geocenter with specified values
	* None
	* Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
	* SLR: satellite laser ranging coefficients from CSR  
	* SLF: Sutterley and Velicogna coefficients, Remote Sensing (2019)
- `GIA`: GIA model type
	* `None`
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
- `GIA_FILE`: path to specific GIA file to be read
- `DATAFORM`: input data format for mascon files and files to be removed from the GRACE/GRACE-FO data
	* `'ascii'`
	* `'netCDF4'`
	* `'HDF5'`
- `DIRECTORY`: Directory to output data (will create directory if non-existent)
- `MASCON_INDEX`: file index listing the full path to each mascon file to fit to the GRACE data
- `FIT_METHOD`: method of fitting mascons coefficients
     * 1: convert coefficients to mass
     * 2: keep coefficients as normalized geoid
- `MEAN`: Remove a mean field to isolate the time-variable gravity field
- `MEAN_FILE`: use a file to remove as static field (default: mean of imported month)
- `MEANFORM`: Data format for input `MEAN_FILE`  
	* `'ascii'`
	* `'netCDF4'`
	* `'HDF5'`
     * `'gfc'`
- `REMOVE_INDEX`: Remove sets of spherical harmonics using a file index (can be multiple indices)
- `REDISTRIBUTE_REMOVED`: Redistribute total mass of removed harmonics over the ocean
- `MASCON_OCEAN`: remove uniformly distributed mascon mass over ocean
- `RECONSTRUCT`: remove the reconstructed time series for a region to get the statistical leakage
- `POLE_TIDE`: correct GSM C<sub>21</sub> and S<sub>21</sub> for pole tide ([Wahr et al., 2015](https://doi.org/10.1002/2015JB011986))
- `ATM`: correct Atmosphere with [ECMWF "jump" corrections](https://doi.org/10.1093/gji/ggv276)
