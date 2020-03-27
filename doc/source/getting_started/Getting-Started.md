Getting Started
===============

- [Register at NASA Earthdata and add `PO.DAAC Drive OPS` as application](./NASA-Earthdata.md)
- Run `podaac_grace_sync.py` program to acquire GRACE/GRACE-FO and auxiliary data  
```bash
python podaac_grace_sync.py --user=<username> --directory=<path_to_grace_directory>
```
- Move Load Love Numbers file from PREM into GRACE/GRACE-FO working directory  
- If correcting for Glacial Isostatic Adjustment: copy full path to data file  
    * These files can ascii files direct from many modeling groups or a reformatted netCDF4/HDF5 file  
- Run Jupyter notebook `GRACE-Spatial-Maps.ipynb` to create monthly maps  
    * This program uses [Jupyter widgets](https://ipywidgets.readthedocs.io/en/latest/) to select [datasets](./GRACE-Data-File-Formats.md) and processing parameters  
    * Can output monthly spatial maps to netCDF4 or HDF5 in specified units
    * Will create an animation of the GRACE/GRACE-FO monthly data in the specified units  
