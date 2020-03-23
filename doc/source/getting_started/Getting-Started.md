Getting Started
===============

1. [Register at NASA Earthdata and add `PO.DAAC Drive OPS` as application](./NASA-Earthdata.md)
2. Run `podaac_grace_sync.py` program to acquire GRACE/GRACE-FO and auxiliary data  
```bash
python podaac_grace_sync.py --user=<username> --directory=<path_to_grace_directory>
```
3. Move Load Love Numbers file from PREM into GRACE/GRACE-FO working directory  
4. Run Jupyter notebook `GRACE-Spatial-Maps.ipynb` to create monthly maps  
  - This program uses [Jupyter widgets](https://ipywidgets.readthedocs.io/en/latest/) to select [datasets](./GRACE-Data-File-Formats.md) and processing parameters  
  - Can output monthly spatial maps to netCDF4 or HDF5 in specified units
  - Will create an animation of the GRACE/GRACE-FO monthly data in the specified units  
