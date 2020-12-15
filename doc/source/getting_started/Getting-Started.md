Getting Started
===============

- [Register at NASA Earthdata and add `PO.DAAC Drive OPS` as application](./NASA-Earthdata.md)
- Get PO.DAAC WebDAV credentials by running [`podaac_webdav.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_webdav.py) or logging onto [PO.DAAC Drive](https://podaac-tools.jpl.nasa.gov/drive)  
- Run [`podaac_grace_sync.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/podaac_grace_sync.py) program with your WebDAV credentials to acquire GRACE/GRACE-FO and auxiliary data  
```bash
python podaac_grace_sync.py --user=<username> --directory=<path_to_grace_directory> --release RL06
```
- Run [`cnes_grace_sync.py`](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/cnes_grace_sync.py) program to acquire CNES/GRGS GRACE/GRACE-FO data  
```bash
python cnes_grace_sync.py --directory=<path_to_grace_directory> --release RL05
```
- Get geocenter results from [Sutterley and Velicogna (2019)](https://doi.org/10.3390/rs11182108) if wanting to use that dataset
```python
>>> import gravity_toolkit.utilities
>>> geocenter_directory = os.path.join(grace_dir,'geocenter')
>>> gravity_toolkit.utilities.from_figshare(geocenter_directory,verbose=True)
```
- If correcting for Glacial Isostatic Adjustment: have full path to data file known  
    * These can be ascii files direct from many modeling groups or a reformatted ascii/netCDF4/HDF5 file  
- If correcting for other geophysical processes such as terrestrial water storage: have full path known  
    * These can be a single netCDF4 or HDF5 file or an index of ascii/netCDF4/HDF5 files  
- Run Jupyter notebook `GRACE-Spatial-Maps.ipynb` to create monthly maps  
    * This program uses [Jupyter widgets](https://ipywidgets.readthedocs.io/en/latest/) to select [datasets](./GRACE-Data-File-Formats.md) and processing parameters  
    * Can also sync the data from within the Jupyter Notebook using [magics](https://ipython.readthedocs.io/en/stable/interactive/magics.html)
    * Can output monthly spatial maps to netCDF4 or HDF5 in specified units
    * Will create an animation of the GRACE/GRACE-FO monthly data in the specified units  
