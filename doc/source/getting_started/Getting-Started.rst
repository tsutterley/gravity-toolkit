===============
Getting Started
===============

- `Register at NASA Earthdata and add PO.DAAC Drive OPS as application <./NASA-Earthdata.md>`_
- Run `podaac_cumulus.py <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/podaac_cumulus.py>`_ program with your NASA Earthdata credentials to acquire GRACE/GRACE-FO and auxiliary data

.. code-block:: bash

    python podaac_cumulus.py --user=<username> --directory=<path_to_grace_directory> --release RL06

- Run `cnes_grace_sync.py <https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/cnes_grace_sync.py>`_ program to acquire CNES/GRGS GRACE/GRACE-FO data

.. code-block:: bash

    python cnes_grace_sync.py --directory=<path_to_grace_directory> --release RL05

- Get geocenter results from `Sutterley and Velicogna (2019) <https://doi.org/10.3390/rs11182108>`_ if wanting to use that dataset

.. code-block:: python

    import gravity_toolkit.utilities
    gravity_toolkit.utilities.from_figshare(path_to_grace_directory,verbose=True)

- If correcting for Glacial Isostatic Adjustment: have full path to data file known

    * These can be ascii files direct from many modeling groups or a reformatted ascii/netCDF4/HDF5 file

- If correcting for other geophysical processes such as terrestrial water storage: have full path known

    * These can be a single netCDF4 or HDF5 file or an index of ascii/netCDF4/HDF5 files

- Run Jupyter notebook `GRACE-Spatial-Maps.ipynb <https://github.com/tsutterley/gravity-toolkit/blob/main/notebooks/GRACE-Spatial-Maps.ipynb>`_ to create monthly maps

    * This program uses `Jupyter widgets <https://ipywidgets.readthedocs.io/en/latest/>`_ to select `datasets <./GRACE-Data-File-Formats.html>`_ and processing parameters
    * Can also sync the data from within the Jupyter Notebook using `magics <https://ipython.readthedocs.io/en/stable/interactive/magics.html>`_
    * Can output monthly spatial maps to netCDF4 or HDF5 in specified units
    * Will create an animation of the GRACE/GRACE-FO monthly data in the specified units
