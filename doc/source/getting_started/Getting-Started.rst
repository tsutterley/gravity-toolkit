===============
Getting Started
===============



This documentation is intended to explain how to compute spatial and time series
estimates using GRACE/GRACE-FO time-variable gravity measurements.
``gravity-toolkit`` is a Python-based geophysical software that reads
GRACE/GRACE-FO time-variable gravity solutions for estimating regional mass change.
A suite of geophysical corrections can be applied to the gravity solutions to
optimize the GRACE/GRACE-FO data for particular applications.
This software was developed with the goal of supporting science applications for
time-variable gravity.
``gravity-toolkit`` provides data access utilities for ascii, netCDF4, HDF5 and gfc file formats.
``gravity-toolkit`` also provides some very high-level plotting programs through the
use of `Jupyter Notebooks <../user_guide/Examples.html>`_.

.. graphviz::
    :caption: Data Processing Framework
    :align: center

    digraph {
        G [label="GRACE/GRACE-FO\ntime-variable gravity"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        A [label="Non-tidal Ocean and\nAtmospheric Variation"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        I [label="Glacial Isostatic\nAdjustment"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        W [label="Terrestrial Water\nStorage"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        R [label="gravity-toolkit"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        S [label="Spatial Maps"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Spatial-Maps.html"]
        T [label="Time Series Analysis"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Time-Series-Analysis.html"]
        D [label="Geocenter Variation"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"
            URL="Geocenter-Variations.html"]
        G -> R [arrowsize=0.8]
        A -> R [arrowsize=0.8]
        I -> R [arrowsize=0.8]
        W -> R [arrowsize=0.8]
        R -> S [arrowsize=0.8]
        R -> T [arrowsize=0.8]
        R -> D [arrowsize=0.8]
    }

Steps to Get Started
####################

- `Register at NASA Earthdata <./NASA-Earthdata.html>`_
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

- Run Jupyter notebook `Examples <../user_guide/Examples.html>`_ 

    * These programs use `Jupyter widgets <https://ipywidgets.readthedocs.io/en/latest/>`_ to select `datasets <./GRACE-Data-File-Formats.html>`_ and processing parameters
    * Can also sync the data from within the Jupyter Notebook using `magics <https://ipython.readthedocs.io/en/stable/interactive/magics.html>`_
