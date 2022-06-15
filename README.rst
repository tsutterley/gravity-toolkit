====================
read-GRACE-harmonics
====================

|Language|
|License|
|Documentation Status|
|Binder|
|Pangeo|
|zenodo|

.. |Language| image:: https://img.shields.io/badge/python-v3.8-green.svg
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/github/license/tsutterley/read-grace-harmonics
   :target: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/LICENSE

.. |PyPI Version| image:: https://img.shields.io/pypi/v/gravity-toolkit.svg
   :target: https://pypi.python.org/pypi/gravity-toolkit/

.. |Documentation Status| image:: https://readthedocs.org/projects/read-grace-harmonics/badge/?version=latest
   :target: https://read-grace-harmonics.readthedocs.io/en/latest/?badge=latest

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/tsutterley/read-GRACE-harmonics/main

.. |Pangeo| image:: https://img.shields.io/static/v1.svg?logo=Jupyter&label=PangeoBinderAWS&message=us-west-2&color=orange
   :target: https://aws-uswest2-binder.pangeo.io/v2/gh/tsutterley/read-GRACE-harmonics/main?urlpath=lab

.. |zenodo| image:: https://zenodo.org/badge/107323776.svg
   :target: https://zenodo.org/badge/latestdoi/107323776

Python tools for obtaining and working with Level-2 spherical harmonic coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

Resources
#########

- `NASA GRACE mission site <https://www.nasa.gov/mission_pages/Grace/index.html>`_
- `NASA GRACE-FO mission site <https://www.nasa.gov/missions/grace-fo>`_
- `JPL GRACE Tellus site <https://grace.jpl.nasa.gov/>`_
- `JPL GRACE-FO site <https://gracefo.jpl.nasa.gov/>`_
- `UTCSR GRACE site <http://www.csr.utexas.edu/grace/>`_
- `GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC) <https://podaac.jpl.nasa.gov/grace>`_
- `GRACE at the GFZ Information System and Data Center <http://isdc.gfz-potsdam.de/grace-isdc/>`_

Dependencies
############

- `numpy: Scientific Computing Tools For Python <https://www.numpy.org>`_
- `scipy: Scientific Tools for Python <https://docs.scipy.org/doc/>`_
- `dateutil: powerful extensions to datetime <https://dateutil.readthedocs.io/en/stable/>`_
- `PyYAML: YAML parser and emitter for Python <https://github.com/yaml/pyyaml>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_
- `future: Compatibility layer between Python 2 and Python 3 <https://python-future.org/>`_
- `matplotlib: Python 2D plotting library <https://matplotlib.org/>`_
- `cartopy: Python package designed for geospatial data processing <https://scitools.org.uk/cartopy/docs/latest/>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/>`_
- `h5py: Python interface for Hierarchal Data Format 5 (HDF5) <https://www.h5py.org/>`_
- `geoid-toolkit: Python utilities for calculating geoid heights from static gravity field coefficients <https://github.com/tsutterley/geoid-toolkit/>`_

References
##########

    T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self‐Consistent Ice Mass Balance
    and Regional Sea Level From Time‐Variable Gravity", *Earth and Space Science*, 7,
    (2020). `doi:10.1029/2019EA000860 <https://doi.org/10.1029/2019EA000860>`_

    T. C. Sutterley and I. Velicogna, "Improved estimates of geocenter variability
    from time-variable gravity and ocean model outputs", *Remote Sensing*, 11(18),
    2108, (2019). `doi:10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

    J. Wahr, S. C. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates",
    *Geophysical Research Letters*, 33(6), L06401, (2006).
    `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_

    J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity
    field: Hydrological and oceanic effects and their possible detection using
    GRACE", *Journal of Geophysical Research: Solid Earth*, 103(B12), (1998).
    `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

    D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified
    earth, and a further analysis of postglacial rebound", *Geophysical Journal
    International*, 120(2), (1995).
    `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

Data Repositories
#################

    T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Ice Mass and Regional Sea Level
    Estimates from Time-Variable Gravity", (2020).
    `doi:10.6084/m9.figshare.9702338 <https://doi.org/10.6084/m9.figshare.9702338>`_

    T. C. Sutterley and I. Velicogna, "Geocenter Estimates from Time-Variable
    Gravity and Ocean Model Outputs", (2019).
    `doi:10.6084/m9.figshare.7388540 <https://doi.org/10.6084/m9.figshare.7388540>`_

Download
########

| The program homepage is:
| https://github.com/tsutterley/read-GRACE-harmonics
| A zip archive of the latest version is available directly at:
| https://github.com/tsutterley/read-GRACE-harmonics/archive/main.zip

Disclaimer
##########

This project contains work and contributions from the `scientific community <./CONTRIBUTORS.rst>`_.
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.
It is provided here for your convenience but *with no guarantees whatsoever*.

License
#######

The content of this project is licensed under the `Creative Commons Attribution 4.0 Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ and the source code is licensed under the `MIT license <LICENSE>`_.
