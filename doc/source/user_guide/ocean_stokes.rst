===============
ocean_stokes.py
===============

- Reads a land-sea mask and converts to a series of spherical harmonics
- `netCDF4 land-sea mask files <https://doi.org/10.6084/m9.figshare.9702338>`_ from [Sutterley2020]_
    * updated 1.0, 0.5 and 0.25 degree masks from `ORNL as part of ISLSCP <https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html>`_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.ocean_stokes import ocean_stokes
    ocean_Ylms = ocean_stokes(LANDMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ocean_stokes.py)

.. autofunction:: gravity_toolkit.ocean_stokes

.. [Sutterley2020] T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self-Consistent Ice Mass Balance and Regional Sea Level From Time-Variable Gravity", *Earth and Space Science*, 7(3), (2020). `doi: 10.1029/2019EA000860 <https://doi.org/10.1029/2019EA000860>`_
