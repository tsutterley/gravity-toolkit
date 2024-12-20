============
ocean_stokes
============

- Reads a land-sea mask and converts to a series of spherical harmonics
- `netCDF4 land-sea mask files <https://doi.org/10.6084/m9.figshare.9702338>`_ from :cite:p:`Sutterley:2020js`
    * updated 1.0, 0.5 and 0.25 degree masks from `ORNL as part of ISLSCP <https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html>`_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.ocean_stokes import ocean_stokes
    ocean_Ylms = ocean_stokes(LANDMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/ocean_stokes.py)

.. autofunction:: gravity_toolkit.ocean_stokes

.. autofunction:: gravity_toolkit.ocean_stokes.find_isolated_points
