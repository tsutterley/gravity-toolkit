===============
ocean_stokes.py
===============

- Reads a land-sea mask and converts to a series of spherical harmonics

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.ocean_stokes import ocean_stokes
    ocean_Ylms = ocean_stokes(LANDMASK, LMAX, MMAX=MMAX, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/ocean_stokes.py)

Arguments
#########

- ``LANDMASK``: `netCDF4 mask file <https://doi.org/10.6084/m9.figshare.9702338>`_ [Sutterley2020]_
-
    * updated 1.0, 0.5 and 0.25 degree masks from `ORNL as part of ISLSCP <https://daac.ornl.gov/ISLSCP_II/guides/combined_ancillary_xdeg.html>`_
- ``LMAX``:  maximum spherical harmonic degree of the output harmonics

Keyword arguments
#################

- ``MMAX``: maximum spherical harmonic order of the output harmonics (default: ``LMAX``)
- ``LOVE``: input load Love numbers up to degree of truncation (``hl``, ``kl``, ``ll``)
- ``VARNAME``: variable name for mask in netCDF4 file
- ``SIMPLIFY``: simplify land mask by removing isolated points

Returns
#######

- ``clm``: Cosine spherical harmonic coefficients (geodesy normalization)
- ``slm``: Sine spherical harmonic coefficients (geodesy normalization)
- ``l``: spherical harmonic degree to ``LMAX``
- ``m``: spherical harmonic order to ``MMAX``

References
##########

.. [Sutterley2020] T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self-Consistent Ice Mass Balance and Regional Sea Level From Time-Variable Gravity", *Earth and Space Science*, 7(3), (2020). `doi: 10.1029/2019EA000860 <https://doi.org/10.1029/2019EA000860>`_
