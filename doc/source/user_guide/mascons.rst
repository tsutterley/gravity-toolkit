==========
mascons.py
==========

Conversion routines for publicly available GRACE/GRACE-FO mascon solutions

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/mascons.py


General Methods
===============

.. method:: gravity_toolkit.mascons.to_gsfc(gdata, lon, lat, lon_center, lat_center, lon_span, lat_span)

    Converts an input gridded field to an output GSFC mascon array [Luthcke2013]_

    Arguments:

        ``gdata``: (lat x lon) array of gridded map

        ``lon``: column vector of defined longitude points

        ``lat``: column vector of defined latitude points

        ``lon_center``: mascon longitudinal center points

        ``lat_center``: mascon latitudinal center points

        ``lon_span``: mascon longitudinal central angles

        ``lat_span``: mascon latitudinal central angles

    Returns:

        ``data``: row vector of mascons

        ``lat_center``: row vector of latitude values for mascon centers

        ``lon_center``: row vector of longitude values for mascon centers


.. method:: gravity_toolkit.mascons.to_jpl(gdata, lon, lat, lon_bound, lat_bound)

    Converts an input gridded field to an output JPL mascon array [Watkins2015]_

    Arguments:

        ``gdata``: (lat x lon) array of gridded map

        ``lon``: column vector of defined longitude points

        ``lat``: column vector of defined latitude points

        ``lon_bound``: mascon longitudinal bounds from coordinate file

        ``lat_bound``: mascon latitudinal bounds from coordinate file

    Returns:

        ``data``: row vector of mascons

        ``mask``: row vector of mask values showing if mascon has no data

        ``lat``: row vector of latitude values for mascons

        ``lon``: row vector of longitude values for mascons


.. method:: gravity_toolkit.mascons.from_gsfc(mscdata, grid_spacing, lon_center, lat_center, lon_span, lat_span, TRANSPOSE=False)

    Converts an input GSFC mascon array to an output gridded field [Luthcke2013]_

    Arguments:

        ``mscdata``: row vector of mascons

        ``grid_spacing``: spacing of the lat/lon grid

        ``lon_center``: mascon longitudinal center points

        ``lat_center``: mascon latitudinal center points

        ``lon_span``: mascon longitudinal central angles

        ``lat_span``: mascon latitudinal central angles

    Keyword arguments:

        ``TRANSPOSE``: transpose output matrix (lon/lat)

    Returns:

        ``mdata``: distributed mass grid

.. method:: gravity_toolkit.mascons.from_jpl(mscdata, grid_spacing, lon_bound, lat_bound, TRANSPOSE=False)

    Converts an input JPL mascon array to an output gridded field [Watkins2015]_

    Arguments:

        ``mscdata``: row vector of mascons

        ``grid_spacing``: spacing of lat/lon grid

        ``lon_bound``: mascon longitudinal bounds from coordinate file

        ``lat_bound``: mascon latitudinal bounds from coordinate file

    Keyword arguments:

        ``TRANSPOSE``: transpose output matrix (lon/lat)

    Returns

        ``mdata``: distributed mass grid

References
##########

.. [Luthcke2013] S. B. Luthcke, T. J. Sabaka, B. D. Loomis, A. A. Arendt, J. J. McCarthy, and J. Camp, "Antarctica, Greenland and Gulf of Alaska land-ice evolution from an iterated GRACE global mascon solution", *Journal of Glaciology*, 59(216), (2013). `doi: 10.3189/2013JoG12J147 <https://doi.org/10.3189/2013JoG12J147>`_

.. [Watkins2015] M. M. Watkins, D. N. Wiese, D.-N. Yuan, C. Boening, and F. W. Landerer, "Improved methods for observing Earth's time variable mass distribution with GRACE using spherical cap mascons". *Journal of Geophysical Research: Solid Earth*, 120(4), 2648--2671, (2015). `doi: 10.1002/2014JB011547 <https://doi.org/10.1002/2014JB011547>`_
