========
tools.py
========

`User interface <https://ipywidgets.readthedocs.io/en/latest/>`_ and plotting tools for use in `Jupyter notebooks <https://jupyter.org/>`_

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/tools.py

General Attributes and Methods
==============================

.. class:: widgets(object)

    .. attribute:: object.directory

        Text widget for setting working data directory

    .. attribute:: object.directory_button

        Button widget for setting working data directory with `Tkinter file dialog <https://docs.python.org/3/library/dialog.html>`_

    .. attribute:: object.update

        Checkbox widget for updating GRACE/GRACE-FO data in directory

    .. attribute:: object.center

        Dropdown menu widget for setting processing center

    .. attribute:: object.release

        Dropdown menu widget for setting GRACE/GRACE-FO data release

    .. attribute:: object.product

        Dropdown menu widget for setting GRACE/GRACE-FO data product

    .. attribute:: object.months

        Selection widget for setting GRACE/GRACE-FO months

    .. attribute:: object.lmax

        Text entry widget for setting spherical harmonic degree

    .. attribute:: object.mmax

        Text entry widget for setting spherical harmonic order

    .. attribute:: object.geocenter

        Dropdown menu widget for setting geocenter data product

    .. attribute:: object.C20

        Dropdown menu widget for setting *C*\ :sub:`20` data product

    .. attribute:: object.CS21

        Dropdown menu widget for setting *C*\ :sub:`21` and *S*\ :sub:`21` data product

    .. attribute:: object.CS22

        Dropdown menu widget for setting *C*\ :sub:`22` and *S*\ :sub:`22` data product

    .. attribute:: object.C30

        Dropdown menu widget for setting *C*\ :sub:`30` data product

    .. attribute:: object.C50

        Dropdown menu widget for setting *C*\ :sub:`50` data product

    .. attribute:: object.pole_tide

        Checkbox widget for correcting for Pole Tide Drift [Wahr2015]_

    .. attribute:: object.atm

        Checkbox widget for correcting ECMWF Atmospheric Jumps [Fagiolini2015]_

    .. attribute:: object.GIA_file

        Text entry widget for setting GIA correction file

    .. attribute:: object.GIA_button

        Button widget for setting GIA correction file with `Tkinter file dialog <https://docs.python.org/3/library/dialog.html>`_

    .. attribute:: object.GIA

        Dropdown menu for setting GIA model file type

    .. attribute:: object.remove_file

        Text entry widget for setting spherical harmonic files to be removed

    .. attribute:: object.remove_button

        Button widget for setting remove files with `Tkinter file dialog <https://docs.python.org/3/library/dialog.html>`_

    .. attribute:: object.remove_format

        Dropdown menu for setting remove file type

    .. attribute:: object.redistribute_removed

        Checkbox widget for redestributing removed file mass over the ocean

    .. attribute:: object.mask

        Text entry widget for setting land-sea mask file for ocean redistribution

    .. attribute:: object.mask_button

        Button widget for setting land-sea mask files with `Tkinter file dialog <https://docs.python.org/3/library/dialog.html>`_

    .. attribute:: object.gaussian

        Text entry widget for setting Gaussian Smoothing Radius in kilometers

    .. attribute:: object.destripe

        Checkbox widget for destriping spherical harmonics [Swenson2006]_

    .. attribute:: object.spacing

        Text entry widget for setting output spatial degree spacing

    .. attribute:: object.interval

        Dropdown menu widget for setting output degree interval

    .. attribute:: object.units

        Dropdown menu widget for setting output units

    .. attribute:: object.output_format

        Dropdown menu widget for setting output file format


.. class:: colormap(object)

    .. attribute:: object.range

        Slider widget for setting output colormap normalization

    .. attribute:: object.step

        Slider widget for setting output colormap discretization

    .. attribute:: object.name

        Dropdown widget for setting output `colormap <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_

    .. attribute:: object.reverse

        Checkbox widget for reversing the output colormap

.. function:: from_cpt(filename, use_extremes=True)

    Reads GMT color palette table files and registers the colormap to be recognizable by ``plt.cm.get_cmap()``

.. function:: custom_colormap(N, map_name)

    Calculates a custom colormap and registers it to be recognizable by ``plt.cm.get_cmap()``

    Arguments:

        - ``N``: number of slices in initial HSV color map

        - ``map_name``: name of color map

            * ``'Joughin'``: [Joughin2018]_ standard velocity colormap

            * ``'Rignot'``: [Rignot2011]_ standard velocity colormap

            * ``'Seroussi'``:  [Seroussi2011]_ velocity divergence colormap


References
##########

.. [Fagiolini2015] E. Fagiolini, F. Flechtner, M. Horwath, and H. Dobslaw, "Correction of inconsistencies in ECMWF's operational analysis data during de-aliasing of GRACE gravity models", *Geophysical Journal International*, 202(3), 2150--2158, (2015). `doi: 10.1093/gji/ggv276 <https://doi.org/10.1093/gji/ggv276>`_

.. [Joughin2018] I. Joughin, B. E. Smith, and I. Howat, "Greenland Ice Mapping Project: ice flow velocity variation at sub-monthly to decadal timescales", *The Cryosphere*, 12, 2211--2227, (2018). `doi: 10.5194/tc-12-2211-2018 <https://doi.org/10.5194/tc-12-2211-2018>`_

.. [Rignot2011] E. Rignot J. Mouginot, and B. Scheuchl, "Ice Flow of the Antarctic Ice Sheet", *Science*, 333(6048), 1427--1430, (2011). `doi: 10.1126/science.1208336 <https://doi.org/10.1126/science.1208336>`_

.. [Seroussi2011] H. Seroussi, M. Morlighem, E. Rignot, E. Larour, D. Aubry, H. Ben Dhia, and S. S. Kristensen, "Ice flux divergence anomalies on 79north Glacier, Greenland", *Geophysical Research Letters*, 38(L09501), (2011). `doi: 10.1029/2011GL047338 <https://doi.org/10.1029/2011GL047338>`_

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Wahr2015] J. Wahr, R. S. Nerem, and S. V. Bettadpur, "The pole tide and its effect on GRACE time‐variable gravity measurements: Implications for estimates of surface mass variations". *Journal of Geophysical Research: Solid Earth*, 120(6), 4597--4615, (2015). `doi: 10.1002/2015JB011986 <https://doi.org/10.1002/2015JB011986>`_

