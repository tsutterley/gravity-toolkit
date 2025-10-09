=============================
gravity-toolkit Documentation
=============================

Welcome to the documentation for ``gravity-toolkit``, a set of Python tools for working with time-variable gravity fields.

This documentation is intended to explain how to obtain and work with the Level-2 spherical harmonic
coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and
the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions.

Introduction
------------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Installation
      :text-align: center
      :link: ./getting_started/Install.html

      :material-outlined:`download;5em`

    .. grid-item-card::  Getting Started
      :text-align: center
      :link: ./getting_started/Getting-Started.html

      :material-outlined:`hiking;5em`

    .. grid-item-card::  Background
      :text-align: center
      :link: ./background/Background.html

      :material-outlined:`library_books;5em`

    .. grid-item-card::  Examples
      :text-align: center
      :link: ./user_guide/Examples.html

      :material-outlined:`apps;5em`

Contribute
----------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Guidelines
      :text-align: center
      :link: ./getting_started/Contributing.html

      :material-outlined:`groups;5em`

    .. grid-item-card::  Code of Conduct
      :text-align: center
      :link: ./getting_started/Code-of-Conduct.html

      :material-outlined:`gavel;5em`

    .. grid-item-card::  Discussions
      :text-align: center
      :link: https://github.com/tsutterley/gravity-toolkit/discussions

      :material-outlined:`forum;5em`

    .. grid-item-card::  Citation Information
      :text-align: center
      :link: ./project/Citations.html

      :material-outlined:`alternate_email;5em`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    getting_started/Install.rst
    getting_started/Getting-Started.rst
    getting_started/NASA-Earthdata.rst
    getting_started/GRACE-Data-File-Formats.rst
    getting_started/Contributing.rst
    getting_started/Code-of-Conduct.rst
    getting_started/Resources.rst

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Background

    background/Background.rst
    background/Spatial-Maps.rst
    background/Time-Series-Analysis.rst
    background/Geocenter-Variations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: User Guide

    user_guide/Examples.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: API Reference

    api_reference/associated_legendre.rst
    api_reference/clenshaw_summation.rst
    api_reference/degree_amplitude.rst
    api_reference/destripe_harmonics.rst
    api_reference/fourier_legendre.rst
    api_reference/gauss_weights.rst
    api_reference/gen_averaging_kernel.rst
    api_reference/gen_disc_load.rst
    api_reference/gen_harmonics.rst
    api_reference/gen_point_load.rst
    api_reference/gen_spherical_cap.rst
    api_reference/gen_stokes.rst
    api_reference/geocenter.rst
    api_reference/grace_date.rst
    api_reference/grace_find_months.rst
    api_reference/grace_input_months.rst
    api_reference/grace_months_index.rst
    api_reference/harmonic_gradients.rst
    api_reference/harmonic_summation.rst
    api_reference/harmonics.rst
    api_reference/legendre.rst
    api_reference/legendre_polynomials.rst
    api_reference/mascons.rst
    api_reference/ocean_stokes.rst
    api_reference/read_gfc_harmonics.rst
    api_reference/read_GIA_model.rst
    api_reference/read_GRACE_harmonics.rst
    api_reference/read_love_numbers.rst
    api_reference/read_SLR_harmonics.rst
    api_reference/sea_level_equation.rst
    api_reference/SLR/C20.rst
    api_reference/SLR/CS2.rst
    api_reference/SLR/C30.rst
    api_reference/SLR/C40.rst
    api_reference/SLR/C50.rst
    api_reference/spatial.rst
    api_reference/time.rst
    api_reference/time_series/amplitude.rst
    api_reference/time_series/fit.rst
    api_reference/time_series/lomb_scargle.rst
    api_reference/time_series/piecewise.rst
    api_reference/time_series/regress.rst
    api_reference/time_series/savitzky_golay.rst
    api_reference/time_series/smooth.rst
    api_reference/tools.rst
    api_reference/units.rst
    api_reference/utilities.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Access

    api_reference/access/cnes_grace_sync.rst
    api_reference/access/esa_costg_swarm_sync.rst
    api_reference/access/gfz_icgem_costg_ftp.rst
    api_reference/access/gfz_isdc_dealiasing_sync.rst
    api_reference/access/gfz_isdc_grace_sync.rst
    api_reference/access/itsg_graz_grace_sync.rst
    api_reference/access/podaac_cumulus.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Dealiasing

    api_reference/dealiasing/aod1b_geocenter.rst
    api_reference/dealiasing/aod1b_oblateness.rst
    api_reference/dealiasing/dealiasing_global_uplift.rst
    api_reference/dealiasing/dealiasing_monthly_mean.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Geocenter

    api_reference/geocenter/calc_degree_one.rst
    api_reference/geocenter/monte_carlo_degree_one.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Use Cases

    api_reference/scripts/calc_mascon.rst
    api_reference/scripts/calc_harmonic_resolution.rst
    api_reference/scripts/calc_sensitivity_kernel.rst
    api_reference/scripts/combine_harmonics.rst
    api_reference/scripts/convert_harmonics.rst
    api_reference/scripts/grace_mean_harmonics.rst
    api_reference/scripts/grace_raster_grids.rst
    api_reference/scripts/grace_spatial_error.rst
    api_reference/scripts/grace_spatial_maps.rst
    api_reference/scripts/mascon_reconstruct.rst
    api_reference/scripts/piecewise_grace_maps.rst
    api_reference/scripts/regress_grace_maps.rst
    api_reference/scripts/run_sea_level_equation.rst
    api_reference/scripts/scale_grace_maps.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Mapping

    api_reference/mapping/plot_AIS_grid_maps.rst
    api_reference/mapping/plot_AIS_grid_3maps.rst
    api_reference/mapping/plot_AIS_grid_4maps.rst
    api_reference/mapping/plot_AIS_grid_movie.rst
    api_reference/mapping/plot_AIS_GrIS_maps.rst
    api_reference/mapping/plot_AIS_regional_maps.rst
    api_reference/mapping/plot_AIS_regional_movie.rst
    api_reference/mapping/plot_global_grid_maps.rst
    api_reference/mapping/plot_global_grid_3maps.rst
    api_reference/mapping/plot_global_grid_4maps.rst
    api_reference/mapping/plot_global_grid_5maps.rst
    api_reference/mapping/plot_global_grid_9maps.rst
    api_reference/mapping/plot_global_grid_movie.rst
    api_reference/mapping/plot_GrIS_grid_maps.rst
    api_reference/mapping/plot_GrIS_grid_3maps.rst
    api_reference/mapping/plot_GrIS_grid_5maps.rst
    api_reference/mapping/plot_GrIS_grid_movie.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Utilities

    api_reference/utilities/make_grace_index.rst
    api_reference/utilities/quick_mascon_plot.rst
    api_reference/utilities/quick_mascon_regress.rst
    api_reference/utilities/run_grace_date.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Project Details

    project/Contributors.rst
    project/Licenses.rst
    project/Testing.rst
    project/Citations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Bibliography

    project/Bibliography.rst
