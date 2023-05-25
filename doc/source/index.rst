===============
gravity-toolkit
===============

Python tools for obtaining and working with Level-2 spherical harmonic
coefficients from the NASA/DLR Gravity Recovery and Climate Experiment (GRACE)
and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO)
missions

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    getting_started/Install.rst
    getting_started/Background.rst
    getting_started/Getting-Started.rst
    getting_started/NASA-Earthdata.rst
    getting_started/GRACE-Data-File-Formats.rst
    getting_started/Spatial-Maps.rst
    getting_started/Time-Series-Analysis.rst
    getting_started/Geocenter-Variations.rst
    getting_started/Contributing.rst
    getting_started/Resources.rst
    getting_started/Citations.rst

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
    :caption: Utilities

    api_reference/cnes_grace_sync.rst
    api_reference/esa_costg_swarm_sync.rst
    api_reference/gfz_icgem_costg_ftp.rst
    api_reference/gfz_isdc_dealiasing_ftp.rst
    api_reference/gfz_isdc_grace_ftp.rst
    api_reference/itsg_graz_grace_sync.rst
    api_reference/make_grace_index.rst
    api_reference/podaac_cumulus.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Use Cases

    api_reference/aod1b_geocenter.rst
    api_reference/aod1b_oblateness.rst
    api_reference/calc_degree_one.rst
    api_reference/calc_mascon.rst
    api_reference/calc_harmonic_resolution.rst
    api_reference/calc_sensitivity_kernel.rst
    api_reference/combine_harmonics.rst
    api_reference/convert_harmonics.rst
    api_reference/dealiasing_global_uplift.rst
    api_reference/dealiasing_monthly_mean.rst
    api_reference/grace_mean_harmonics.rst
    api_reference/grace_spatial_error.rst
    api_reference/grace_spatial_maps.rst
    api_reference/mascon_reconstruct.rst
    api_reference/monte_carlo_degree_one.rst
    api_reference/quick_mascon_regress.rst
    api_reference/piecewise_grace_maps.rst
    api_reference/regress_grace_maps.rst
    api_reference/run_grace_date.rst
    api_reference/run_sea_level_equation.rst
    api_reference/scale_grace_maps.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Graphing

    api_reference/plot_AIS_grid_maps.rst
    api_reference/plot_AIS_grid_3maps.rst
    api_reference/plot_AIS_grid_4maps.rst
    api_reference/plot_AIS_grid_movie.rst
    api_reference/plot_AIS_GrIS_maps.rst
    api_reference/plot_AIS_regional_maps.rst
    api_reference/plot_AIS_regional_movie.rst
    api_reference/plot_global_grid_maps.rst
    api_reference/plot_global_grid_3maps.rst
    api_reference/plot_global_grid_4maps.rst
    api_reference/plot_global_grid_5maps.rst
    api_reference/plot_global_grid_9maps.rst
    api_reference/plot_global_grid_movie.rst
    api_reference/plot_GrIS_grid_maps.rst
    api_reference/plot_GrIS_grid_3maps.rst
    api_reference/plot_GrIS_grid_movie.rst
    api_reference/quick_mascon_plot.rst
