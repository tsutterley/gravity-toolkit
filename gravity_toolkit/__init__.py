"""
A spherical harmonics toolkit for Python
========================================

gravity_toolkit contains Python tools for working with Level-2
spherical harmonic coefficients from the NASA/DLR GRACE and
NASA/GFZ GRACE-FO missions

The package works using scientific Python packages (numpy and scipy)
combined with data storage in ascii, netCDF4 and HDF5 and mapping with
matplotlib and cartopy

It aims to be a simple and efficient solution for using harmonic data from
the GRACE/GRACE-FO missions and to support their science applications

Documentation is available at https://gravity-toolkit.readthedocs.io
"""
import gravity_toolkit.geocenter
import gravity_toolkit.mascons
import gravity_toolkit.time
import gravity_toolkit.tools
import gravity_toolkit.utilities
import gravity_toolkit.version
from gravity_toolkit import SLR
from gravity_toolkit import time_series
from gravity_toolkit.associated_legendre import (
    associated_legendre,
    plm_colombo,
    plm_holmes,
    plm_mohlenkamp
)
from gravity_toolkit.clenshaw_summation import clenshaw_summation
from gravity_toolkit.degree_amplitude import degree_amplitude
from gravity_toolkit.destripe_harmonics import destripe_harmonics
from gravity_toolkit.fourier_legendre import (
    fourier_legendre,
    legendre_gradient
)
from gravity_toolkit.gauss_weights import gauss_weights
from gravity_toolkit.gen_averaging_kernel import gen_averaging_kernel
from gravity_toolkit.gen_disc_load import gen_disc_load
from gravity_toolkit.gen_harmonics import gen_harmonics
from gravity_toolkit.gen_point_load import gen_point_load
from gravity_toolkit.gen_spherical_cap import gen_spherical_cap
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.geocenter import geocenter
from gravity_toolkit.grace_date import grace_date
from gravity_toolkit.grace_find_months import grace_find_months
from gravity_toolkit.grace_input_months import (
    grace_input_months,
    read_ecmwf_corrections
)
from gravity_toolkit.grace_months_index import grace_months_index
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.harmonic_gradients import (
    harmonic_gradients,
    geostrophic_currents
)
from gravity_toolkit.harmonic_summation import (
    harmonic_summation,
    harmonic_transform,
    stokes_summation
)
from gravity_toolkit.legendre_polynomials import legendre_polynomials
from gravity_toolkit.legendre import legendre
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
from gravity_toolkit.read_GIA_model import (
    read_GIA_model,
    gia
)
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
from gravity_toolkit.read_love_numbers import (
    read_love_numbers,
    load_love_numbers,
    love_numbers
)
from gravity_toolkit.read_SLR_harmonics import (
    read_SLR_harmonics,
    convert_weekly
)
from gravity_toolkit.sea_level_equation import sea_level_equation
from gravity_toolkit.spatial import (
    spatial,
    scaling_factors
)
from gravity_toolkit.units import units
# get version number
__version__ = gravity_toolkit.version.version
