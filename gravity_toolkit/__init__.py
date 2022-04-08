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

Documentation is available at https://read-grace-harmonics.readthedocs.io
"""
import gravity_toolkit.geocenter
import gravity_toolkit.mascons
import gravity_toolkit.time
import gravity_toolkit.utilities
import gravity_toolkit.version
from gravity_toolkit.clenshaw_summation import clenshaw_summation
from gravity_toolkit.degree_amplitude import degree_amplitude
from gravity_toolkit.destripe_harmonics import destripe_harmonics
from gravity_toolkit.fourier_legendre import fourier_legendre
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
from gravity_toolkit.grace_input_months import grace_input_months, read_ecmwf_corrections
from gravity_toolkit.grace_months_index import grace_months_index
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.legendre_polynomials import legendre_polynomials
from gravity_toolkit.legendre import legendre
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.piecewise_regress import piecewise_regress
from gravity_toolkit.plm_colombo import plm_colombo
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp
from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
from gravity_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from gravity_toolkit.read_love_numbers import read_love_numbers,load_love_numbers
from gravity_toolkit.read_SLR_C20 import read_SLR_C20
from gravity_toolkit.read_SLR_CS2 import read_SLR_CS2
from gravity_toolkit.read_SLR_C30 import read_SLR_C30
from gravity_toolkit.read_SLR_C50 import read_SLR_C50
from gravity_toolkit.read_SLR_harmonics import read_SLR_harmonics, convert_weekly
from gravity_toolkit.savitzky_golay import savitzky_golay
from gravity_toolkit.sea_level_equation import sea_level_equation
from gravity_toolkit.spatial import spatial
from gravity_toolkit.tsamplitude import tsamplitude
from gravity_toolkit.tsregress import tsregress
from gravity_toolkit.tssmooth import tssmooth
from gravity_toolkit.units import units
# get version number
__version__ = gravity_toolkit.version.version