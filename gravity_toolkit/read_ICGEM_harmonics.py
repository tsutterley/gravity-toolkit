#!/usr/bin/env python
u"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (07/2020)

Read gfc files and extract gravity model spherical harmonics from the GFZ ICGEM

GFZ International Centre for Global Earth Models (ICGEM)
    http://icgem.gfz-potsdam.de/

INPUTS:
    model_file: GFZ ICGEM gfc spherical harmonic data file

OPTIONS:
    FLAG: string denoting data lines (default gfc)

OUTPUTS:
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for the gravity model
    radius: semi-major axis of the Earth for the gravity model
    max_degree: maximum degree and order for the gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model (mean_tide, zero_tide, tide_free)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    read_ICGEM_harmonics.py Reads the coefficients for a given gravity model file
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system

UPDATE HISTORY:
    Updated 03/2021: convert to a simple wrapper function to geoid toolkit
    Updated 07/2020: added function docstrings
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import warnings

# attempt imports
try:
    import geoid_toolkit.read_ICGEM_harmonics
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("geoid_toolkit not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use geoid toolkit instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return geoid_toolkit.read_ICGEM_harmonics(*args,**kwargs)
