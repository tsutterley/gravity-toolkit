#!/usr/bin/env python
u"""
test_units.py (01/2023)
Verify spherical harmonic and spatial unit factors

UPDATE HISTORY:
    Written 01/2023
"""
import pytest
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: test spherical harmonic units
@pytest.mark.parametrize("LMAX", np.random.randint(60,696,size=1))
def test_harmonic_units(LMAX):
    # extract arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF',
        FORMAT='class')
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    # cmwe, centimeters water equivalent
    assert np.all(factors.get('cmwe') == factors.cmwe)
    # mmGH, mm geoid height
    assert np.all(factors.get('mmGH') == factors.mmGH)
    # mmCU, mm elastic crustal deformation
    assert np.all(factors.get('mmCU') == factors.mmCU)
    # microGal, microGal gravity perturbations
    assert np.all(factors.get('microGal') == factors.microGal)
    # mbar, equivalent surface pressure
    assert np.all(factors.get('mbar') == factors.mbar)
    # cmVCU, cm viscoelastic  crustal uplift (GIA ONLY)
    assert np.all(factors.get('cmVCU') == factors.cmVCU)

# PURPOSE: test harmonic units with different love number formats
@pytest.mark.parametrize("LMAX", np.random.randint(60,696,size=1))
def test_harmonic_love_numbers(LMAX):
    # extract arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF',
        FORMAT='tuple')
    LOVE = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF',
        FORMAT='class')
    factors_tuple = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
    factors_class = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    # cmwe, centimeters water equivalent
    assert np.all(factors_tuple.get('cmwe') == factors_class.get('cmwe'))
    # mmCU, mm elastic crustal deformation
    assert np.all(factors_tuple.get('mmCU') == factors_class.get('mmCU'))
    # mbar, equivalent surface pressure
    assert np.all(factors_tuple.get('mbar') == factors_class.get('mbar'))

# PURPOSE: test spatial units
@pytest.mark.parametrize("LMAX", np.random.randint(60,696,size=1))
def test_spatial_units(LMAX):
    # extract arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF',
        FORMAT='class')
    factors = gravtk.units(lmax=LMAX).spatial(*LOVE)
    # cmwe, centimeters water equivalent
    assert np.all(factors.get('cmwe') == factors.cmwe)
    # mmwe, millimeters water equivalent
    assert np.all(factors.get('mmwe') == factors.mmwe)
