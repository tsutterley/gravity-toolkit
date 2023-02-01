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
    hl,kl,ll = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF')
    factors = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
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

# PURPOSE: test spatial units
@pytest.mark.parametrize("LMAX", np.random.randint(60,696,size=1))
def test_spatial_units(LMAX):
    # extract arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF')
    factors = gravtk.units(lmax=LMAX).spatial(hl,kl,ll)
    # cmwe, centimeters water equivalent
    assert np.all(factors.get('cmwe') == factors.cmwe)
    # mmwe, millimeters water equivalent
    assert np.all(factors.get('mmwe') == factors.mmwe)
