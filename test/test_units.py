#!/usr/bin/env python
u"""
test_units.py (03/2023)
Verify spherical harmonic and spatial unit factors

UPDATE HISTORY:
    Updated 03/2023: include explicit comparison with some units
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
    cmwe = factors.rho_e*factors.rad_e*(2.0*factors.l+1.0)/(1.0+LOVE.kl)/3.0
    assert np.all(factors.get('cmwe') == factors.cmwe)
    assert np.all(cmwe == factors.cmwe)
    # mmGH, mm geoid height
    mmGH = np.ones((LMAX + 1))*(10.0*factors.rad_e)
    assert np.all(factors.get('mmGH') == factors.mmGH)
    assert np.all(mmGH == factors.mmGH)
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
    cmwe = 3.0*(1.0+LOVE.kl)/(1.0+2.0*factors.l)/(4.0*np.pi*factors.rad_e*factors.rho_e)
    assert np.all(factors.get('cmwe') == factors.cmwe)
    assert np.all(cmwe == factors.cmwe)
    # mmwe, millimeters water equivalent
    mmwe = 3.0*(1.0+LOVE.kl)/(1.0+2.0*factors.l)/(40.0*np.pi*factors.rad_e*factors.rho_e)
    assert np.all(factors.get('mmwe') == factors.mmwe)
    assert np.all(mmwe == factors.mmwe)
