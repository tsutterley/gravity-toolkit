#!/usr/bin/env python
u"""
test_legendre.py (11/2021)
"""
import numpy as np
import gravity_toolkit

#-- PURPOSE: test unnormalized Legendre polynomials
def test_unnormalized(l=3, x=[-1.0, -0.9, -0.8]):
    obs = gravity_toolkit.legendre(l, x)
    expected = np.array([
        [-1.00000, -0.47250, -0.08000],
        [0.00000, -1.99420, -1.98000],
        [0.00000, -2.56500, -4.32000],
        [0.00000, -1.24229, -3.24000]
    ])
    assert np.isclose(obs, expected, atol=1e-05).all()

#-- PURPOSE: test fully-normalized Legendre polynomials
def test_normalized(l=3, x=[-1.0, -0.9, -0.8]):
    obs = gravity_toolkit.legendre(l, x, NORMALIZE=True)
    expected = np.array([
        [-2.64575, -1.25012, -0.21166],
        [-0.00000, 2.15398, 2.13864],
        [0.00000, -0.87611, -1.47556],
        [-0.00000, 0.17323, 0.45180]
    ])
    assert np.isclose(obs, expected, atol=1e-05).all()

#-- PURPOSE: test fully-normalized zonal Legendre polynomials
def test_zonal(l=3, x=[-1.0, -0.9, -0.8]):
    obs,_ = gravity_toolkit.legendre_polynomials(l, x)
    expected = np.array([
        [1.00000, 1.00000, 1.00000],
        [-1.73205, -1.55885, -1.38564],
        [2.23607, 1.59879, 1.02859],
        [-2.64575, -1.25012, -0.21166],
    ])
    assert np.isclose(obs, expected, atol=1e-05).all()

#-- PURPOSE: compare fully-normalized Legendre polynomials
def test_plms(l=240, x=0.1):
    obs = gravity_toolkit.legendre(l, x, NORMALIZE=True)
    #-- calculate associated Legendre polynomials
    holmes,_ = gravity_toolkit.plm_holmes(l, x)
    colombo,_ = gravity_toolkit.plm_colombo(l, x)
    mohlenkamp = gravity_toolkit.plm_mohlenkamp(l, x)
    #-- compare Legendre polynomials
    assert np.isclose(obs, holmes[l,:]).all()
    assert np.isclose(holmes, colombo).all()
    assert np.isclose(holmes, mohlenkamp).all()
