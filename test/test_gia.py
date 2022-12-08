#!/usr/bin/env python
u"""
test_gia.py (12/2022)
Tests the that GIA model readers are equivalent
"""
import os
import gzip
import time
import pytest
import shutil
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: Download ICE-6G GIA model
@pytest.fixture(scope="module", autouse=True)
def download_GIA_model():
    # output GIA file
    GIA_FILE = 'Stokes_trend_High_Res.txt'
    # download GIA model
    HOST = ['https://www.atmosp.physics.utoronto.ca','~peltier','datasets',
        'Ice6G_C_VM5a','ICE-6G_High_Res_Stokes_trend.txt.gz']
    fid = gravtk.utilities.from_http(HOST, verbose=True)
    # decompress GIA model from virtual BytesIO object
    with gzip.open(fid, 'rb') as f_in, open(GIA_FILE, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # run tests
    yield
    # clean up
    os.remove(GIA_FILE)

# PURPOSE: read ICE-6G GIA test outputs
def test_GIA_model_read():
    # output GIA file and type
    GIA_FILE = 'Stokes_trend_High_Res.txt'
    GIA = 'ICE6G-D'
    # read GIA model
    Ylms = gravtk.read_GIA_model(GIA_FILE, GIA=GIA)
    # assert input GIA values
    assert Ylms['clm'][2,0] == 1.43961238E-11
    assert Ylms['clm'][3,0] == 1.52009079E-12
    assert Ylms['slm'][3,1] == -8.05198489E-12
    # assert parameters
    assert Ylms['title'] == 'ICE6G-D_High_Res'

# PURPOSE: read ICE-6G GIA model and test harmonic outputs
def test_GIA_model_harmonics():
    # output GIA file and type
    GIA_FILE = 'Stokes_trend_High_Res.txt'
    GIA = 'ICE6G-D'
    # degree of truncation
    LMAX,MMAX = (60, 30)
    # read GIA model
    Ylms = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
    # assert input GIA values
    assert Ylms.clm[2,0] == 1.43961238E-11
    assert Ylms.clm[3,0] == 1.52009079E-12
    assert Ylms.slm[3,1] == -8.05198489E-12
    # assert parameters
    assert Ylms.title == 'ICE6G-D_High_Res'
    # assert truncation
    assert Ylms.lmax == LMAX
    assert Ylms.l[-1] == LMAX
    assert Ylms.mmax == MMAX
    assert Ylms.m[-1] == MMAX

# PURPOSE: read ICE-6G GIA model and compare drift estimates
def test_GIA_model_drift_estimate():
    # output GIA file and type
    GIA_FILE = 'Stokes_trend_High_Res.txt'
    GIA = 'ICE6G-D'
    # degree and order of truncation
    LMAX,MMAX = (60, 30)
    # synthetic time estimate
    now = time.gmtime()
    tdec = np.arange(2002, now.tm_year+1, 1.0/12.0)
    epoch = 2003.3
    # read GIA model
    GIA_Ylms_rate = gravtk.read_GIA_model(GIA_FILE, GIA=GIA, LMAX=LMAX, MMAX=MMAX)
    # calculate the monthly mass change from GIA
    GIA_Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    GIA_Ylms.time = np.copy(tdec)
    GIA_Ylms.month = gravtk.time.calendar_to_grace(tdec)
    GIA_Ylms.month = gravtk.time.adjust_months(GIA_Ylms.month)
    # allocate for output harmonics
    GIA_Ylms.clm = np.zeros((GIA_Ylms.lmax+1, GIA_Ylms.mmax+1, len(tdec)))
    GIA_Ylms.slm = np.zeros((GIA_Ylms.lmax+1, GIA_Ylms.mmax+1, len(tdec)))
    # assert input GIA values
    # monthly GIA calculated by gia_rate*time elapsed
    # finding change in GIA each month
    for i,t in enumerate(tdec):
        GIA_Ylms.clm[:,:,i] = GIA_Ylms_rate['clm']*(t - epoch)
        GIA_Ylms.slm[:,:,i] = GIA_Ylms_rate['slm']*(t - epoch)
    # read GIA model and calculate drift from harmonics class
    Ylms = gravtk.gia(lmax=LMAX).from_GIA(
        GIA_FILE, GIA=GIA, mmax=MMAX).drift(tdec, epoch=epoch)
    # assert that spherical harmonics are equal
    assert np.all(GIA_Ylms.clm == Ylms.clm)
    assert np.all(GIA_Ylms.slm == Ylms.slm)
    # assert time variables are equal
    assert np.all(GIA_Ylms.time == Ylms.time)
    assert np.all(GIA_Ylms.month == Ylms.month)
