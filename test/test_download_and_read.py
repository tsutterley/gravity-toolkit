#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
Tests the read program to verify that coefficients are being extracted
"""
import warnings
import pytest
import gravity_toolkit.utilities
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics

#-- PURPOSE: Download a GRACE file from PO.DAAC and check that read program runs
def test_podaac_download_and_read(username,webdav):
    HOST=['https://podaac-tools.jpl.nasa.gov','drive','files','allData','grace',
        'L2','CSR','RL06','GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    #-- download and read as virtual file object
    FILE = gravity_toolkit.utilities.from_podaac(HOST,username=username,
        password=webdav,verbose=True)
    Ylms = read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

#-- PURPOSE: Download a GRACE file from GFZ and check that read program runs
def test_gfz_ftp_download_and_read(username,webdav):
    HOST=['isdcftp.gfz-potsdam.de','grace','Level-2','CSR','RL06',
        'GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    #-- download and read as virtual file object
    FILE = gravity_toolkit.utilities.from_ftp(HOST,verbose=True)
    Ylms = read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)
