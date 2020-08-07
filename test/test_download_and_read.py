#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
"""
import warnings
import pytest
import gravity_toolkit.utilities
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics

#-- PURPOSE: Download a GRACE file from PO.DAAC and check that read program runs
def test_download_and_read(username,webdav):
    HOST=['https://podaac-tools.jpl.nasa.gov','drive','files','allData','grace',
        'L2','CSR','RL06','GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    gravity_toolkit.utilities.from_podaac(HOST,username=username,
        password=webdav,local=HOST[-1],verbose=True)
    Ylms = read_GRACE_harmonics(HOST[-1], 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
