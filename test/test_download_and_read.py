#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
Tests the read program to verify that coefficients are being extracted
"""
import os
import pytest
import inspect
import warnings
import gravity_toolkit.utilities
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics
from read_GRACE_geocenter.read_GRACE_geocenter import read_GRACE_geocenter

#-- PURPOSE: Download a GRACE file from PO.DAAC and check that read program runs
def test_podaac_download_and_read(username,webdav):
    HOST=['https://podaac-tools.jpl.nasa.gov','drive','files','allData','grace',
        'L2','CSR','RL06','GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    #-- download and read as virtual file object
    FILE = gravity_toolkit.utilities.from_drive(HOST,username=username,
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

#-- parameterize processing center and data release
@pytest.mark.parametrize("PROC", ['CSR','GFZ','JPL'])
@pytest.mark.parametrize("DREL", ['RL06'])
def test_geocenter_download_and_read(PROC, DREL):
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    gravity_toolkit.utilities.from_figshare(filepath,verbose=True)
    MODEL = dict(RL04='OMCT', RL05='OMCT', RL06='MPIOM')
    args = (PROC,DREL,MODEL[DREL],'SLF_iter')
    FILE = '{0}_{1}_{2}_{3}.txt'.format(*args)
    #-- assert that file exists
    assert os.access(os.path.join(filepath,'geocenter',FILE),os.F_OK)
    #-- test geocenter read program
    DEG1 = read_GRACE_geocenter(os.path.join(filepath,'geocenter',FILE))
    keys = ['time', 'JD', 'month', 'C10', 'C11', 'S11','header']
    assert all((key in DEG1.keys()) for key in keys)
