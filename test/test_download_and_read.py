#!/usr/bin/env python
u"""
test_download_and_read.py (11/2021)
Tests the read program to verify that coefficients are being extracted
"""
import os
import pytest
import shutil
import inspect
import warnings
import posixpath
import gravity_toolkit.geocenter
import gravity_toolkit.utilities
from gravity_toolkit.read_gfc_harmonics import read_gfc_harmonics
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
def test_gfz_ftp_download_and_read():
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

#-- PURPOSE: Download a GRACE-FO COST-G file from the GFZ ICGEM
def test_gfz_icgem_costg_download_and_read():
    #-- attempt to download from ftp server
    try:
        HOST=['icgem.gfz-potsdam.de','02_COST-G','Grace-FO',
            'GSM-2_2018152-2018181_GRFO_COSTG_BF01_0100.gfc']
        FILE = gravity_toolkit.utilities.from_ftp(HOST,verbose=True)
    except:
        pass
    #-- attempt to download from http server
    try:
        HOST=['http://icgem.gfz-potsdam.de','getseries','02_COST-G',
            'Grace-FO','GSM-2_2018152-2018181_GRFO_COSTG_BF01_0100.gfc']
        FILE = gravity_toolkit.utilities.from_http(HOST,verbose=True)
    except:
        return
    #-- read as virtual file object
    Ylms = read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2458270.5, end=2458299.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484165436067e-03)

#-- PURPOSE: Download a Swarm file from ESA and check that read program runs
def test_esa_swarm_download_and_read():
    #-- build url for Swarm file
    HOST='https://swarm-diss.eo.esa.int'
    swarm_file='SW_OPER_EGF_SHA_2__20131201T000000_20131231T235959_0101.ZIP'
    parameters = gravity_toolkit.utilities.urlencode({'file':
        posixpath.join('swarm','Level2longterm','EGF',swarm_file)})
    remote_file = [HOST,'?do=download&{0}'.format(parameters)]
    #-- download and read as virtual file object
    gravity_toolkit.utilities.from_http(remote_file,
        local=swarm_file,verbose=True)
    Ylms = read_gfc_harmonics(swarm_file)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm']
    test = dict(start=2456627.5, end=2456658.499988426)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.48416530506600003e-03)
    #-- clean up
    os.remove(swarm_file)

#-- PURPOSE: Download a GRACE ITSG GRAZ file and check that read program runs
def test_itsg_graz_download_and_read():
    HOST=['http://ftp.tugraz.at','outgoing','ITSG','GRACE',
        'ITSG-Grace_operational','monthly','monthly_n60',
        'ITSG-Grace_operational_n60_2018-06.gfc']
    #-- download and read as virtual file object
    gravity_toolkit.utilities.from_http(HOST,local=HOST[-1],verbose=True)
    Ylms = read_gfc_harmonics(HOST[-1])
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm']
    test = dict(start=2458270.5, end=2458300.499988426)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.4841694727612e-03)
    #-- clean up
    os.remove(HOST[-1])

#-- PURPOSE: Download Sutterley and Velicogna (2019) geocenter files
@pytest.fixture(scope="module", autouse=True)
def download_geocenter():
    #-- download geocenter files to filepath
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    gravity_toolkit.utilities.from_figshare(filepath,verbose=True)
    #-- run tests
    yield
    #-- clean up
    shutil.rmtree(os.path.join(filepath,'geocenter'))

#-- parameterize processing center and data release
@pytest.mark.parametrize("PROC", ['CSR','GFZ','JPL'])
@pytest.mark.parametrize("DREL", ['RL06'])
#-- PURPOSE: read Sutterley and Velicogna (2019) geocenter files
def test_geocenter_read(PROC, DREL):
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    MODEL = dict(RL04='OMCT', RL05='OMCT', RL06='MPIOM')
    args = (PROC,DREL,MODEL[DREL],'SLF_iter')
    FILE = '{0}_{1}_{2}_{3}.txt'.format(*args)
    #-- assert that file exists
    geocenter_file = os.path.join(filepath,'geocenter',FILE)
    assert os.access(geocenter_file, os.F_OK)
    #-- test geocenter read program
    DEG1 = read_GRACE_geocenter(geocenter_file)
    keys = ['time', 'JD', 'month', 'C10', 'C11', 'S11','header']
    assert all((key in DEG1.keys()) for key in keys)
    #-- test geocenter class
    DATA = gravity_toolkit.geocenter().from_UCI(geocenter_file)
    for key in ['time', 'month', 'C10', 'C11', 'S11']:
        val = getattr(DATA, key)
        assert all(val == DEG1[key])
