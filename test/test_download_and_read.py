#!/usr/bin/env python
u"""
test_download_and_read.py (11/2021)
Tests the read program to verify that coefficients are being extracted
"""
import pytest
import pathlib
import posixpath
import gravity_toolkit as gravtk

# PURPOSE: Download a GRACE file from PO.DAAC Cumulus and check that read program runs
def test_podaac_cumulus_download_and_read(username,password):
    # find the path to the data files
    ids, urls, mtimes = gravtk.utilities.cmr(
        mission='grace', center='CSR', release='RL06', level='L2',
        product='GSM', start_date='2002-04-01', end_date='2002-04-30',
        provider='POCLOUD', endpoint='data')
    # attempt to download the GRACE file
    try:
        # build opener for data client access
        URS = 'urs.earthdata.nasa.gov'
        opener = gravtk.utilities.attempt_login(URS,
            username=username, password=password,
            authorization_header=False, verbose=True)
        # download and read as virtual file object
        FILE = gravtk.utilities.from_http(urls[0], context=None, verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    except EOFError as exc:
        pytest.xfail("NASA Earthdata Login Error")
    # read as virtual file object
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

# PURPOSE: Download a GRACE file from GFZ and check that read program runs
def test_gfz_http_download_and_read():
    HOST=['https://isdc-data.gfz.de','grace','Level-2','CSR','RL06',
        'GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    # attempt to download the GRACE file
    try:
        # download and read as virtual file object
        FILE = gravtk.utilities.from_http(HOST,verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    # read as virtual file object
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

# PURPOSE: Download a GRACE file from GFZ and check that read program runs
@pytest.mark.skip(reason="Deprecated GFZ FTP server")
def test_gfz_ftp_download_and_read():
    HOST=['isdcftp.gfz-potsdam.de','grace','Level-2','CSR','RL06',
        'GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    # attempt to download the GRACE file
    try:
        # download and read as virtual file object
        FILE = gravtk.utilities.from_ftp(HOST, verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    # read as virtual file object
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

# PURPOSE: Download a GRACE-FO COST-G file from the GFZ ICGEM
def test_gfz_icgem_costg_download_and_read():
    HOST=['https://icgem.gfz.de','getseries','02_COST-G_',
        'Grace-FO_RL02','GSM-2_2018152-2018181_GRFO_COSTG_BF01_0200.gfc']
    # attempt to download the GRACE file
    try:
        # download and read as virtual file object
        FILE = gravtk.utilities.from_http(HOST,verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    # read as virtual file object
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2458270.5, end=2458299.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484165368910e-03)

# PURPOSE: Download a Swarm file from ESA and check that read program runs
def test_esa_swarm_download_and_read():
    # build url for Swarm file
    HOST='https://swarm-diss.eo.esa.int'
    swarm_file='SW_OPER_EGF_SHA_2__20131201T000000_20131231T235959_0101.ZIP'
    parameters = gravtk.utilities.urlencode({'file':
        posixpath.join('swarm','Level2longterm','EGF',swarm_file)})
    remote_file = [HOST,'?do=download&{0}'.format(parameters)]    
    # attempt to download the Swarm file
    try:
        # download as local file object
        gravtk.utilities.from_http(remote_file,
            local=swarm_file, verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    # read the local file
    swarm_file = pathlib.Path(swarm_file).absolute()
    Ylms = gravtk.read_gfc_harmonics(swarm_file)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm']
    test = dict(start=2456627.5, end=2456658.499988426)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.48416530506600003e-03)
    # clean up
    swarm_file.unlink()

# PURPOSE: Download a GRACE ITSG GRAZ file and check that read program runs
def test_itsg_graz_download_and_read():
    HOST=['http://ftp.tugraz.at','pub','ITSG','GRACE',
        'ITSG-Grace_operational','monthly','monthly_n60',
        'ITSG-Grace_operational_n60_2018-06.gfc']
    # attempt to download the GRAZ file
    try:
        # download as local file object
        gravtk.utilities.from_http(HOST, local=HOST[-1], verbose=True)
    except gravtk.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    # read the local file
    itsg_file = pathlib.Path(HOST[-1]).absolute()
    Ylms = gravtk.read_gfc_harmonics(itsg_file)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm']
    test = dict(start=2458270.5, end=2458300.499988426)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.4841694727612e-03)
    # clean up
    itsg_file.unlink()
