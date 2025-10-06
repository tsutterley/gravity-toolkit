#!/usr/bin/env python
u"""
test_download_and_read.py (11/2021)
Tests the read program to verify that coefficients are being extracted
"""
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
    # build opener for data client access
    URS = 'urs.earthdata.nasa.gov'
    opener = gravtk.utilities.attempt_login(URS,
        username=username, password=password,
        authorization_header=False, verbose=True)
    # download and read as virtual file object
    FILE = gravtk.utilities.from_http(urls[0], context=None, verbose=True)
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

# PURPOSE: Download a GRACE file from GFZ and check that read program runs
def test_gfz_ftp_download_and_read():
    HOST=['isdcftp.gfz-potsdam.de','grace','Level-2','CSR','RL06',
        'GSM-2_2002095-2002120_GRAC_UTCSR_BA01_0600.gz']
    # download and read as virtual file object
    FILE = gravtk.utilities.from_ftp(HOST,verbose=True)
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2452369.5, end=2452394.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484169355584e-03)

# PURPOSE: Download a GRACE-FO COST-G file from the GFZ ICGEM
def test_gfz_icgem_costg_download_and_read():
    # attempt to download from ftp server
    try:
        HOST=['icgem.gfz-potsdam.de','02_COST-G','Grace-FO',
            'GSM-2_2018152-2018181_GRFO_COSTG_BF01_0100.gfc']
        FILE = gravtk.utilities.from_ftp(HOST,verbose=True)
    except:
        pass
    # attempt to download from http server
    try:
        HOST=['http://icgem.gfz-potsdam.de','getseries','02_COST-G',
            'Grace-FO','GSM-2_2018152-2018181_GRFO_COSTG_BF01_0100.gfc']
        FILE = gravtk.utilities.from_http(HOST,verbose=True)
    except:
        return
    # read as virtual file object
    Ylms = gravtk.read_GRACE_harmonics(FILE, 60)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm', 'header']
    test = dict(start=2458270.5, end=2458299.5)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.484165436067e-03)

# PURPOSE: Download a Swarm file from ESA and check that read program runs
def test_esa_swarm_download_and_read():
    # build url for Swarm file
    HOST='https://swarm-diss.eo.esa.int'
    swarm_file='SW_OPER_EGF_SHA_2__20131201T000000_20131231T235959_0101.ZIP'
    parameters = gravtk.utilities.urlencode({'file':
        posixpath.join('swarm','Level2longterm','EGF',swarm_file)})
    remote_file = [HOST,'?do=download&{0}'.format(parameters)]
    # download and read as virtual file object
    gravtk.utilities.from_http(remote_file,
        local=swarm_file,verbose=True)
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
    HOST=['http://ftp.tugraz.at','outgoing','ITSG','GRACE',
        'ITSG-Grace_operational','monthly','monthly_n60',
        'ITSG-Grace_operational_n60_2018-06.gfc']
    # download and read as virtual file object
    gravtk.utilities.from_http(HOST, local=HOST[-1], verbose=True)
    itsg_file = pathlib.Path(HOST[-1]).absolute()
    Ylms = gravtk.read_gfc_harmonics(itsg_file)
    keys = ['time', 'start', 'end', 'clm', 'slm', 'eclm', 'eslm']
    test = dict(start=2458270.5, end=2458300.499988426)
    assert all((key in Ylms.keys()) for key in keys)
    assert all((Ylms[key] == val) for key,val in test.items())
    assert (Ylms['clm'][2,0] == -0.4841694727612e-03)
    # clean up
    itsg_file.unlink()
