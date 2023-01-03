#!/usr/bin/env python
u"""
test_love_numbers.py (01/2023)
UPDATE HISTORY:
    Updated 01/2023: single implicit import of gravity toolkit
    Updated 12/2020: add linear interpolation to degree 1000
        add tests for Gegout and Wang load Love number sets
    Written 08/2020
"""
import os
import warnings
import pytest
import gravity_toolkit as gravtk

# PURPOSE: Define Load Love Numbers for lower degree harmonics
def get_love_numbers():
    """
    Gets a list of Load Love Numbers for degrees 0 to 3
    """
    hl = [-0.13273, -0.26052933922888, -0.99015777857079, -1.0499804506108]
    kl = [0.0, 0.02743231435436, -0.30252982142510, -0.19413374466744]
    ll = [0.0, 0.13026466961444, 0.023882296795977, 0.069842389427609]
    return dict(hl=hl,kl=kl,ll=ll)

# PURPOSE: Check that Load Love Numbers match expected for reference frame
def test_love_numbers():
    # valid low degree Love numbers for reference frame CF
    VALID = get_love_numbers()
    # path to load Love numbers file
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','love_numbers'])
    # read load Love numbers and convert to reference frame CF
    TEST = gravtk.read_love_numbers(love_numbers_file,
        LMAX=1000, HEADER=2, FORMAT='dict', REFERENCE='CF')
    assert all((v==t).all() for key in ['hl','kl','ll']
        for v,t in zip(VALID[key],TEST[key]))
    assert (TEST['l'].max() == 1000)

# PURPOSE: Check that Gegout (2005) Load Love Numbers can be read
def test_Gegout_love_numbers():
    # path to load Love numbers file
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','Load_Love2_CE.dat'])
    COLUMNS = ['l','hl','ll','kl']
    # read load Love numbers and convert to reference frame CM
    TEST = gravtk.read_love_numbers(love_numbers_file,
        HEADER=3, COLUMNS=COLUMNS, FORMAT='dict', REFERENCE='CM')
    assert (TEST['l'].max() == 1024)

# PURPOSE: Check that Wang et al. (2012) Load Love Numbers can be read
def test_Wang_love_numbers():
    # path to load Love numbers file (truncated from degree 46341)
    love_numbers_file = gravtk.utilities.get_data_path(
        ['data','PREM-LLNs-truncated.dat'])
    COLUMNS = ['l','hl','ll','kl','nl','nk']
    # read load Love numbers and convert to reference frame CE
    TEST = gravtk.read_love_numbers(love_numbers_file,
        HEADER=1, COLUMNS=COLUMNS, FORMAT='dict', REFERENCE='CE')
    assert (TEST['l'].max() == 5000)
