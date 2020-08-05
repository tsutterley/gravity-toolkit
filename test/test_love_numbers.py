#!/usr/bin/env python
u"""
test_love_numbers.py (08/2020)
"""
import os
import inspect
import warnings
import pytest
import gravity_toolkit.read_love_numbers

#-- PURPOSE: Define Load Love Numbers for lower degree harmonics
def get_love_numbers():
    """
    Gets a list of Load Love Numbers for degrees 0 to 3
    """
    hl = [-0.13273, -0.26052933922888, -0.99015777857079, -1.0499804506108]
    kl = [0.0, 0.02743231435436, -0.30252982142510, -0.19413374466744]
    ll = [0.0, 0.13026466961444, 0.023882296795977 , 0.069842389427609]
    return dict(hl=hl,kl=kl,ll=ll)

#-- PURPOSE: Check that Load Love Numbers match expected for reference frame
def test_love_numbers():
    # valid low degree Love numbers for reference frame CF
    VALID = get_love_numbers()
    # path to load Love numbers file
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    love_numbers_file = os.path.join(filepath,'..','love_numbers')
    # read load Love numbers and convert to reference frame CF
    TEST = gravity_toolkit.read_love_numbers(os.path.relpath(love_numbers_file),
        FORMAT='dict',REFERENCE='CF')
    assert all((v==t).all() for key in ['hl','kl','ll'] for v,t in zip(VALID[key],TEST[key]))
