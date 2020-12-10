#!/usr/bin/env python
u"""
read_love_numbers.py
Written by Tyler Sutterley (12/2020)

Reads sets of load Love numbers from PREM and applies isomorphic parameters
Can linearly extrapolate load love numbers beyond maximum degree of dataset

INPUTS:
    love_numbers_file: Elastic load Love numbers file
        computed using Preliminary Reference Earth Model (PREM) outputs

OUTPUTS:
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement

OPTIONS:
    LMAX: truncate or interpolate to maximum spherical harmonic degree
    HEADER: number of header lines to be skipped
    COLUMNS: column names of ascii file
        l: spherical harmonic degree
        hl: vertical displacement
        kl: gravitational potential
        ll: horizontal displacement
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure
        CL: Center of Surface Lateral Figure
        CH: Center of Surface Height Figure
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth (default)
    FORMAT: format of output variables
        'dict': dictionary with variable keys as listed above
        'tuple': tuple with variable order hl,kl,ll
        'zip': aggregated variable sets

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 12/2020: generalized ascii read for outputs from Gegout and Wang
        added linear interpolation of love numbers to a specified LMAX
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 03/2020: added reference frame transformations within function
    Updated 03/2020 for public release
"""
import os
import re
import numpy as np

#-- PURPOSE: read load love numbers from PREM
def read_love_numbers(love_numbers_file, LMAX=None, HEADER=2, 
    COLUMNS=['l','hl','kl','ll'], REFERENCE='CE', FORMAT='tuple'):
    """
    Reads PREM load Love numbers file and applies isomorphic parameters

    Arguments
    ---------
    love_numbers_file: Elastic load Love numbers file

    Keyword arguments
    -----------------
    LMAX: truncate or interpolate to maximum spherical harmonic degree
    HEADER: number of header lines to be skipped
    COLUMNS: column names of ascii file
        l: spherical harmonic degree
        hl: vertical displacement
        kl: gravitational potential
        ll: horizontal displacement
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure
        CL: Center of Surface Lateral Figure
        CH: Center of Surface Height Figure
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth (default)
    FORMAT: format of output variables
        'dict': dictionary with variable keys as listed above
        'tuple': tuple with variable order hl,kl,ll
        'zip': aggregated variable sets

    Returns
    -------
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement
    """

    #-- check that load love number data file is present in file system
    if not os.access(os.path.expanduser(love_numbers_file), os.F_OK):
        #-- raise error if love_numbers file is not found in path
        raise IOError('{0} not found'.format(love_numbers_file))

    #-- Input load love number data file and read contents
    with open(os.path.expanduser(love_numbers_file),'r') as f:
        file_contents = f.read().splitlines()

    #-- compile regular expression operator to find numerical instances
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- extract maximum spherical harmonic degree from final line in file
    if LMAX is None:
        LMAX = np.int(rx.findall(file_contents[-1])[COLUMNS.index('l')])

    #-- output love numbers
    love = {}
    #-- spherical harmonic degree
    love['l'] = np.arange(LMAX+1)
    #-- vertical displacement hl
    love['hl'] = np.zeros((LMAX+1))
    #-- gravitational potential kl
    love['kl'] = np.zeros((LMAX+1))
    #-- horizontal displacement ll
    love['ll'] = np.zeros((LMAX+1))
    #-- for each line in the file (skipping the 2 header lines)
    for file_line in file_contents[HEADER:]:
        #-- find numerical instances in line
        #-- replacing fortran double precision exponential
        love_numbers = rx.findall(file_line.replace('D','E'))            
        #-- spherical harmonic degree
        l = np.int(love_numbers[COLUMNS.index('l')])
        #-- truncate to spherical harmonic degree LMAX
        if (l <= LMAX):
            #-- convert love numbers to float
            #-- vertical displacement hl
            love['hl'][l] = np.float(love_numbers[COLUMNS.index('hl')])
            #-- gravitational potential kl
            love['kl'][l] = np.float(love_numbers[COLUMNS.index('kl')])
            #-- horizontal displacement ll
            love['ll'][l] = np.float(love_numbers[COLUMNS.index('ll')])

    #-- LMAX of load love numbers from Han and Wahr (1995) is 696
    #-- From Wahr (2007), can linearly extrapolate the load numbers
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696 (or LMAX of dataset)
    for lint in range(l,LMAX+1):
        #-- linearly extrapolating load love numbers
        love['hl'][lint] = 2.0*love['hl'][lint-1] - love['hl'][lint-2]
        love['kl'][lint] = 2.0*love['kl'][lint-1] - love['kl'][lint-2]
        love['ll'][lint] = 2.0*love['ll'][lint-1] - love['ll'][lint-2]

    #-- calculate isomorphic parameters for different reference frames
    #-- From Blewitt (2003), Wahr (1998), Trupin (1992) and Farrell (1972)
    if (REFERENCE.upper() == 'CF'):
        #-- Center of Surface Figure
        alpha = (love['hl'][1] + 2.0*love['ll'][1])/3.0
    elif (REFERENCE.upper() == 'CL'):
        #-- Center of Surface Lateral Figure
        alpha = love['ll'][1].copy()
    elif (REFERENCE.upper() == 'CH'):
        #-- Center of Surface Height Figure
        alpha = love['hl'][1].copy()
    elif (REFERENCE.upper() == 'CM'):
        #-- Center of Mass of Earth System
        alpha = 1.0
    elif (REFERENCE.upper() == 'CE'):
        #-- Center of Mass of Solid Earth
        alpha = 0.0
    else:
        raise Exception('Invalid Reference Frame {0}'.format(REFERENCE))
    #-- apply isomorphic parameters
    love['hl'][1] -= alpha
    love['kl'][1] -= alpha
    love['ll'][1] -= alpha

    #-- return love numbers in output format
    if (FORMAT == 'dict'):
        return love
    elif (FORMAT == 'tuple'):
        return (love['hl'], love['kl'], love['ll'])
    elif (FORMAT == 'zip'):
        return zip(love['hl'], love['kl'], love['ll'])
