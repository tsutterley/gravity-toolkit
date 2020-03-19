#!/usr/bin/env python
u"""
read_love_numbers.py
Written by Tyler Sutterley (03/2020)

Reads sets of load Love numbers output from PREM

INPUTS:
    love_numbers_file: Elastic load Love numbers computed using Preliminary
        Reference Earth Model (PREM) outputs as described by Han and Wahr (1995)

OUTPUTS:
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement

OPTIONS:
    HEADER: file contains header text to be skipped (default: True)
    FORMAT: format of output variables
        'dict': dictionary with variable keys as listed above
        'tuple': tuple with variable order hl,kl,ll
        'zip': aggregated variable sets

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

NOTES:
    love_numbers file must be in the base directory
    for l=1 coordinate system center is the center of mass of the system
    to change to a center of figure reference frame:
        replace kl[1] with -(hl[1]+2.0*ll[1])/3.0
        following Wahr (1998) and Trupin (1992)

UPDATE HISTORY:
    Updated 03/2020 for public release
"""
import os
import re
import numpy as np

#-- PURPOSE: read load love numbers from PREM
def read_love_numbers(love_numbers_file, HEADER=True, FORMAT='tuple'):

    #-- check that load love number data file is present in file system
    if not os.access(os.path.expanduser(love_numbers_file), os.F_OK):
        #-- raise error if love_numbers file is not found in path
        raise IOError('{0} not found'.format(love_numbers_file))

    #-- Input load love number data file and read contents
    with open(os.path.expanduser(love_numbers_file),'r') as f:
        file_contents = f.read().splitlines()

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text
    while HEADER:
        #-- file line at count
        line = file_contents[count]
        #-- find the final line within the header text
        #-- to set HEADER flag to False when found
        HEADER = not bool(re.match('\*\*\*',line))
        #-- add 1 to counter
        count += 1

    #-- compile regular expression operator to find numerical instances
    regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- maximum spherical harmonic degree in file
    #-- from the final line
    LMAX = np.int(rx.findall(file_contents[-1])[0])

    #-- vertical displacement hl
    hl = np.zeros((LMAX+1))
    #-- gravitational potential kl
    kl = np.zeros((LMAX+1))
    #-- horizontal displacement ll
    ll = np.zeros((LMAX+1))
    #-- for each line in the file (skipping the 2 header lines)
    for file_line in file_contents[count:]:
        #-- find numerical instances in line
        #-- Replacing IDL double precision exponential with
        #-- standard E exponential for kl and ll
        love_numbers = rx.findall(file_line.replace('D','E'))
        #-- spherical harmonic degree
        l = np.int(love_numbers[0])
        #-- convert love numbers to float
        hl[l] = np.float(love_numbers[1])
        kl[l] = np.float(love_numbers[2])
        ll[l] = np.float(love_numbers[3])

    #-- return love numbers in output format (default python dictionary)
    if (FORMAT == 'dict'):
        return {'kl':kl, 'hl':hl, 'll':ll}
    elif (FORMAT == 'tuple'):
        return (hl, kl, ll)
    elif (FORMAT == 'zip'):
        return zip(hl, kl, ll)
