#!/usr/bin/env python
u"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (07/2017)

Read gfc files and extract gravity model spherical harmonics from the GFZ ICGEM

GFZ International Centre for Global Earth Models (ICGEM)
    http://icgem.gfz-potsdam.de/

INPUTS:
    model_file: GFZ ICGEM gfc spherical harmonic data file

OPTIONS:
    FLAG: string denoting data lines (default gfc)

OUTPUTS:
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for the gravity model
    radius: semi-major axis of the Earth for the gravity model
    max_degree: maximum degree and order for the gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model (mean_tide, zero_tide, tide_free)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import numpy as np

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(model_file, FLAG='gfc'):
    #-- read input data
    with open(os.path.expanduser(model_file),'r') as f:
        file_contents = f.read().splitlines()
    #-- python dictionary with model input and headers
    model_input = {}
    #-- extract parameters from header
    header_parameters = ['modelname','earth_gravity_constant','radius',
        'max_degree','errors','norm','tide_system']
    parameters_regex = '(' + '|'.join(header_parameters) + ')'
    header = [l for l in file_contents if re.match(parameters_regex,l)]
    for line in header:
        #-- split the line into individual components
        line_contents = line.split()
        model_input[line_contents[0]] = line_contents[1]
    #-- set maximum spherical harmonic order
    LMAX = np.int(model_input['max_degree'])
    #-- allocate for each Coefficient
    model_input['clm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['slm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['eclm'] = np.zeros((LMAX+1,LMAX+1))
    model_input['eslm'] = np.zeros((LMAX+1,LMAX+1))
    #-- reduce file_contents to input data using data marker flag
    input_data = [l for l in file_contents if re.match(FLAG,l)]
    #-- for each line of data in the gravity file
    for line in input_data:
        #-- split the line into individual components replacing fortran d
        line_contents = re.sub('d','e',line,flags=re.IGNORECASE).split()
        #-- degree and order for the line
        l1 = np.int(line_contents[1])
        m1 = np.int(line_contents[2])
        #-- read spherical harmonic coefficients
        model_input['clm'][l1,m1] = np.float(line_contents[3])
        model_input['slm'][l1,m1] = np.float(line_contents[4])
        model_input['eclm'][l1,m1] = np.float(line_contents[5])
        model_input['eslm'][l1,m1] = np.float(line_contents[6])
    #-- return the spherical harmonics and parameters
    return model_input
