#!/usr/bin/env python
u"""
read_ICGEM_harmonics.py
Written by Tyler Sutterley (07/2020)

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
    Updated 12/2020: added GRAZ information extraction
    Updated 07/2020: added function docstrings
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import io
import numpy as np

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_ICGEM_harmonics(model_file, FLAG='gfc'):
    """
    Extract gravity model spherical harmonics from GFZ/GRAZ ICGEM gfc files

    Arguments
    ---------
    model_file: GFZ/GRAZ ICGEM gfc spherical harmonic data file

    Keyword arguments
    -----------------
    FLAG: string denoting data lines

    Returns
    -------
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations of type errors
    eslm: sine spherical harmonic standard deviations of type errors
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for gravity model
    radius: semi-major axis of the Earth for gravity model
    max_degree: maximum degree and order for gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model
    """
    #-- python dictionary with model input and headers
    model_input = {}
    if 'ITSG' in model_file:
        #-- parse filename
        PFX, SAT, trunc, year, month, SFX = parse_file(model_file)
        #-- convert string to integer
        year, month = int(year), int(month)

        #-- calculate mid-month date taking into account if measurements are
        #-- on different years
        if (year % 4) == 0:  #-- Leap Year
            dpy = 366.0
            dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:  #-- Standard Year
            dpy = 365.0
            dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        start_day = np.sum(dpm[:month - 1]) + 1
        end_day = np.sum(dpm[:month])

        # -- Calculation of Mid-month value
        mid_day = np.mean([start_day, end_day])
        # -- Calculating the mid-month date in decimal form
        model_input['time'] = year + mid_day / dpy
        # -- Calculating the Julian dates of the start and end date
        model_input['start'] = np.float(367.0 * year - np.floor(7.0 * year / 4.0) -
                                           np.floor(3.0 * (np.floor((year - 8.0 / 7.0) / 100.0) + 1.0) / 4.0) +
                                           np.floor(275.0 / 9.0) + start_day + 1721028.5)
        model_input['end'] = np.float(367.0 * year - np.floor(7.0 * year / 4.0) -
                                         np.floor(3.0 * (np.floor((year - 8.0 / 7.0) / 100.0) + 1.0) / 4.0) +
                                         np.floor(275.0 / 9.0) + end_day + 1721028.5)

    #-- read input data
    with open(os.path.expanduser(model_file),'r') as f:
        file_contents = f.read().splitlines()

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

#-- PURPOSE: extract parameters from filename
def parse_file(input_file):
    """
    Extract parameters from filename

    Arguments
    ---------
    input_file: GRACE/GRACE-FO Level-2 spherical harmonic data file
    """
    #-- compile numerical expression operator for parameters from files
    # -- GRAZ: Institute of Geodesy from GRAZ University of Technology
    regex_pattern = (r'(.*?)-({0})_(.*?)_(\d+)-(\d+)'
                     r'(\.gz|\.gfc|\.txt)').format(r'Grace_operational|Grace2018')
    rx = re.compile(regex_pattern, re.VERBOSE)
    #-- extract parameters from input filename
    if isinstance(input_file, io.IOBase):
        return rx.findall(input_file.filename).pop()
    else:
        return rx.findall(os.path.basename(input_file)).pop()