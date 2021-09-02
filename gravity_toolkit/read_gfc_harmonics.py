#!/usr/bin/env python
u"""
read_gfc_harmonics.py
Written by Tyler Sutterley (09/2021)
Contributions by Hugo Lecomte

Reads gfc files and extracts spherical harmonics for SWARM and
    GRAZ GRACE/GRACE-FO data
Parses date of GRACE/GRACE-FO data from filename

GRAZ: https://www.tugraz.at/institute/ifg/downloads/gravity-field-models
SWARM: https://earth.esa.int/eogateway/missions/swarm

INPUTS:
    input_file: full path to gfc spherical harmonic data file

OPTIONS:
    TIDE: tide system of output gravity fields
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential
    FLAG: string denoting data lines

OUTPUTS:
    time: mid-month date in decimal form
    start: Julian dates of the start date
    end: Julian dates of the start date
    l: spherical harmonic degree to maximum degree of model
    m: spherical harmonic order to maximum degree of model
    clm: cosine spherical harmonics of input data
    slm: sine spherical harmonics of input data
    eclm: cosine spherical harmonic standard deviations
    eslm: sine spherical harmonic standard deviations
    modelname: name of the gravity model
    earth_gravity_constant: GM constant of the Earth for the gravity model
    radius: semi-major axis of the Earth for the gravity model
    max_degree: maximum degree and order for the gravity model
    errors: error type of the gravity model
    norm: normalization of the spherical harmonics
    tide_system: tide system of gravity model (mean_tide, zero_tide, tide_free)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system

UPDATE HISTORY:
    Updated 09/2021: forked from read_ICGEM_harmonics in geoid toolkit
        use gravity toolkit time modules and reorganize structure
        output spherical harmonic degree and order in dict
    Updated 05/2021: Add GRAZ/SWARM/COST-G ICGEM file
    Updated 03/2021: made degree of truncation LMAX a keyword argument
    Updated 07/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import io
import zipfile
import numpy as np
import gravity_toolkit.time
from geoid_toolkit.calculate_tidal_offset import calculate_tidal_offset

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_gfc_harmonics(input_file, TIDE=None, FLAG='gfc'):
    """
    Extract gravity model spherical harmonics from gfc files

    Arguments
    ---------
    input_file: full path to gfc spherical harmonic data file

    Keyword arguments
    -----------------
    TIDE: tide system of output gravity fields
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential
    FLAG: string denoting data lines

    Returns
    -------
    time: mid-month date in decimal form
    start: Julian dates of the start date
    end: Julian dates of the start date
    l: spherical harmonic degree to maximum degree of model
    m: spherical harmonic order to maximum degree of model
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
    #-- regular expression operators for ITSG data and models
    itsg_products = []
    itsg_products.append(r'atmosphere')
    itsg_products.append(r'dealiasing')
    itsg_products.append(r'Grace2014')
    itsg_products.append(r'Grace2016')
    itsg_products.append(r'Grace2018')
    itsg_products.append(r'Grace_operational')
    itsg_pattern = (r'(AOD1B_RL\d+|model|ITSG)[-_]({0})(_n\d+)?_'
        r'(\d+)-(\d+)(\.gfc)').format(r'|'.join(itsg_products))
    #-- regular expression operators for SWARM data and models
    swarm_data = r'(SW)_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)(\.gfc|.ZIP)'
    swarm_model = r'(GAA|GAB|GAC|GAD)_Swarm_(\d+)_(\d{2})_(\d{4}).(\.gfc|.ZIP)'
    #-- extract parameters for each data center and product
    if re.match(itsg_pattern, os.path.basename(input_file)):
        #-- compile numerical expression operator for parameters from files
        #-- GRAZ: Institute of Geodesy from GRAZ University of Technology
        rx = re.compile(itsg_pattern, re.VERBOSE | re.IGNORECASE)
        #-- extract parameters from input filename
        PFX,PRD,trunc,year,month,SFX = rx.findall(input_file).pop()
        #-- number of days in each month for the calendar year
        dpm = gravity_toolkit.time.calendar_days(int(year))
        #-- create start and end date lists
        start_date = [int(year),int(month),1,0,0,0]
        end_date = [int(year),int(month),dpm[int(month)-1],23,59,59]
    elif re.match(swarm_data, os.path.basename(input_file)):
        #-- compile numerical expression operator for parameters from files
        #-- SWARM: data from SWARM satellite
        rx = re.compile(swarm_data, re.VERBOSE | re.IGNORECASE)
        #-- extract parameters from input filename
        SAT,tmp,PROD,starttime,endtime,RL,SFX = rx.findall(input_file).pop()
        start_date,_ = gravity_toolkit.time.parse_date_string(starttime)
        end_date,_ = gravity_toolkit.time.parse_date_string(endtime)
        #-- number of days in each month for the calendar year
        dpm = gravity_toolkit.time.calendar_days(start_date[0])
    elif re.match(swarm_model, os.path.basename(input_file)):
        #-- compile numerical expression operator for parameters from files
        #-- SWARM: dealiasing products for SWARM data
        rx = re.compile(swarm_data, re.VERBOSE | re.IGNORECASE)
        #-- extract parameters from input filename
        PROD,trunc,month,year,SFX = rx.findall(input_file).pop()
        #-- number of days in each month for the calendar year
        dpm = gravity_toolkit.time.calendar_days(int(year))
        #-- create start and end date lists
        start_date = [int(year),int(month),1,0,0,0]
        end_date = [int(year),int(month),dpm[int(month)-1],23,59,59]

    #-- python dictionary with model input and headers
    model_input = {}

    #-- start and end day of the year
    start_day = np.sum(dpm[:start_date[1]-1]) + start_date[2] + \
        start_date[3]/24.0 + start_date[4]/1440.0 + start_date[5]/86400.0
    end_day = np.sum(dpm[:end_date[1]-1]) + end_date[2] + \
        end_date[3]/24.0 + end_date[4]/1440.0 + end_date[5]/86400.0
    #-- end date taking into account measurements taken on different years
    end_cyclic = (end_date[0]-start_date[0])*np.sum(dpm) + end_day
    #-- calculate mid-month value
    mid_day = np.mean([start_day, end_cyclic])
    #-- Calculating the mid-month date in decimal form
    model_input['time'] = start_date[0] + mid_day/np.sum(dpm)
    #-- Calculating the Julian dates of the start and end date
    model_input['start'] = 2400000.5 + \
        gravity_toolkit.time.convert_calendar_dates(start_date[0],
        start_date[1],start_date[2],hour=start_date[3],minute=start_date[4],
        second=start_date[5],epoch=(1858,11,17,0,0,0))
    model_input['end'] = 2400000.5 + \
        gravity_toolkit.time.convert_calendar_dates(end_date[0],
        end_date[1],end_date[2],hour=end_date[3],minute=end_date[4],
        second=end_date[5],epoch=(1858,11,17,0,0,0))

    #-- read data from compressed or gfc file
    if re.search('ZIP',SFX,re.IGNORECASE):
        #-- extract zip file with gfc file
        with zipfile.ZipFile(os.path.expanduser(input_file)) as zs:
            #-- find gfc file within zipfile
            gfc, = [io.BytesIO(zs.read(s)) for s in zs.namelist()
                if s.endswith('gfc')]
            #-- read input gfc data file
            file_contents = gfc.read().decode('ISO-8859-1').splitlines()
    else:
        #-- read input gfc data file
        with open(os.path.expanduser(input_file),'r') as f:
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
    #-- set degree of truncation from model
    LMAX = int(model_input['max_degree'])
    #-- output dimensions
    model_input['l'] = np.arange(LMAX+1)
    model_input['m'] = np.arange(LMAX+1)
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
        l1 = int(line_contents[1])
        m1 = int(line_contents[2])
        #-- fill gravity fields for degree and order
        model_input['clm'][l1,m1] = np.float64(line_contents[3])
        model_input['slm'][l1,m1] = np.float64(line_contents[4])
        #-- check if model contains errors
        try:
            model_input['eclm'][l1,m1] = np.float64(line_contents[5])
            model_input['eslm'][l1,m1] = np.float64(line_contents[6])
        except:
            pass
    #-- calculate the tidal offset if changing the tide system
    if TIDE in ('mean_tide','zero_tide','tide_free'):
        #-- earth parameters
        GM = np.float64(model_input['earth_gravity_constant'])
        R = np.float64(model_input['radius'])
        model_input['clm'][2,0] += calculate_tidal_offset(TIDE,GM,R,'WGS84',
            REFERENCE=model_input['tide_system'])
        #-- update attribute for tide system
        model_input['tide_system'] = TIDE
    #-- return the spherical harmonics and parameters
    return model_input
