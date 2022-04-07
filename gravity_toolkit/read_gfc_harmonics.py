#!/usr/bin/env python
u"""
read_gfc_harmonics.py
Written by Tyler Sutterley (04/2022)
Contributions by Hugo Lecomte

Reads gfc files and extracts spherical harmonics for Swarm and
    GRAZ GRACE/GRACE-FO data
Parses date of GRACE/GRACE-FO/Swarm data from filename

GRAZ: https://www.tugraz.at/institute/ifg/downloads/gravity-field-models
Swarm: https://earth.esa.int/eogateway/missions/swarm

INPUTS:
    input_file: full path to gfc spherical harmonic data file

OPTIONS:
    TIDE: tide system of output gravity fields
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed
    FLAG: string denoting data lines

OUTPUTS:
    time: mid-month date in decimal form
    start: Julian dates of the start date
    end: Julian dates of the start date
    l: spherical harmonic degree to maximum degree of data
    m: spherical harmonic order to maximum degree of data
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
    tide_system: tide system of gravity model
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
        zero_tide: permanent direct tidal potential removed

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    read_ICGEM_harmonics.py: reads gravity model coefficients from GFZ ICGEM
    calculate_tidal_offset.py: calculates the C20 offset for a tidal system

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: forked from read_ICGEM_harmonics in geoid toolkit
        use gravity toolkit time modules and reorganize structure
    Updated 05/2021: Add GRAZ/Swarm/COST-G ICGEM file
    Updated 03/2021: made degree of truncation LMAX a keyword argument
    Updated 07/2020: added function docstrings
    Updated 07/2019: split read and wrapper funciton into separate files
    Updated 07/2017: include parameters to change the tide system
    Written 12/2015
"""
import os
import re
import numpy as np
import gravity_toolkit.time
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics

#-- PURPOSE: read spherical harmonic coefficients of a gravity model
def read_gfc_harmonics(input_file, TIDE=None, FLAG='gfc'):
    """
    Extract gravity model spherical harmonics from Gravity
    Field Coefficient (gfc) files

    Parameters
    ----------
    input_file: str
        full path to gfc spherical harmonic data file
    TIDE: string
        Permanent tide system of output gravity fields [Losch2003]_

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
            - ``'zero_tide'``: permanent direct tidal potential removed
    FLAG: str, default 'gfc'
        Flag denoting data lines

    Returns
    -------
    time: float
        mid-month date in decimal form
    start: float
        Julian dates of the start date
    end: float
        Julian dates of the start date
    l: int
        spherical harmonic degree to maximum degree of data
    m: int
        spherical harmonic order to maximum degree of data
    clm: float
        cosine spherical harmonics of input data
    slm: float
        sine spherical harmonics of input data
    eclm: float
        cosine spherical harmonic standard deviations of type errors
    eslm: float
        sine spherical harmonic standard deviations of type errors
    modelname: str
        name of the gravity model
    earth_gravity_constant: str
        GM constant of the Earth for gravity model
    radius: str
        semi-major axis of the Earth for gravity model
    max_degree: str
        maximum degree and order for gravity model
    errors: str
        error type of the gravity model
    norm: str
        normalization of the spherical harmonics
    tide_system: str
        Permanent tide system of gravity model (``'mean_tide'``, ``'zero_tide'``, ``'tide_free'``)

    Reference
    ---------
    .. [Losch2003] M. Losch and V. Seufer,
        "How to Compute Geoid Undulations (Geoid Height Relative
        to a Given Reference Ellipsoid) from Spherical Harmonic
        Coefficients for Satellite Altimetry Applications", (2003).
        `eprint ID: 11802 <http://mitgcm.org/~mlosch/geoidcookbook.pdf>`_
    """
    #-- regular expression operators for ITSG data and models
    itsg_products = []
    itsg_products.append(r'atmosphere')
    itsg_products.append(r'dealiasing')
    itsg_products.append(r'oceanBottomPressure')
    itsg_products.append(r'ocean')
    itsg_products.append(r'Grace2014')
    itsg_products.append(r'Grace2016')
    itsg_products.append(r'Grace2018')
    itsg_products.append(r'Grace_operational')
    itsg_pattern = (r'(AOD1B_RL\d+|model|ITSG)[-_]({0})(_n\d+)?_'
        r'(\d+)-(\d+)(\.gfc)').format(r'|'.join(itsg_products))
    #-- regular expression operators for Swarm data and models
    swarm_data = r'(SW)_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)(\.gfc|\.ZIP)'
    swarm_model = r'(GAA|GAB|GAC|GAD)_Swarm_(\d+)_(\d{2})_(\d{4})(\.gfc|\.ZIP)'
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
        #-- Swarm: data from Swarm satellite
        rx = re.compile(swarm_data, re.VERBOSE | re.IGNORECASE)
        #-- extract parameters from input filename
        SAT,tmp,PROD,starttime,endtime,RL,SFX = rx.findall(input_file).pop()
        start_date,_ = gravity_toolkit.time.parse_date_string(starttime)
        end_date,_ = gravity_toolkit.time.parse_date_string(endtime)
        #-- number of days in each month for the calendar year
        dpm = gravity_toolkit.time.calendar_days(start_date[0])
    elif re.match(swarm_model, os.path.basename(input_file)):
        #-- compile numerical expression operator for parameters from files
        #-- Swarm: dealiasing products for Swarm data
        rx = re.compile(swarm_data, re.VERBOSE | re.IGNORECASE)
        #-- extract parameters from input filename
        PROD,trunc,month,year,SFX = rx.findall(input_file).pop()
        #-- number of days in each month for the calendar year
        dpm = gravity_toolkit.time.calendar_days(int(year))
        #-- create start and end date lists
        start_date = [int(year),int(month),1,0,0,0]
        end_date = [int(year),int(month),dpm[int(month)-1],23,59,59]

    #-- python dictionary with model input and headers
    ZIP = bool(re.search('ZIP',SFX,re.IGNORECASE))
    model_input = read_ICGEM_harmonics(input_file, TIDE=TIDE,
        FLAG=FLAG, ZIP=ZIP)

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

    #-- return the spherical harmonics and parameters
    return model_input
