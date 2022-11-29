#!/usr/bin/env python
u"""
grace_find_months.py
Written by Tyler Sutterley (11/2022)

Parses date index file from grace_date program
Finds the months available for a GRACE/GRACE-FO/Swarm product
Finds the all months missing from the product

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data
    PROC: Data processing center or satellite mission
        CSR: University of Texas Center for Space Research
        GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        JPL: Jet Propulsion Laboratory
        CNES: French Centre National D'Etudes Spatiales
        GRAZ: Institute of Geodesy from GRAZ University of Technology
        COSTG: Combination Service for Time-variable Gravity Fields
        Swarm: Time-variable gravity data from Swarm satellites
    DREL: GRACE/GRACE-FO/Swarm data release

OPTIONS:
    DSET: GRACE/GRACE-FO/Swarm dataset (GSM, GAC, GAD, GAB, GAA)

OUTPUTS:
    start: First month in a GRACE/GRACE-FO dataset
    end: Last month in a GRACE/GRACE-FO dataset
    missing: missing months in a GRACE/GRACE-FO dataset
    months: all available months in a GRACE/GRACE-FO dataset
    time: center dates of all available months in a GRACE/GRACE-FO dataset

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    grace_date.py: reads GRACE index file and calculates dates for each month

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 07/2020: added function docstrings
    Updated 03/2020: check that GRACE/GRACE-FO date file exists
    Updated 10/2019: using local() function to set subdirectories
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 05/2016: using __future__ print function
    Updated 11/2015: simplified to use sets for find missing months
    Updated 09/2014: add CNES versions as RL03 is monthly data
    Updated 09/2013: missing periods for for CNES
    Written 05/2013
"""
from __future__ import print_function

import os
import numpy as np
from gravity_toolkit.grace_date import grace_date

def grace_find_months(base_dir, PROC, DREL, DSET='GSM'):
    """
    Parses date index file

    Finds the months available for a GRACE/GRACE-FO/Swarm product

    Finds the all months missing from the product

    Parameters
    ----------
    base_dir: str
        working data directory
    PROC: str
        GRACE data processing center

            - ``'CSR'``: University of Texas Center for Space Research
            - ``'GFZ'``: German Research Centre for Geosciences (GeoForschungsZentrum)
            - ``'JPL'``: Jet Propulsion Laboratory
            - ``'CNES'``: French Centre National D'Etudes Spatiales
            - ``'GRAZ'``: Institute of Geodesy from GRAZ University of Technology
            - ``'COSTG'``: Combination Service for Time-variable Gravity Fields
            - ``'Swarm'``: Time-variable gravity data from Swarm satellites
    DREL: str
        GRACE/GRACE-FO/Swarm data release

    DSET: str, default 'GSM'
        GRACE/GRACE-FO/Swarm dataset

            - ``'GAA'``: non-tidal atmospheric correction
            - ``'GAB'``: non-tidal oceanic correction
            - ``'GAC'``: combined non-tidal atmospheric and oceanic correction
            - ``'GAD'``: ocean bottom pressure product
            - ``'GSM'``: corrected monthly static gravity field product

    Returns
    -------
    start: int
        First month in a GRACE/GRACE-FO dataset
    end: int
        Last month in a GRACE/GRACE-FO dataset
    missing: list
        missing months in a GRACE/GRACE-FO dataset
    months: list
        all available months in a GRACE/GRACE-FO dataset
    time: list
        center dates of all available months in a GRACE/GRACE-FO dataset
    """

    #  Directory of exact product (using date index from GSM)
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)

    # check that GRACE/GRACE-FO date file exists
    date_file = os.path.join(grace_dir, f'{PROC}_{DREL}_DATES.txt')
    if not os.access(date_file, os.F_OK):
        grace_date(base_dir, PROC=PROC, DREL=DREL, DSET=DSET, OUTPUT=True)

    # read GRACE/GRACE-FO date ascii file from grace_date.py
    # skip the header row and extract dates (decimal format) and months
    date_input = np.loadtxt(date_file, skiprows=1)
    tdec = date_input[:,0]
    months = date_input[:,1].astype(np.int64)

    # array of all possible months (or in case of CNES RL01/2: 10-day sets)
    all_months = np.arange(1,months.max(),dtype=np.int64)
    # missing months (values in all_months but not in months)
    missing = sorted(set(all_months)-set(months))
    # If CNES RL01/2: simply convert into numpy array
    # else: remove months 1-3 and convert into numpy array
    if ((PROC == 'CNES') & (DREL in ('RL01','RL02'))):
        missing = np.array(missing,dtype=np.int64)
    else:
        missing = np.array(missing[3:],dtype=np.int64)

    return {'time':tdec, 'start':months[0], 'end':months[-1], 'months':months,
        'missing':missing}
