#!/usr/bin/env python
u"""
grace_find_months.py
Written by Tyler Sutterley (05/2023)

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
    Updated 05/2023: use formatting for reading from date file
        use pathlib to define and operate on paths
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
import pathlib
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
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    grace_dir = base_dir.joinpath(PROC, DREL, DSET)

    # check that GRACE/GRACE-FO date file exists
    grace_date_file = grace_dir.joinpath(f'{PROC}_{DREL}_DATES.txt')
    if not grace_date_file.exists():
        grace_date(base_dir, PROC=PROC, DREL=DREL, DSET=DSET, OUTPUT=True)

    # names and formats of GRACE/GRACE-FO date ascii file
    names = ('t','mon','styr','stday','endyr','endday','total')
    formats = ('f','i','i','i','i','i','i')
    dtype = np.dtype({'names':names, 'formats':formats})
    # read GRACE/GRACE-FO date ascii file
    # skip the header row and extract dates (decimal format) and months
    date_input = np.loadtxt(grace_date_file, skiprows=1, dtype=dtype)
    # date info dictionary
    var_info = {}
    var_info['time'] = date_input['t']
    var_info['months'] = date_input['mon']
    var_info['start'] = np.min(date_input['mon'])
    var_info['end'] = np.max(date_input['mon'])

    # array of all possible months (or in case of CNES RL01/2: 10-day sets)
    all_months = np.arange(1, var_info['end'], dtype=np.int64)
    # missing months (values in all_months but not in months)
    var_info['missing'] = sorted(set(all_months) - set(date_input['mon']))
    # If CNES RL01/2: simply convert into numpy array
    # else: remove months 1-3 and convert into numpy array
    if ((PROC == 'CNES') & (DREL in ('RL01','RL02'))):
        var_info['missing'] = np.array(var_info['missing'], dtype=np.int64)
    else:
        var_info['missing'] = np.array(var_info['missing'][3:], dtype=np.int64)

    # return the date information dictionary
    return var_info
