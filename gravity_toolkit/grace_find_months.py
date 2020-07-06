#!/usr/bin/env python
u"""
grace_find_months.py
Written by Tyler Sutterley (07/2020)

Parses date index file from grace_date.py
Finds the months available for a GRACE/GRACE-FO product
Finds the all months missing from the product

INPUTS:
    base_dir: Working data directory for GRACE/GRACE-FO data
    PROC: GRACE/GRACE-FO data processing center (CSR, CNES, JPL, GFZ)
    DREL: GRACE/GRACE-FO data release (RL04, RL05, RL06)

OPTIONS:
    DSET: GRACE dataset (GSM, GAC, GAD, GAB, GAA)

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
    Parses date index file from grace_date.py
    Finds the months available for a GRACE/GRACE-FO product
    Finds the all months missing from the product

    Arguments
    ---------
    base_dir: working data directory
    PROC: GRACE data processing center
        CSR: University of Texas Center for Space Research
        GFZ: German Research Centre for Geosciences (GeoForschungsZentrum)
        JPL: Jet Propulsion Laboratory
        CNES: French Centre National D'Etudes Spatiales
    DREL: GRACE/GRACE-FO data release

    Keyword arguments
    -----------------
    DSET: GRACE/GRACE-FO dataset
        GAA: non-tidal atmospheric correction
        GAB: non-tidal oceanic correction
        GAC: combined non-tidal atmospheric and oceanic correction
        GAD: ocean bottom pressure product
        GSM: corrected monthly static gravity field product

    Returns
    -------
    start: First month in a GRACE/GRACE-FO dataset
    end: Last month in a GRACE/GRACE-FO dataset
    missing: missing months in a GRACE/GRACE-FO dataset
    months: all available months in a GRACE/GRACE-FO dataset
    time: center dates of all available months in a GRACE/GRACE-FO dataset
    """

    #--  Directory of exact product (using date index from GSM)
    grace_dir = os.path.join(base_dir, PROC, DREL, DSET)

    #-- check that GRACE/GRACE-FO date file exists
    date_file = os.path.join(grace_dir,'{0}_{1}_DATES.txt'.format(PROC, DREL))
    if not os.access(date_file, os.F_OK):
        grace_date(base_dir,PROC=PROC,DREL=DREL,DSET=DSET,OUTPUT=True)

    #-- read GRACE/GRACE-FO date ascii file from grace_date.py
    #-- skip the header row and extract dates (decimal format) and months
    date_input = np.loadtxt(date_file, skiprows=1)
    tdec = date_input[:,0]
    months = date_input[:,1].astype(np.int)

    #-- array of all possible months (or in case of CNES RL01/2: 10-day sets)
    all_months = np.arange(1,months.max(),dtype=np.int)
    #-- missing months (values in all_months but not in months)
    missing = sorted(set(all_months)-set(months))
    #-- If CNES RL01/2: simply convert into numpy array
    #-- else: remove months 1-3 and convert into numpy array
    if ((PROC == 'CNES') & (DREL in ('RL01','RL02'))):
        missing = np.array(missing,dtype=np.int)
    else:
        missing = np.array(missing[3:],dtype=np.int)

    return {'time':tdec, 'start':months[0], 'end':months[-1], 'months':months,
        'missing':missing}
