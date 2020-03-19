#!/usr/bin/env python
u"""
grace_find_months.py
Written by Tyler Sutterley (03/2020)

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
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

PROGRAM DEPENDENCIES:
    grace_date.py: reads GRACE index file and calculates dates for each month

UPDATE HISTORY:
    Updated 03/2020 for public release
"""
from __future__ import print_function

import os
import numpy as np
from gravity_toolkit.grace_date import grace_date

def grace_find_months(base_dir, PROC, DREL, DSET='GSM'):
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
