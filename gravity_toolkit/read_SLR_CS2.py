#!/usr/bin/env python
u"""
read_SLR_CS2.py
Written by Hugo Lecomte (11/2020)

Reads monthly degree 2,x spherical harmonic data files from SLR

Dataset distributed by CSR
    http://download.csr.utexas.edu/pub/slr/degree_2/
        C21_S21_RL06.txt or C22_S22_RL06.txt

REFERENCE:
    Dahle, C., Murböck, M., Flechtner, F. , Dobslaw, H., Michalak, G.,
    Neumayer, K. H., Abrykosov, O., Reinhold, A., König, R., Sulzbach, R.
    and Förste C., "The GFZ GRACE RL06 Monthly Gravity Field Time Series:
    Processing Details,and Quality Assessment", Remote Sensing, 11(18), 2116, 2019.
        https://doi.org/10.3390/rs11182116

CALLING SEQUENCE:
    SLR_2m = read_SLR_CS2(SLR_file)

INPUTS:
    SLR_file:
        CSR 2,1: C21_S21_RL06.txt
        CSR 2,2: C22_S22_RL06.txt

OUTPUTS:
    datac: SLR degree 2 order x cosine stokes coefficients (C2x)
    datas: SLR degree 2 order x sine stokes coefficients (S2x)
    errorc: SLR degree 2 order x cosine stokes coefficient error (eC2x)
    errors: SLR degree 2 order x sine stokes coefficient error (eS2x)
    month: GRACE/GRACE-FO month of measurement (Apr. 2002 = 004)
    time: date of SLR measurement

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Written 11/2020
"""
import os
import re
import numpy as np

#-- PURPOSE: read Degree 2,x data from Satellite Laser Ranging (SLR)
def read_SLR_CS2(SLR_file):
    """
    Reads CS2,x spherical harmonic coefficients from SLR measurements

    Arguments
    ---------
    SLR_file: Satellite Laser Ranging file

    Returns
    -------
    datac: SLR degree 2 order x cosine stokes coefficients (C2x)
    datas: SLR degree 2 order x sine stokes coefficients (S2x)
    errorc: SLR degree 2 order x cosine stokes coefficient error (eC2x)
    errors: SLR degree 2 order x sine stokes coefficient error (eS2x)
    month: GRACE/GRACE-FO month of measurement
    time: date of SLR measurement
    """

    #-- check that SLR file exists
    if not os.access(os.path.expanduser(SLR_file), os.F_OK):
        raise IOError('SLR file not found in file system')
    #-- output dictionary with input data
    dinput = {}

    if bool(re.search('C2\d_S2\d_RL',SLR_file)):

        #-- SLR 2x RL06 file from CSR
        #-- automatically skip the header denoted with '#'
        content = np.genfromtxt(os.path.expanduser(SLR_file))

        #-- number of months within the file
        n_mon = content.shape[0]
        date_conv = content[:,0]
        #-- remove the monthly mean of the AOD model
        C2x_input = content[:,1] - content[:,5]*10**-10
        eC2x_input = content[:,3]*10**-10
        # -- remove the monthly mean of the AOD model
        S2x_input = content[:,2] - content[:,6]*10**-10
        eS2x_input = content[:,4]*10**-10
        mon = np.zeros((n_mon),dtype=np.int)

        #-- for every line convert the date into month number:
        for t in range(content.shape[0]):
            # -- GRACE/GRACE-FO month of SLR solutions
            mon[t] = 1 + t

        #-- convert to output variables and truncate if necessary
        dinput['time'] = date_conv
        dinput['datac'] = C2x_input
        dinput['errorc'] = eC2x_input
        dinput['datas'] = S2x_input
        dinput['errors'] = eS2x_input
        dinput['month'] = mon

    else:
        raise FileNotFoundError("Invalid file given to read_SLR_2x:", SLR_file)

    #-- return the input CS2x data, year-decimal date, and GRACE/GRACE-FO month
    return dinput
