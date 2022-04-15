#!/usr/bin/env python
u"""
read_love_numbers.py
Written by Tyler Sutterley (04/2022)

Reads sets of load Love numbers from PREM and applies isomorphic parameters
Linearly interpolates load love numbers for missing degrees
Linearly extrapolates load love numbers beyond maximum degree of dataset

INPUTS:
    love_numbers_file: Elastic load Love numbers file
        computed using Preliminary Reference Earth Model (PREM) outputs

OUTPUTS:
    kl: Love number of Gravitational Potential
    hl: Love number of Vertical Displacement
    ll: Love number of Horizontal Displacement

OPTIONS:
    LMAX: truncate or interpolate to maximum spherical harmonic degree
    HEADER: number of header lines to be skipped
    COLUMNS: column names of ascii file
        l: spherical harmonic degree
        hl: vertical displacement
        kl: gravitational potential
        ll: horizontal displacement
    REFERENCE: Reference frame for calculating degree 1 love numbers
        CF: Center of Surface Figure
        CL: Center of Surface Lateral Figure
        CH: Center of Surface Height Figure
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth (default)
    FORMAT: format of output variables
        'dict': dictionary with variable keys as listed above
        'tuple': tuple with variable order hl,kl,ll
        'zip': aggregated variable sets

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files

REFERENCES:
    G. Blewitt, "Self‐consistency in reference frames, geocenter definition,
        and surface loading of the solid Earth",
        Journal of Geophysical Research: Solid Earth, 108(B2), 2103, (2003)
    W. E. Farrell, "Deformation of the Earth by surface loads",
        Reviews of Geophysics, 10(3), 761--797, (1972)
    A. S. Trupin, M. F. Meier, and J. Wahr, "Effect of melting glaciers
        on the Earth's rotation and gravitational field: 1965-1984"
        Geophysical Journal International, 108(1), (1992)
    J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's
        gravity field: Hydrological and oceanic effects and their possible
        detection using GRACE", Journal of Geophysical Research,
        103(B12), 30205-30229, (1998)

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        added wrapper function for reading load Love numbers from file
        added function for checking if BytesIO object and extracting contents
        include utf-8 encoding in reads to be windows compliant
    Updated 07/2021: added check if needing to interpolate love numbers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use file not found exceptions
    Updated 12/2020: generalized ascii read for outputs from Gegout and Wang
        added linear interpolation of love numbers to a specified LMAX
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 07/2020: added function docstrings
    Updated 03/2020: added reference frame transformations within function
    Updated 03/2020 for public release
    Updated 07/2017: added FORMAT option to change output format
    Updated 06/2016: added check for love_numbers file within file system
        added while loop for skipping header text and option HEADER
    Updated 03/2015: using regular expressions and generic read
        Updated comments
    Updated 10/2013: minor changes to use the numpy genfromtxt function
    Updated 05/2013: python updates and comment updates
    Written 01/2012
"""
import os
import io
import re
import numpy as np
from gravity_toolkit.utilities import get_data_path

#-- PURPOSE: read load love numbers from PREM
def read_love_numbers(love_numbers_file, LMAX=None, HEADER=2,
    COLUMNS=['l','hl','kl','ll'], REFERENCE='CE', FORMAT='tuple'):
    """
    Reads PREM load Love numbers file and applies isomorphic
    parameters [Dziewonski1981]_ [Blewett2003]_

    Parameters
    ----------
    love_numbers_file: str
        Elastic load Love numbers file
    LMAX: int or NoneType, default None
        Truncate or interpolate to maximum spherical harmonic degree
    HEADER: int, default 2
        Number of header lines to be skipped
    COLUMNS: list
        Column names of ascii file

            - ``'l'``: spherical harmonic degree
            - ``'hl'``: vertical displacement
            - ``'kl'``: gravitational potential
            - ``'ll'``: horizontal displacement
    REFERENCE: str, default 'CE'
        Reference frame of degree 1 love numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    FORMAT: str, default 'tuple'
        Format of output variables

            - ``'dict'``: dictionary with variable keys as listed above
            - ``'tuple'``: tuple with variable order (``hl``, ``kl``, ``ll``)
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    kl: float
        Love number of Gravitational Potential
    hl: float
        Love number of Vertical Displacement
    ll: float
        Love number of Horizontal Displacement

    References
    ----------
    .. [Blewett2003] G. Blewitt, "Self‐consistency in reference frames, geocenter
        definition, and surface loading of the solid Earth",
        *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003).
        `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

    .. [Dziewonski1981] A. M. Dziewonski and D. L. Anderson,
        "Preliminary reference Earth model",
        *Physics of the Earth and Planetary Interiors*, 25(4), 297--356, (1981).
        `doi: 10.1016/0031-9201(81)90046-7 <https://doi.org/10.1016/0031-9201(81)90046-7>`_

    .. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya,
        "Practical numerical computation of love numbers and applications",
        Workshop of the COST Action ES0701, (2010).
        `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

    .. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a
        realistically stratified earth, and a further analysis of postglacial
        rebound", *Geophysical Journal International*, 120(2), 287--311, (1995).
        `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

    .. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan,
        "Time variability of the Earth's gravity field: Hydrological and
        oceanic effects and their possible detection using GRACE",
        *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998).
        `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

    .. [Wang2012] H. Wang et al., "Load Love numbers and Green's
        functions for elastic Earth models PREM, iasp91, ak135, and
        modified models with refined crustal structure from Crust 2.0",
        *Computers & Geosciences*, 49, 190--199, (2012).
        `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
    """
    #-- Input load love number data file and read contents
    file_contents = extract_love_numbers(love_numbers_file)

    #-- compile regular expression operator to find numerical instances
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- extract maximum spherical harmonic degree from final line in file
    if LMAX is None:
        LMAX = np.int64(rx.findall(file_contents[-1])[COLUMNS.index('l')])

    #-- dictionary of output love numbers
    love = {}
    #-- spherical harmonic degree
    love['l'] = np.arange(LMAX+1)
    #-- vertical displacement hl
    #-- gravitational potential kl
    #-- horizontal displacement ll
    for n in ('hl','kl','ll'):
        love[n] = np.zeros((LMAX+1))
    #-- check if needing to interpolate between degrees
    flag = np.ones((LMAX+1),dtype=bool)
    #-- for each line in the file (skipping header lines)
    for file_line in file_contents[HEADER:]:
        #-- find numerical instances in line
        #-- replacing fortran double precision exponential
        love_numbers = rx.findall(file_line.replace('D','E'))
        #-- spherical harmonic degree
        l = np.int64(love_numbers[COLUMNS.index('l')])
        #-- truncate to spherical harmonic degree LMAX
        if (l <= LMAX):
            #-- convert love numbers to float
            #-- vertical displacement hl
            #-- gravitational potential kl
            #-- horizontal displacement ll
            for n in ('hl','kl','ll'):
                love[n][l] = np.float64(love_numbers[COLUMNS.index(n)])
            #-- set interpolation flag for degree
            flag[l] = False

    #-- if needing to linearly interpolate love numbers
    if np.any(flag):
        #-- linearly interpolate each load love number following Wahr (1998)
        for n in ('hl','kl','ll'):
            love[n][flag] = np.interp(love['l'][flag],
                love['l'][~flag], love[n][~flag])

    #-- if needing to linearly extrapolate love numbers
    #-- NOTE: use caution if extrapolating far beyond the
    #-- maximum degree of the love numbers dataset
    for lint in range(l,LMAX+1):
        #-- linearly extrapolate each load love number
        for n in ('hl','kl','ll'):
            love[n][lint] = 2.0*love[n][lint-1] - love[n][lint-2]

    #-- calculate isomorphic parameters for different reference frames
    #-- From Blewitt (2003), Wahr (1998), Trupin (1992) and Farrell (1972)
    if (REFERENCE.upper() == 'CF'):
        #-- Center of Surface Figure
        alpha = (love['hl'][1] + 2.0*love['ll'][1])/3.0
    elif (REFERENCE.upper() == 'CL'):
        #-- Center of Surface Lateral Figure
        alpha = love['ll'][1].copy()
    elif (REFERENCE.upper() == 'CH'):
        #-- Center of Surface Height Figure
        alpha = love['hl'][1].copy()
    elif (REFERENCE.upper() == 'CM'):
        #-- Center of Mass of Earth System
        alpha = 1.0
    elif (REFERENCE.upper() == 'CE'):
        #-- Center of Mass of Solid Earth
        alpha = 0.0
    else:
        raise Exception('Invalid Reference Frame {0}'.format(REFERENCE))
    #-- apply isomorphic parameters
    for n in ('hl','kl','ll'):
        love[n][1] -= alpha

    #-- return love numbers in output format
    if (FORMAT == 'dict'):
        return love
    elif (FORMAT == 'tuple'):
        return (love['hl'], love['kl'], love['ll'])
    elif (FORMAT == 'zip'):
        return zip(love['hl'], love['kl'], love['ll'])

#-- PURPOSE: read input file and extract contents
def extract_love_numbers(love_numbers_file):
    """
    Read load love number file and extract contents

    Parameters
    ----------
    love_numbers_file: str
        Elastic load Love numbers file
    """
    #-- check if input love numbers are a string or bytesIO object
    if isinstance(love_numbers_file, str):
        #-- tilde expansion of load love number data file
        love_numbers_file = os.path.expanduser(love_numbers_file)
        #-- check that load love number data file is present in file system
        if not os.access(love_numbers_file, os.F_OK):
            raise FileNotFoundError('{0} not found'.format(love_numbers_file))
        #-- Input load love number data file and read contents
        with open(love_numbers_file, mode='r', encoding='utf8') as f:
            return f.read().splitlines()
    elif isinstance(love_numbers_file, io.IOBase):
        #-- read contents from load love number data
        return love_numbers_file.read().decode('utf8').splitlines()

#-- PURPOSE: read load love numbers for a range of spherical harmonic degrees
def load_love_numbers(LMAX, LOVE_NUMBERS=0, REFERENCE='CF', FORMAT='tuple'):
    """
    Wrapper function for reading PREM load Love numbers for a
    range of spherical harmonic degrees and applying
    isomorphic parameters [Blewett2003]_

    Parameters
    ----------
    LMAX: int
        maximum spherical harmonic degree
    LOVE_NUMBERS: int, default 0
        Treatment of the Load Love numbers

            - ``0``: [Han1995]_ values from PREM
            - ``1``: [Gegout2010]_ values from PREM
            - ``2``: [Wang2012]_ values from PREM
    REFERENCE: str
        Reference frame for calculating degree 1 love numbers [Blewett2003]_

            - ``'CF'``: Center of Surface Figure (default)
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    FORMAT: str, default 'tuple'
        Format of output variables

            - ``'dict'``: dictionary with variable keys as listed above
            - ``'tuple'``: tuple with variable order (``hl``, ``kl``, ``ll``)
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    kl: float
        Love number of Gravitational Potential
    hl: float
        Love number of Vertical Displacement
    ll: float
        Love number of Horizontal Displacement

    References
    ----------
    .. [Blewett2003] G. Blewitt, "Self‐consistency in reference frames, geocenter
        definition, and surface loading of the solid Earth",
        *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003).
        `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

    .. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya,
        "Practical numerical computation of love numbers and applications",
        Workshop of the COST Action ES0701, (2010).
        `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

    .. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a
        realistically stratified earth, and a further analysis of postglacial
        rebound", *Geophysical Journal International*, 120(2), 287--311, (1995).
        `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

    .. [Wang2012] H. Wang et al., "Load Love numbers and Green's
        functions for elastic Earth models PREM, iasp91, ak135, and
        modified models with refined crustal structure from Crust 2.0",
        *Computers & Geosciences*, 49, 190--199, (2012).
        `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
    """
    #-- load love numbers file
    if (LOVE_NUMBERS == 0):
        #-- PREM outputs from Han and Wahr (1995)
        #-- https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
        love_numbers_file = get_data_path(['data','love_numbers'])
        header = 2
        columns = ['l','hl','kl','ll']
    elif (LOVE_NUMBERS == 1):
        #-- PREM outputs from Gegout (2005)
        #-- http://gemini.gsfc.nasa.gov/aplo/
        love_numbers_file = get_data_path(['data','Load_Love2_CE.dat'])
        header = 3
        columns = ['l','hl','ll','kl']
    elif (LOVE_NUMBERS == 2):
        #-- PREM outputs from Wang et al. (2012)
        #-- https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREM-LLNs-truncated.dat'])
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    #-- LMAX of load love numbers from Han and Wahr (1995) is 696.
    #-- from Wahr (2007) linearly interpolating kl works
    #-- however, as we are linearly extrapolating out, do not make
    #-- LMAX too much larger than 696
    #-- read arrays of kl, hl, and ll Love Numbers
    return read_love_numbers(love_numbers_file, LMAX=LMAX, HEADER=header,
        COLUMNS=columns, REFERENCE=REFERENCE, FORMAT=FORMAT)
