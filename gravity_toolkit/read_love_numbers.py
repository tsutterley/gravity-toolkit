#!/usr/bin/env python
u"""
read_love_numbers.py
Written by Tyler Sutterley (11/2024)

Reads sets of load Love numbers from PREM and applies isomorphic parameters
Linearly interpolates load Love/Shida numbers for missing degrees
Linearly extrapolates load Love/Shida numbers beyond maximum degree of dataset

INPUTS:
    love_numbers_file: Elastic load Love/Shida numbers file
        computed using Preliminary Reference Earth Model (PREM) outputs

OUTPUTS:
    hl: Love number of Vertical Displacement
    kl: Love number of Gravitational Potential
    ll: Love (Shida) number of Horizontal Displacement

OPTIONS:
    LMAX: truncate or interpolate to maximum spherical harmonic degree
    HEADER: number of header lines to be skipped
    COLUMNS: column names of ascii file
        l: spherical harmonic degree
        hl: vertical displacement
        kl: gravitational potential
        ll: horizontal displacement
    REFERENCE: Reference frame for calculating degree 1 Love/Shida numbers
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
    G. Blewitt, "Self-consistency in reference frames, geocenter definition,
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
    Updated 11/2024: allow reading where degree is infinite
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 08/2023: add string representation of the love_numbers class
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: improve typing for variables in docstrings
    Updated 02/2023: fix degree zero case and add load love number formatter
        added options for hard and soft PREM sediment cases
        added data class for load love numbers with attributes for model
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: use logging for debugging level verbose output
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
import io
import re
import logging
import pathlib
import numpy as np
from gravity_toolkit.utilities import get_data_path

# default maximum degree and order in case of infinite
_default_max_degree = 100000

# PURPOSE: read load Love/Shida numbers from PREM
def read_love_numbers(love_numbers_file, LMAX=None, HEADER=2,
    COLUMNS=['l','hl','kl','ll'], REFERENCE='CE', FORMAT='tuple'):
    """
    Reads PREM load Love/Shida numbers file and applies isomorphic
    parameters :cite:p:`Dziewonski:1981bz,Blewitt:2003bz`
    :cite:p:`Wahr:1998hy`

    Parameters
    ----------
    love_numbers_file: str
        Elastic load Love/Shida numbers file
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
        Reference frame of degree 1 Love/Shida numbers

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
            - ``'class'``: ``love_numbers`` class

    Returns
    -------
    hl: np.ndarray
        Love number of Vertical Displacement
    kl: np.ndarray
        Love number of Gravitational Potential
    ll: np.ndarray
        Love (Shida) number of Horizontal Displacement
    """
    # Input load Love/Shida number data file and read contents
    file_contents = extract_love_numbers(love_numbers_file)

    # compile regular expression operator to find numerical instances
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?)|(?:inf))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    # extract maximum spherical harmonic degree from final line in file
    degree = rx.findall(file_contents[-1])[COLUMNS.index('l')]
    if LMAX is None and (degree == 'inf'):
        LMAX = np.copy(_default_max_degree)
    elif LMAX is None:
        LMAX = int(degree)

    # dictionary of output Love/Shida numbers
    love = {}
    # spherical harmonic degree
    love['l'] = np.arange(LMAX+1)
    # vertical displacement hl
    # gravitational potential kl
    # horizontal displacement ll (Shida number)
    for n in ('hl','kl','ll'):
        love[n] = np.zeros((LMAX+1))
    # check if needing to interpolate between degrees
    flag = np.ones((LMAX+1),dtype=bool)
    # for each line in the file (skipping header lines)
    for file_line in file_contents[HEADER:]:
        # find numerical instances in line
        # replacing fortran double precision exponential
        love_numbers = rx.findall(file_line.replace('D','E'))
        # spherical harmonic degree
        degree = love_numbers[COLUMNS.index('l')]
        l = _default_max_degree if (degree == 'inf') else int(degree)
        # truncate to spherical harmonic degree LMAX
        if (l <= LMAX):
            # convert Love/Shida numbers to float
            # vertical displacement hl
            # gravitational potential kl
            # horizontal displacement ll (Shida number)
            for n in ('hl','kl','ll'):
                love[n][l] = np.float64(love_numbers[COLUMNS.index(n)])
            # set interpolation flag for degree
            flag[l] = False

    # return Love/Shida numbers in output format
    if (LMAX == 0):
        return love_number_formatter(love, FORMAT=FORMAT)

    # if needing to linearly interpolate Love/Shida numbers
    if np.any(flag):
        # linearly interpolate following Wahr (1998)
        for n in ('hl','kl','ll'):
            love[n][flag] = np.interp(love['l'][flag],
                love['l'][~flag], love[n][~flag])

    # if needing to linearly extrapolate Love/Shida numbers
    # NOTE: use caution if extrapolating far beyond the
    # maximum degree of the Love/Shida numbers dataset
    for lint in range(l,LMAX+1):
        # linearly extrapolate to maximum degree
        for n in ('hl','kl','ll'):
            love[n][lint] = 2.0*love[n][lint-1] - love[n][lint-2]

    # calculate isomorphic parameters for different reference frames
    # From Blewitt (2003), Wahr (1998), Trupin (1992) and Farrell (1972)
    if (REFERENCE.upper() == 'CF'):
        # Center of Surface Figure
        alpha = (love['hl'][1] + 2.0*love['ll'][1])/3.0
    elif (REFERENCE.upper() == 'CL'):
        # Center of Surface Lateral Figure
        alpha = love['ll'][1].copy()
    elif (REFERENCE.upper() == 'CH'):
        # Center of Surface Height Figure
        alpha = love['hl'][1].copy()
    elif (REFERENCE.upper() == 'CM'):
        # Center of Mass of Earth System
        alpha = 1.0
    elif (REFERENCE.upper() == 'CE'):
        # Center of Mass of Solid Earth
        alpha = 0.0
    else:
        raise Exception(f'Invalid Reference Frame {REFERENCE}')
    # apply isomorphic parameters
    for n in ('hl','kl','ll'):
        love[n][1] -= alpha

    # return Love/Shida numbers in output format
    return love_number_formatter(love, FORMAT=FORMAT)

# PURPOSE: return load Love/Shida numbers in a particular format
def love_number_formatter(love, FORMAT='tuple'):
    """
    Converts a dictionary of Load Love/Shida Numbers
    to a particular output format

    Parameters
    ----------
    love: dict
        Load Love/Shida numbers
    FORMAT: str, default 'tuple'
        Format of output variables

            - ``'dict'``: dictionary with variable keys
            - ``'tuple'``: tuple with variable order (``hl``, ``kl``, ``ll``)
            - ``'zip'``: aggregated variable sets
            - ``'class'``: ``love_numbers`` class

    Returns
    -------
    hl: np.ndarray
        Love number of Vertical Displacement
    kl: np.ndarray
        Love number of Gravitational Potential
    ll: np.ndarray
        Love (Shida) number of Horizontal Displacement
    """
    if (FORMAT == 'dict'):
        return love
    elif (FORMAT == 'tuple'):
        return (love['hl'], love['kl'], love['ll'])
    elif (FORMAT == 'zip'):
        return zip(love['hl'], love['kl'], love['ll'])
    elif (FORMAT == 'class'):
        return love_numbers().from_dict(love)

# PURPOSE: read input file and extract contents
def extract_love_numbers(love_numbers_file):
    """
    Read load Love/Shida number file and extract contents

    Parameters
    ----------
    love_numbers_file: str, bytesIO or pathlib.Path
        Elastic load Love/Shida numbers file
    """
    # check if input Love/Shida numbers are a string or bytesIO object
    if isinstance(love_numbers_file, (str, pathlib.Path)):
        # tilde expansion of load love number data file
        love_numbers_file = pathlib.Path(love_numbers_file).expanduser().absolute()
        # check that load Love/Shida number data file is present in file system
        if not love_numbers_file.exists():
            raise FileNotFoundError(f'{str(love_numbers_file)} not found')
        # Input load Love/Shida number data file and read contents
        with love_numbers_file.open(mode='r', encoding='utf8') as f:
            return f.read().splitlines()
    elif isinstance(love_numbers_file, io.IOBase):
        # read contents from load Love/Shida number data
        return love_numbers_file.read().decode('utf8').splitlines()
    else:
        raise ValueError('Invalid Love/Shida numbers file input')

# PURPOSE: read load Love/Shida numbers for a range of spherical harmonic degrees
def load_love_numbers(LMAX, LOVE_NUMBERS=0, REFERENCE='CF', FORMAT='tuple'):
    """
    Wrapper function for reading PREM load Love/Shida numbers for a
    range of spherical harmonic degrees and applying
    isomorphic parameters :cite:p:`Blewitt:2003bz`

    Parameters
    ----------
    LMAX: int
        maximum spherical harmonic degree
    LOVE_NUMBERS: int, default 0
        Treatment of the Load Love/Shida numbers

            - ``0``: :cite:p:`Han:1995go` values from PREM
            - ``1``: :cite:p:`Gegout:2010gc` values from PREM
            - ``2``: :cite:p:`Wang:2012gc` values from PREM
            - ``3``: :cite:p:`Wang:2012gc` values from PREM with hard sediment
            - ``4``: :cite:p:`Wang:2012gc` values from PREM with soft sediment
    REFERENCE: str
        Reference frame for calculating degree 1 Love/Shida numbers :cite:p:`Blewitt:2003bz`

            - ``'CF'``: Center of Surface Figure (default)
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    FORMAT: str, default 'tuple'
        Format of output variables

            - ``'dict'``: dictionary with variable keys as listed above
            - ``'tuple'``: tuple with variable order (``hl``, ``kl``, ``ll``)
            - ``'zip'``: aggregated variable sets
            - ``'class'``: ``love_numbers`` class

    Returns
    -------
    hl: np.ndarray
        Love number of Vertical Displacement
    kl: np.ndarray
        Love number of Gravitational Potential
    ll: np.ndarray
        Love (Shida) number of Horizontal Displacement
    """
    # load Love/Shida numbers file
    if (LOVE_NUMBERS == 0):
        # PREM outputs from Han and Wahr (1995)
        # https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
        love_numbers_file = get_data_path(['data','love_numbers'])
        model = 'PREM'
        citation = 'Han and Wahr (1995)'
        header = 2
        columns = ['l','hl','kl','ll']
    elif (LOVE_NUMBERS == 1):
        # PREM outputs from Gegout (2005)
        # http://gemini.gsfc.nasa.gov/aplo/
        love_numbers_file = get_data_path(['data','Load_Love2_CE.dat'])
        model = 'PREM'
        citation = 'Gegout et al. (2010)'
        header = 3
        columns = ['l','hl','ll','kl']
    elif (LOVE_NUMBERS == 2):
        # PREM outputs from Wang et al. (2012)
        # https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREM-LLNs-truncated.dat'])
        model = 'PREM'
        citation = 'Wang et al. (2012)'
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    elif (LOVE_NUMBERS == 3):
        # PREM hard outputs from Wang et al. (2012)
        # case with 0.46 kilometers thick hard sediment
        # https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREMhard-LLNs-truncated.dat'])
        model = 'PREMhard'
        citation = 'Wang et al. (2012)'
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    elif (LOVE_NUMBERS == 4):
        # PREM soft outputs from Wang et al. (2012)
        # case with 0.52 kilometers thick soft sediment
        # https://doi.org/10.1016/j.cageo.2012.06.022
        love_numbers_file = get_data_path(['data','PREMsoft-LLNs-truncated.dat'])
        model = 'PREMsoft'
        citation = 'Wang et al. (2012)'
        header = 1
        columns = ['l','hl','ll','kl','nl','nk']
    else:
        raise ValueError(f'Unknown Love Numbers Type {LOVE_NUMBERS:d}')
    # validate as pathlib object
    love_numbers_file = pathlib.Path(love_numbers_file).expanduser().absolute()
    # log load Love/Shida numbers file if debugging
    logging.debug(f'Reading Love/Shida numbers file: {str(love_numbers_file)}')
    # LMAX of load Love/Shida numbers from Han and Wahr (1995) is 696.
    # from Wahr (2007) linearly interpolating kl works
    # however, as we are linearly extrapolating out, do not make
    # LMAX too much larger than 696
    # read arrays of kl, hl, and ll Love/Shida Numbers
    love = read_love_numbers(love_numbers_file, LMAX=LMAX, HEADER=header,
        COLUMNS=columns, REFERENCE=REFERENCE, FORMAT=FORMAT)
    # append model and filename attributes to class
    if (FORMAT == 'class'):
        love.filename = love_numbers_file.name
        love.reference = REFERENCE
        love.model = model
        love.citation = citation
    # return the load love numbers
    return love

class love_numbers(object):
    """
    Data class for Load Love/Shida numbers

    Attributes
    ----------
    lmax: int
        maximum degree of the Load Love/Shida Numbers
    l: np.ndarray
        Spherical harmonic degrees
    hl: np.ndarray or List
        Love number of Vertical Displacement
    kl: np.ndarray or List
        Love number of Gravitational Potential
    ll: np.ndarray or List
        Love (Shida) number of Horizontal Displacement
    reference: str
        Reference frame for degree 1 Love/Shida numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth

    model: str
        Reference Earth Model
    citation: str
        Citation for Reference Earth Model
    filename: str
        input filename of Load Love/Shida Numbers
    """
    np.seterr(invalid='ignore')
    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault('lmax',None)
        # set default class attributes
        self.hl=[]
        self.kl=[]
        self.ll=[]
        self.lmax=kwargs['lmax']
        # calculate spherical harmonic degree (0 is falsy)
        self.l=np.arange(self.lmax+1) if (self.lmax is not None) else None
        self.reference=None
        self.model=None
        self.citation=None
        self.filename=None

    def from_dict(self, d):
        """
        Convert a dict object to a ``love_numbers`` object

        Parameters
        ----------
        d: dict
            dictionary object to be converted
        """
        # retrieve each Load Love/Shida Number
        for key in ('hl','kl','ll'):
            setattr(self, key, d.get(key))
        self.lmax = len(self.hl) - 1
        # calculate spherical harmonic degree
        self.update_dimensions()
        return self

    def to_dict(self):
        """
        Convert a ``love_numbers`` object to a dict object

        Returns
        -------
        d: dict
            output dictionary object
        """
        # retrieve each Load Love/Shida Number
        d = {}
        for key in ('hl','kl','ll'):
            d[key] = getattr(self, key)
        return d

    def to_tuple(self):
        """
        Convert a ``love_numbers`` object to a tuple object

        Returns
        -------
        t: tuple
            output tuple object
        """
        # return Load Love/Shida Numbers
        return (self.hl, self.kl, self.ll)

    def transform(self, reference):
        """
        Calculate and apply calculate isomorphic parameters to
        transform from the Center of Mass of the Solid Earth
        Reference Frame :cite:p:`Blewitt:2003bz`

        Parameters
        ----------
        reference: str
            Output reference frame for degree 1 Love/Shida numbers

                - ``'CF'``: Center of Surface Figure
                - ``'CL'``: Center of Surface Lateral Figure
                - ``'CH'``: Center of Surface Height Figure
                - ``'CM'``: Center of Mass of Earth System
                - ``'CE'``: Center of Mass of Solid Earth
        """
        # calculate isomorphic parameters for different reference frames
        # From Blewitt (2003), Wahr (1998), Trupin (1992) and Farrell (1972)
        if (reference.upper() == 'CF'):
            # Center of Surface Figure
            alpha = (self.hl[1] + 2.0*self.ll[1])/3.0
        elif (reference.upper() == 'CL'):
            # Center of Surface Lateral Figure
            alpha = self.ll[1].copy()
        elif (reference.upper() == 'CH'):
            # Center of Surface Height Figure
            alpha = self.hl[1].copy()
        elif (reference.upper() == 'CM'):
            # Center of Mass of Earth System
            alpha = 1.0
        elif (reference.upper() == 'CE'):
            # Center of Mass of Solid Earth
            alpha = 0.0
        else:
            raise Exception(f'Invalid Reference Frame {reference}')
        # apply isomorphic parameters
        self.hl[1] -= alpha
        self.kl[1] -= alpha
        self.ll[1] -= alpha
        # set reference attribute
        self.reference = reference
        return self

    def update_dimensions(self):
        """
        Update the dimensions of the ``love_numbers`` object
        """
        # calculate spherical harmonic degree (0 is falsy)
        self.l=np.arange(self.lmax+1) if (self.lmax is not None) else None
        return self

    def __str__(self):
        """String representation of the ``love_numbers`` object
        """
        properties = ['gravity_toolkit.love_numbers']
        properties.append(f"    citation: {self.citation}")
        properties.append(f"    earth_model: {self.model}")
        properties.append(f"    max_degree: {self.lmax}")
        properties.append(f"    reference: {self.reference}")
        return '\n'.join(properties)

    def __len__(self):
        """Number of degrees
        """
        return len(self.l)

    def __iter__(self):
        """Iterate over load Love/Shida numbers variables
        """
        yield self.hl
        yield self.kl
        yield self.ll

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
