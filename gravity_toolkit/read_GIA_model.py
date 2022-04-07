#!/usr/bin/env python
u"""
read_GIA_model.py
Written by Tyler Sutterley (04/2022)

Reads GIA data files that can come in various formats depending on the group
Outputs spherical harmonics for the GIA rates and the GIA model parameters
Can also output fully normalized harmonics to netCDF4 or HDF5 formatted files

INPUTS:
    input_file: GIA file to read

OPTIONS:
    GIA: GIA model type to read and output
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    LMAX: maximum degree of spherical harmonics
    MMAX: maximum order of spherical harmonics
    DATAFORM: Spherical harmonic data output format
        None: output only as variables
        netCDF4: output to netCDF4 format (.nc)
        HDF5: output to HDF5 format (.H5)
    MODE: permissions mode of output spherical harmonic files

OUTPUTS:
    clm: cosine spherical harmonic of GIA rate
    slm: sine spherical harmonic of GIA rate
    l: spherical harmonic degree
    m: spherical harmonic order
    title: parameters of GIA model

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (https://www.h5py.org/)

PROGRAM DEPENDENCIES:
    ncdf_stokes.py: writes output spherical harmonic data to netcdf
    hdf5_stokes.py: writes output spherical harmonic data to HDF5
    ncdf_read_stokes.py: reads spherical harmonic data from netcdf
    hdf5_read_stokes.py: reads spherical harmonic data from HDF5

REFERENCES:
    E. R. Ivins, T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer, and
    K. M. Simon, "Antarctic contribution to sea level rise observed by
    GRACE with improved GIA correction." Journal of Geophysical Research:
    Solid Earth, 118(6), 3126-3141 (2013). https://doi.org/10.1002/jgrb.50208

    P. L. Whitehouse, M. J. Bentley, G. A. Milne, M. A. King, and I. D. Thomas,
    "A new glacial isostatic adjustment model for Antarctica: calibrated and
    tested using observations of relative sea-level change and present-day
    uplift rates." Geophysical Journal International, 190(3), 1464-1482 (2012).
    https://doi.org/10.1111/j.1365-246X.2012.05557.x

    M. J. R. Simpson, L. Wake, G. A. Milne, and P. Huybrechts, "The influence of
    decadal- to millennial-scale ice mass changes on present-day vertical land
    motion in Greenland: Implications for the interpretation of GPS
    observations." Journal of Geophysical Research: Solid Earth, 116(B2),
    B02406 (2011). https://doi.org/10.1029/2010JB007776

    W. R. Peltier, D. F. Argus, and R. Drummond, "Space geodesy constrains ice
    age terminal deglaciation: The global ICE‐6G_C (VM5a) model." Journal of
    Geophysical Research: Solid Earth, 120(1), 450-487 (2015).
    https://doi.org/10.1002/2014JB011176

    X. Wu, M. B. Heflin, H. Schotman, B. L. A. Vermeersen, D. Dong, R. S. Gross,
    E. R. Ivins, A. W. Moore, S. E. Owen, "Simultaneous estimation of global
    present-day water transport and glacial isostatic adjustment." Nature
    Geoscience, 3(9), 642-646 (2010). https://doi.org/10.1038/ngeo938

    G. A, J. Wahr, S. Zhong, "Computations of the viscoelastic response of a 3-D
    compressible Earth to surface loading: an application to Glacial Isostatic
    Adjustment in Antarctica and Canada." Geophysical Journal International,
    192(2), 557-572 (2013). https://doi.org/10.1093/gji/ggs030

    L. Caron, E. R. Ivins, E. Larour, S. Adhikari, J. Nilsson, and G. Blewitt,
    "GIA Model Statistics for GRACE Hydrology, Cryosphere, and Ocean Science."
    Geophysical Research Letters, 45(5), 2203-2212 (2018).
    https://doi.org/10.1002/2017GL076644

    W. R. Peltier, D. F. Argus, and R. Drummond, "Comment on 'An Assessment of
    the ICE-6G_C (VM5a) Glacial Isostatic Adjustment Model' by Purcell et al."
    Journal of Geophysical Research: Solid Earth, 123(2), 2019-2028 (2018).
    https://doi.org/10.1002/2016JB013844

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: use regular expressions to find ICE6G-D header positions
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 04/2020: include spherical harmonic degree and order in output dict
        added option to truncate to spherical harmonic order
    Updated 03/2020: updated for public release.  added reformatted ascii option
    Updated 08/2019: added ICE-6G Version D
    Updated 07/2019: added Geruo ICE-6G models and Caron JPL assimilation
    Updated 06/2018: using python3 compatible octal and input
    Updated 02/2017: added MODE to set output file permissions
    Updated 05-06/2016: using __future__ print function, use format for output
    Updated 08/2015: changed sys.exit to raise ValueError
    Updated 02/2015: update to reading the original file index index_orig
        and using regular expressions for some cases
    Updated 11/2014: added output option HDF5
    Updated 06/2014: changed message to sys.exit
    Updated 05/2014: added test IJ05 files
    Updated 09/2013: changed Wu parameter to 2010 (previously was in the name)
        this is in case the inversion is updated
    Updated 08/2013: updated comments (expanded)
        merged IJ05 with G13, W12a, SM09
        simplified ICE-6G code
        changed Wu code to convert to numerical array first
            then unwrapping the numerical array (vs. string)
    Updated 05/2013: converted to python and updated SM09 with latest best
        updated SM09 parameters to be automated from file input
    Updated 12/2012: changed the naming scheme for Simpson and Whitehouse
    Updated 09/2012: combined several GIA read programs into this standard
"""
from __future__ import print_function

import os
import re
import numpy as np
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.hdf5_stokes import hdf5_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.hdf5_read_stokes import hdf5_read_stokes

def read_GIA_model(input_file, GIA=None, LMAX=60, MMAX=None,
    DATAFORM=None, MODE=0o775):
    """
    Reads Glacial Isostatic Adjustment (GIA) data files

    Parameters
    ----------
    input_file: str
        full path to input GIA file
    GIA: str or NoneType, default None
        GIA model type to read and output

            - ``'IJ05-R2'``: Ivins R2 GIA Models [Ivins2013]_
            - ``'W12a'``: Whitehouse GIA Models [Whitehouse2012]_
            - ``'SM09'``: Simpson/Milne GIA Models [Simpson2009]_
            - ``'ICE6G'``: ICE-6G GIA Models [Peltier2015]_
            - ``'Wu10'``: Wu (2010) GIA Correction [Wu2010]_
            - ``'AW13-ICE6G'``: Geruo A ICE-6G GIA Models [A2013]_
            - ``'Caron'``: Caron JPL GIA Assimilation [Caron2018]_
            - ``'ICE6G-D'``: ICE-6G Version-D GIA Models [Peltier2018]_
            - ``'ascii'``: reformatted GIA in ascii format
            - ``'netCDF4'``: reformatted GIA in netCDF4 format
            - ``'HDF5'``: reformatted GIA in HDF5 format
    LMAX: int, default 60
        maximum degree of spherical harmonics
    MMAX: int or NoneType, default None
        maximum order of spherical harmonics
    DATAFORM: str or NoneType, default None
        Spherical harmonic data output format

            - ``None``: output only as variables
            - ``'netCDF4'``: output to netCDF4 format (.nc)
            - ``'HDF5'``: output to HDF5 format (.H5)
    MODE: oct, default 0o775
        Permissions mode of output spherical harmonic files

    Returns
    -------
    clm: float
        cosine spherical harmonic coefficients
    slm: float
        sine spherical harmonic coefficients
    l: int
        spherical harmonic degree
    m: int
        spherical harmonic order
    title: str
        parameters of GIA model

    References
    ----------
    .. [A2013] G. A, J. Wahr, S. Zhong,
        "Computations of the viscoelastic response of a 3-D
        compressible Earth to surface loading: an application to
        Glacial Isostatic Adjustment in Antarctica and Canada",
        *Geophysical Journal International*, 192(2), 557-572 (2013).
        `https://doi.org/10.1093/gji/ggs030 <https://doi.org/10.1093/gji/ggs030>`_
    .. [Caron2018] L. Caron, E. R. Ivins, E. Larour, S. Adhikari,
        J. Nilsson, and G. Blewitt, "GIA Model Statistics for GRACE Hydrology,
        Cryosphere, and Ocean Science", *Geophysical Research Letters*,
        45(5), 2203-2212 (2018).
        `https://doi.org/10.1002/2017GL076644 <https://doi.org/10.1002/2017GL076644>`_
    .. [Ivins2013] E. R. Ivins, T. S. James, J. Wahr, E. J. O. Schrama,
        F. W. Landerer, and K. M. Simon, "Antarctic contribution to
        sea level rise observed by GRACE with improved GIA correction",
        *Journal of Geophysical Research: Solid Earth*, 118(6), 3126-3141 (2013).
        `https://doi.org/10.1002/jgrb.50208 <https://doi.org/10.1002/jgrb.50208>`_
    .. [Simpson2009] M. J. R. Simpson, L. Wake, G. A. Milne, and P. Huybrechts,
        "The influence of decadal- to millennial-scale ice mass changes
        on present-day vertical land motion in Greenland: Implications
        for the interpretation of GPS observations",
        *Journal of Geophysical Research: Solid Earth*, 116(B2), B02406 (2011).
        `https://doi.org/10.1029/2010JB007776 <https://doi.org/10.1029/2010JB007776>`_
    .. [Peltier2015] W. R. Peltier, D. F. Argus, and R. Drummond,
        "Space geodesy constrains ice age terminal deglaciation:
        The global ICE‐6G_C (VM5a) model", *Journal of Geophysical Research:
        Solid Earth*, 120(1), 450-487 (2015).
        `https://doi.org/10.1002/2014JB011176 <https://doi.org/10.1002/2014JB011176>`_
    .. [Peltier2018] W. R. Peltier, D. F. Argus, and R. Drummond,
        "Comment on 'An Assessment of the ICE-6G_C (VM5a) Glacial
        Isostatic Adjustment Model' by Purcell et al.",
        *Journal of Geophysical Research: Solid Earth*, 123(2), 2019-2028 (2018).
        `https://doi.org/10.1002/2016JB013844 <https://doi.org/10.1002/2016JB013844>`_
    .. [Whitehouse2012] P. L. Whitehouse, M. J. Bentley, G. A. Milne,
        M. A. King, and I. D. Thomas, "A new glacial isostatic adjustment
        model for Antarctica: calibrated and tested using observations of
        relative sea-level change and present-day uplift rates",
        *Geophysical Journal International*, 190(3), 1464-1482 (2012).
        `https://doi.org/10.1111/j.1365-246X.2012.05557.x <https://doi.org/10.1111/j.1365-246X.2012.05557.x>`_
    .. [Wu2010] X. Wu, M. B. Heflin, H. Schotman, B. L. A. Vermeersen,
        D. Dong, R. S. Gross, E. R. Ivins, A. W. Moore, S. E. Owen,
        "Simultaneous estimation of global present-day water transport and
        glacial isostatic adjustment", *Nature Geoscience*, 3(9), 642-646 (2010).
        `https://doi.org/10.1038/ngeo938 <https://doi.org/10.1038/ngeo938>`_
    """

    #-- allocate for output Ylms
    #-- initially read for spherical harmonic degree up to LMAX
    #-- will truncate to MMAX before exiting program
    gia_Ylms = {}
    gia_Ylms['clm'] = np.zeros((LMAX+1,LMAX+1))
    gia_Ylms['slm'] = np.zeros((LMAX+1,LMAX+1))

    if (GIA == 'IJ05-R2'):#-- Ivins R2 IJ05 Models
        prefix = 'IJ05_R2'
        #-- regular expression file pattern
        file_pattern = r'Stokes.R2_(.*?)_L120'
    elif (GIA == 'W12a'):#-- Whitehouse W12a
        prefix = 'W12a'
        parameters = dict(B='Best', L='Lower', U='Upper')
        #-- regular expression file pattern
        file_pattern = r'grate_(B|L|U).clm'
    elif (GIA == 'SM09'):#-- Simpson Milne SM09
        prefix = 'SM09_Huy2'
        #-- regular expression file pattern
        file_pattern = r'grate_(\d+)p(\d)(\d+).clm'
    elif (GIA == 'ICE6G'):#-- ICE-6G VM5 viscosity profile
        prefix = 'ICE6G'
        #-- regular expression file pattern for test cases
        #file_pattern = r'Stokes_G_Rot_60_I6_A_(.*?)_L90'
        #-- regular expression file pattern for VM5
        file_pattern = r'Stokes_G_Rot_60_I6_A_(.*)'
    elif (GIA == 'Wu10'):#-- Wu Global Isnversion
        gia_Ylms['title'] = 'Wu_2010'
    elif (GIA == 'AW13-ICE6G'):#-- Geruo A ICE-6G model versions
        prefix = 'AW13'
        #-- regular expressions file pattern
        file_pattern = r'stokes\.(ice6g)[\.\_](.*?)(\.txt)?$'
    elif (GIA == 'Caron'): #-- Caron et al. (2018) expected GIA rate
        gia_Ylms['title'] = 'Caron_expt'
    elif (GIA == 'ICE6G-D'):#-- ICE-6G Version-D viscosity profile
        prefix = 'ICE6G-D'
        #-- regular expression file pattern for Version-D
        file_pattern = r'(ICE-6G_)?(.*?)[_]?Stokes_trend[_]?(.*?)\.txt$'

    #-- compile numerical expression operator
    rx = re.compile(r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?')

    #-- Header lines and scale factors for individual models
    if GIA in ('IJ05-R2','ICE6G'):
        #-- IJ05
        start = 0
        #-- scale factor for geodesy normalization
        scale = 1e-11
    elif (GIA == 'ICE6G-D'):
        #-- ICE-6G Version-D
        #-- scale factor for geodesy normalization
        scale = 1.0
    else:
        start = 0
        scale = 1.0

    #-- Reading GIA files (ICE-6G and Wu have more complex formats)
    if GIA in ('IJ05-R2','W12a','SM09','AW13-ICE6G'):
        #-- AW13, IJ05, W12a, SM09
        #-- AW13 notes: file headers
        #-- IJ05 notes: need to scale by 1e-11 for geodesy-normalization
        #-- exponents are denoted with D for double

        #-- opening gia data file and read contents
        with open(os.path.expanduser(input_file),'r') as f:
            gia_data = f.read().splitlines()
        #-- number of lines in file
        gia_lines = len(gia_data)

        #-- Skipping file header for geruo files with header
        for ii in range(start,gia_lines):
            #-- check if contents in line
            flag = bool(rx.search(gia_data[ii].replace('D','E')))
            if flag:
                #-- find numerical instances in line including exponents,
                #-- decimal points and negatives
                #-- Replacing Double Exponent with Standard Exponent
                line = rx.findall(gia_data[ii].replace('D','E'))
                l1 = np.int64(line[0])
                m1 = np.int64(line[1])
                #-- truncate to LMAX
                if (l1 <= LMAX) and (m1 <= LMAX):
                    #-- scaling to geodesy normalization
                    gia_Ylms['clm'][l1,m1] = np.float64(line[2])*scale
                    gia_Ylms['slm'][l1,m1] = np.float64(line[3])*scale

    elif (GIA == 'ICE6G'):
        #-- ICE-6G VM5 notes
        #-- need to scale by 1e-11 for geodesy-normalization
        #-- spherical harmonic degrees listed only on order 0
        #-- spherical harmonic order is not listed in file

        #-- opening gia data file and read contents
        with open(os.path.expanduser(input_file),'r') as f:
            gia_data = f.read().splitlines()

        #-- counter variable
        ii = 0
        for l in range(0, LMAX+1):
            for m in range(0, l+1):
                if ((m % 2) == 0):
                    #-- reading gia line if the order is even
                    #-- find numerical instances in line including exponents,
                    #-- decimal points and negatives
                    line = rx.findall(gia_data[ii])
                    #-- counter to next line
                    ii += 1
                    #-- if m is even: clm column = 1, slm column = 2
                    c = 0
                else: #-- if m is odd: clm column = 3, slm column = 4
                    c = 2
                if ((m == 0) or (m == 1)):
                    #-- l is column 1 if m == 0 or 1
                    #-- degree is not listed for other SHd: column 1 = clm
                    c += 1
                if (len(line) > 0):
                    #-- no empty lines
                    #-- convert to float and scale
                    gia_Ylms['clm'][l,m] = np.float64(line[0+c])*scale
                    gia_Ylms['slm'][l,m] = np.float64(line[1+c])*scale

    elif (GIA == 'Wu10'):
        #-- Wu (2010) notes:
        #-- Need to convert from mm geoid to fully normalized
        rad_e = 6.371e9#-- Average Radius of the Earth [mm]
        #-- The file starts with a header.
        #-- converting to numerical array (note 64 bit floating point)
        gia_data = np.loadtxt(os.path.expanduser(input_file),
            skiprows=1, dtype='f8')

        #-- counter variable to upwrap gia file
        ii = 0

        #-- Order of harmonics in the file:
        #--    1    0   c
        #--    1    1   c
        #--    1    1   s
        #--    2    0   c
        #--    2    1   c
        for l in range(1, LMAX+1):
            for m in range(0, l+1):
                for cs in range(0, 2):
                    #-- unwrapping GIA file and converting to geoid
                    #-- Clm
                    if (cs == 0):
                        gia_Ylms['clm'][l,m] = gia_data[ii]/rad_e
                        ii += 1
                    #-- Slm
                    if (m != 0) and (cs == 1):
                        gia_Ylms['slm'][l,m] = gia_data[ii]/rad_e
                        ii += 1

    elif (GIA == 'Caron'):
        #-- Caron et al. (2018)
        #-- The file starts with a header.
        #-- converting to numerical array (note 64 bit floating point)
        gia_data=np.loadtxt(os.path.expanduser(input_file),skiprows=4,
            dtype={'names':('l','m','Ylms'),'formats':('i','i','f8')})
        #-- Order of harmonics in the file
        #--    0    0   c
        #--    1    1   s
        #--    1    0   c
        #--    1    1   c
        #--    2    2   s
        #--    2    1   s
        #--    2    0   c
        #--    2    1   c
        #--    2    2   c
        for l,m,Ylm in zip(gia_data['l'],gia_data['m'],gia_data['Ylms']):
            #-- unwrapping GIA file
            if (m >= 0) and (l <= LMAX) and (m <= LMAX):#-- Clm
                gia_Ylms['clm'][l,m] = Ylm.copy()
            elif (m < 0) and (l <= LMAX) and (m <= LMAX):#-- Slm
                gia_Ylms['slm'][l,np.abs(m)] = Ylm.copy()

    #-- Reading ICE-6G Version-D  GIA files
    elif (GIA == 'ICE6G-D'):
        #-- opening gia data file and read contents
        with open(os.path.expanduser(input_file),'r') as f:
            gia_data = f.read().splitlines()
        #-- number of lines in file
        gia_lines = len(gia_data)

        #-- find header lines to skip
        h1 = r'^GRACE Approximation for degrees 0 to 2'
        h2 = r'^GRACE Approximation\/Absolute Sea-level Values for degrees \> 2'
        #-- header lines to skip
        header, = [(i+1) for i,l in enumerate(gia_data) if re.match(h1,l)]
        start, = [(i+1) for i,l in enumerate(gia_data) if re.match(h2,l)]

        #-- Calculating number of cos and sin harmonics to read from header
        n_harm = (2**2 + 3*2)//2 + 1
        #-- extract header for GRACE approximation
        for ii in range(header,header+n_harm):
            #-- check if contents in line
            flag = bool(rx.search(gia_data[ii].replace('D','E')))
            if flag:
                #-- find numerical instances in line including exponents,
                #-- decimal points and negatives
                #-- Replacing Double Exponent with Standard Exponent
                line = rx.findall(gia_data[ii].replace('D','E'))
                l1 = np.int64(line[0])
                m1 = np.int64(line[1])
                #-- truncate to LMAX
                if (l1 <= LMAX) and (m1 <= LMAX):
                    #-- scaling to geodesy normalization
                    gia_Ylms['clm'][l1,m1] = np.float64(line[2])*scale
                    gia_Ylms['slm'][l1,m1] = np.float64(line[3])*scale

        #-- Skipping rest of file header
        for ii in range(start,gia_lines):
            #-- check if contents in line
            flag = bool(rx.search(gia_data[ii].replace('D','E')))
            if flag:
                #-- find numerical instances in line including exponents,
                #-- decimal points and negatives
                #-- Replacing Double Exponent with Standard Exponent
                line = rx.findall(gia_data[ii].replace('D','E'))
                l1 = np.int64(line[0])
                m1 = np.int64(line[1])
                #-- truncate to LMAX
                if (l1 <= LMAX) and (m1 <= LMAX):
                    #-- scaling to geodesy normalization
                    gia_Ylms['clm'][l1,m1] = np.float64(line[2])*scale
                    gia_Ylms['slm'][l1,m1] = np.float64(line[3])*scale

    elif (GIA == 'ascii'):
        #-- reading GIA data from reformatted (simplified) ascii files
        dtype = {'names':('l','m','clm','slm'), 'formats':('i','i','f8','f8')}
        Ylms = np.loadtxt(os.path.expanduser(input_file), dtype=dtype)
        for ll,mm,clm,slm in zip(Ylms['l'],Ylms['m'],Ylms['clm'],Ylms['slm']):
            #-- only using coefficients within the spherical harmonic range
            if ((ll <= LMAX) and (mm <= LMAX)):
                gia_Ylms['clm'][ll,mm] = clm.copy()
                gia_Ylms['slm'][ll,mm] = slm.copy()
        #-- copy filename (without extension) for parameters
        gia_Ylms['title'] = os.path.basename(os.path.splitext(input_file)[0])

    elif (GIA == 'netCDF4'):
        #-- reading GIA data from reformatted netCDF4 files
        Ylms = ncdf_read_stokes(os.path.expanduser(input_file),DATE=False)
        #-- truncate to degree and order LMAX
        gia_Ylms['clm'][:,:] = Ylms['clm'][:LMAX+1,:LMAX+1]
        gia_Ylms['slm'][:,:] = Ylms['slm'][:LMAX+1,:LMAX+1]
        #-- copy title for parameters
        gia_Ylms['title'] = Ylms['attributes']['title']

    elif (GIA == 'HDF5'):
        #-- reading GIA data from reformatted HDF5 files
        Ylms = hdf5_read_stokes(os.path.expanduser(input_file),DATE=False)
        #-- truncate to degree and order LMAX
        gia_Ylms['clm'][:,:] = Ylms['clm'][:LMAX+1,:LMAX+1]
        gia_Ylms['slm'][:,:] = Ylms['slm'][:LMAX+1,:LMAX+1]
        #-- copy title for parameters
        gia_Ylms['title'] = Ylms['attributes']['title']

    #-- extract rheology from the file name
    if GIA in ('IJ05-R2','ICE6G'):
        #-- for IJ05 and ICE-6G models:
        #-- adding file specific earth parameters
        parameters, = re.findall(file_pattern,os.path.basename(input_file))
        gia_Ylms['title'] = '{0}_{1}'.format(prefix,parameters)
    elif (GIA == 'W12a'):
        #-- for Whitehouse W12a (BEST, LOWER, UPPER):
        model = re.findall(file_pattern,os.path.basename(input_file)).pop()
        gia_Ylms['title'] = '{0}_{1}'.format(prefix,parameters[model])
    elif (GIA == 'SM09'):
        #-- for SM09:
        #-- making parameters in the file similar to IJ05
        #-- split rheological parameters between lithospheric thickness,
        #-- upper mantle viscosity and lower mantle viscosity
        LTh,UMV,LMV=re.findall(file_pattern,os.path.basename(input_file)).pop()
        #-- formatting rheology parameters similar to IJ05 models
        gia_Ylms['title'] = '{0}_{1}_.{2}_{3}'.format(prefix,LTh,UMV,LMV)
    elif (GIA == 'ICE6G-D'):
        #-- for ICE-6G Version-D models:
        #-- adding file specific earth parameters
        m1,p1,p2 = re.findall(file_pattern,os.path.basename(input_file)).pop()
        gia_Ylms['title'] = '{0}_{1}{2}'.format(prefix,p1,p2)
    elif (GIA == 'AW13-ICE6G'):
        #-- for Geruo ICE-6G cases
        #-- extract the ice history and case flags
        hist,case,sf=re.findall(file_pattern,os.path.basename(input_file)).pop()
        gia_Ylms['title'] = '{0}_{1}_{2}'.format(prefix,hist,case)

    #-- output spherical harmonic degree and order
    gia_Ylms['l'],gia_Ylms['m'] = (np.arange(LMAX+1),np.arange(LMAX+1))
    #-- output harmonics to netCDF4 or HDF5 file
    if (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        output_file = 'stokes_{0}_L{1:d}.nc'.format(gia_Ylms['title'],LMAX)
        ncdf_stokes(gia_Ylms['clm'],gia_Ylms['slm'],gia_Ylms['l'],gia_Ylms['m'],
            0,0,FILENAME=os.path.join(os.path.dirname(input_file),output_file),
            TITLE=gia_Ylms['title'], VERBOSE=False, DATE=False)
        #-- set permissions level of output file
        os.chmod(os.path.join(os.path.dirname(input_file),output_file), MODE)
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        output_file = 'stokes_{0}_L{1:d}.H5'.format(gia_Ylms['title'],LMAX)
        hdf5_stokes(gia_Ylms['clm'],gia_Ylms['slm'],gia_Ylms['l'],gia_Ylms['m'],
            0,0,FILENAME=os.path.join(os.path.dirname(input_file),output_file),
            TITLE=gia_Ylms['title'], VERBOSE=False, DATE=False)
        #-- set permissions level of output file
        os.chmod(os.path.join(os.path.dirname(input_file),output_file), MODE)

    #-- truncate to MMAX if specified
    if MMAX is not None:
        #-- spherical harmonic variables
        gia_Ylms['clm'] = gia_Ylms['clm'][:,:MMAX+1]
        gia_Ylms['slm'] = gia_Ylms['slm'][:,:MMAX+1]
        #-- spherical harmonic order
        gia_Ylms['m'] = gia_Ylms['m'][:MMAX+1]

    #-- return the harmonics and the parameters
    return gia_Ylms
