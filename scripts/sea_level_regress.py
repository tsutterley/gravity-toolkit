#!/usr/bin/env python
u"""
sea_level_regress.py
Written by Tyler Sutterley (06/2023)

Reads in sea level grid files and calculates the trends at
    each grid point following an input regression model

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for mascon files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month for time series regression
    -E X, --end X: ending GRACE/GRACE-FO month for time series regression
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -l X, --lmax X: maximum spherical harmonic degree
    -F X, --format X: input data format for auxiliary files
        ascii
        netCDF4
        HDF5
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        AW13-IJ05: Geruo A IJ05-R2 GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --order X: regression fit polynomial order
    --cycles X: regression fit cyclical terms
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    time_series.regress.py: calculates trend coefficients using least-squares
    time_series.amplitude.py: calculates the amplitude and phase of a harmonic
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 06/2023: append amplitude and phase titles when creating flags
        more tidal aliasing periods using values from Ray and Luthcke (2006)
    Updated 05/2023: split S2 tidal aliasing terms into GRACE and GRACE-FO eras
        output data and error variables into single files
        use fit module for getting tidal aliasing terms
        use pathlib to define and operate on paths
    Updated 03/2023: updated inputs to spatial from_ascii function
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 10/2020: use argparse to set command line parameters
    Updated 06/2020: using spatial data class for input and output operations
    Updated 01/2020: output seasonal amplitude and phase
    Updated 10/2019: changing Y/N flags to True/False
    Updated 01/2019: include 161-day S2 tidal aliasing terms in regression
    Updated 12/2018: added parallel processing with multiprocessing
    Updated 11/2018: using future division for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
        can run for different sets of months using the --start and --end options
        use output_data wrapper function for writing data to file
    Updated 04/2018: changed the names of the output log files
    Updated 03/2018: added option --mode to set the output file permissions
    Updated 02/2018: added parameter EXPANSION to specify the SLF truncation
    Updated 08/2017: running at 0.5x0.5 degrees
    Written 08/2017
"""
from __future__ import print_function, division

import sys
import os
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# program module to run with specified parameters
def sea_level_regress(PROC, DREL, DSET, LMAX,
    START=None,
    END=None,
    MISSING=None,
    GIA=None,
    GIA_FILE=None,
    DATAFORM=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    ORDER=None,
    CYCLES=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    VERBOSE=0,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # GRACE/GRACE-FO months
    months = sorted(set(np.arange(START,END+1)) - set(MISSING))
    nmon = len(months)
    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # input GIA stokes coefficients to get titles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    gia_str = f'_{GIA_Ylms_rate.title.upper()}' if GIA else ''
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''

    # input and output file formats
    input_format = 'SLF_ITERATION_{0}{1}{2}{3}_L{4:d}_{5:03d}.{6}'
    output_format = 'SLF_ITERATION_{0}{1}{2}{3}_L{4:d}_{5}{6}_{7:03d}-{8:03d}.{9}'

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF4 file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # create land function
    land_function = np.zeros((nlat, nlon),dtype=np.float64)
    # combine land and island levels for land function
    indy,indx = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indy,indx] = 1.0
    # calculate ocean function from land function
    ocean_function = 1.0 - land_function
    # bad value and land function as mask
    FILL_VALUE = -9999.0
    MASK = land_function.astype(bool)

    # input data spatial object
    spatial_list = []
    for t in range(0, nmon):
        # sea level file for month
        fi = input_format.format(ITERATION,dset_str,gia_str,ocean_str,
            EXPANSION,months[t],suffix)
        INPUT_FILE = OUTPUT_DIRECTORY.joinpath(fi)
        # read sea level file
        if (DATAFORM == 'ascii'):
            dinput = gravtk.spatial().from_ascii(INPUT_FILE,
                spacing=[dlon,dlat], nlon=nlon, nlat=nlat)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            dinput = gravtk.spatial().from_netCDF4(INPUT_FILE)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            dinput = gravtk.spatial().from_HDF5(INPUT_FILE)
        # append to spatial list
        dinput.replace_invalid(FILL_VALUE,mask=MASK)
        spatial_list.append(dinput)

    # concatenate list to single spatial object
    grid = gravtk.spatial().from_list(spatial_list)
    spatial_list = None

    # Setting output parameters for each fit type
    coef_str = ['x{0:d}'.format(o) for o in range(ORDER+1)]
    unit_suffix = [' yr^{0:d}'.format(-o) if o else '' for o in range(ORDER+1)]
    if (ORDER == 0):# Mean
        fit_longname = ['Mean']
    elif (ORDER == 1):# Trend
        fit_longname = ['Constant','Trend']
    elif (ORDER == 2):# Quadratic
        fit_longname = ['Constant','Linear','Quadratic']

    # amplitude string for cyclical components
    amp_str = []
    # amplitude and phase titles for cyclical components
    amp_title = {}
    ph_title = {}
    # unique tidal aliasing periods from Ray and Luthcke (2006)
    tidal_aliasing = {}
    tidal_aliasing['Q1'] = 9.1
    tidal_aliasing['O1'] = 13.6
    tidal_aliasing['P1'] = 171.2
    tidal_aliasing['S1'] = 322.1
    tidal_aliasing['K1'] = 2725.4
    tidal_aliasing['J1'] = 27.8
    tidal_aliasing['M2'] = 13.5
    tidal_aliasing['S2'] = 161.0
    tidal_aliasing['K2'] = 1362.7
    # extra terms for tidal aliasing components or custom fits
    TERMS = []
    term_index = []
    for i,c in enumerate(CYCLES):
        # check if fitting with semi-annual or annual terms
        if (c == 0.5):
            coef_str.extend(['SS','SC'])
            amp_str.append('SEMI')
            amp_title['SEMI'] = 'Semi-Annual Amplitude'
            ph_title['SEMI'] = 'Semi-Annual Phase'
            fit_longname.extend(['Semi-Annual Sine', 'Semi-Annual Cosine'])
            unit_suffix.extend(['',''])
        elif (c == 1.0):
            coef_str.extend(['AS','AC'])
            amp_str.append('ANN')
            amp_title['ANN'] = 'Annual Amplitude'
            ph_title['ANN'] = 'Annual Phase'
            fit_longname.extend(['Annual Sine', 'Annual Cosine'])
            unit_suffix.extend(['',''])
        # check if fitting with tidal aliasing terms
        for t,period in tidal_aliasing.items():
            if np.isclose(c, (period/365.25)):
                # terms for tidal aliasing during GRACE and GRACE-FO periods
                TERMS.extend(gravtk.time_series.aliasing_terms(grid.time,
                    period=period))
                # labels for tidal aliasing during GRACE period
                coef_str.extend([f'{t}SGRC', f'{t}CGRC'])
                amp_str.append(f'{t}GRC')
                amp_title[f'{t}GRC'] = f'{t} Tidal Alias (GRACE) Amplitude'
                ph_title[f'{t}GRC'] = f'{t} Tidal Alias (GRACE) Phase'
                fit_longname.append(f'{t} Tidal Alias (GRACE) Sine')
                fit_longname.append(f'{t} Tidal Alias (GRACE) Cosine')
                unit_suffix.extend(['',''])
                # labels for tidal aliasing during GRACE-FO period
                coef_str.extend([f'{t}SGFO', f'{t}CGFO'])
                amp_str.append(f'{t}GFO')
                amp_title[f'{t}GFO'] = f'{t} Tidal Alias (GRACE-FO) Amplitude'
                ph_title[f'{t}GFO'] = f'{t} Tidal Alias (GRACE-FO) Phase'
                fit_longname.append(f'{t} Tidal Alias (GRACE-FO) Sine')
                fit_longname.append(f'{t} Tidal Alias (GRACE-FO) Cosine')
                unit_suffix.extend(['',''])
                # index to remove the original tidal aliasing term
                term_index.append(i)
    # remove the original tidal aliasing terms
    CYCLES = np.delete(CYCLES, term_index)

    # Fitting seasonal components
    ncomp = len(coef_str)
    ncycles = 2*len(CYCLES) + len(TERMS)
    # confidence interval for regression fit errors
    CONF = 0.95

    # Allocating memory for output variables
    out = dinput.zeros_like()
    out.data = np.zeros((nlat, nlon, ncomp))
    out.error = np.zeros((nlat, nlon, ncomp))
    out.mask = np.ones((nlat, nlon, ncomp),dtype=bool)
    # Fit Significance
    FS = {}
    # SSE: Sum of Squares Error
    # AIC: Akaike information criterion
    # BIC: Bayesian information criterion
    # R2Adj: Adjusted Coefficient of Determination
    for key in ['SSE','AIC','BIC','R2Adj']:
        FS[key] = dinput.zeros_like()

    # valid values for ocean function
    indy,indx = np.nonzero(ocean_function)
    # calculate the regression coefficients and fit significance
    for i,j in zip(indy,indx):
        # Calculating the regression coefficients
        tsbeta = gravtk.time_series.regress(grid.time, grid.data[i,j,:],
            ORDER=ORDER, CYCLES=CYCLES, TERMS=TERMS, CONF=CONF)
        # save regression components
        for k in range(0, ncomp):
            out.data[i,j,k] = tsbeta['beta'][k]
            out.error[i,j,k] = tsbeta['error'][k]
            out.mask[i,j,k] = False
        # Fit significance terms
        # Degrees of Freedom
        nu = tsbeta['DOF']
        # Converting Mean Square Error to Sum of Squares Error
        FS['SSE'].data[i,j] = tsbeta['MSE']*nu
        FS['AIC'].data[i,j] = tsbeta['AIC']
        FS['BIC'].data[i,j] = tsbeta['BIC']
        FS['R2Adj'].data[i,j] = tsbeta['R2Adj']

    # list of output files
    output_files = []
    # Output spatial files
    for i in range(0,ncomp):
        # output spatial file name
        f1 = (ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            coef_str[i], '', START, END, suffix)
        file1 = OUTPUT_DIRECTORY.joinpath(output_format.format(*f1))
        # full attributes
        UNITS_TITLE = f'centimeters{unit_suffix[i]}'
        LONGNAME = 'Equivalent_Water_Thickness'
        FILE_TITLE = f'Sea_Level_Fingerprint_{fit_longname[i]}'
        # output regression fit and fit error to file
        output = out.index(i, date=False)
        output_data(output, FILENAME=file1, DATAFORM=DATAFORM,
            UNITS=UNITS_TITLE, LONGNAME=LONGNAME, TITLE=FILE_TITLE,
            CONF=CONF, VERBOSE=VERBOSE, MODE=MODE)
        # add output files to list object
        output_files.append(file1)

    # if fitting coefficients with seasonal components
    # output amplitude and phase of cyclical components
    for i,flag in enumerate(amp_str):
        # Indice pointing to the cyclical components
        j = 1 + ORDER + 2*i
        # Allocating memory for output amplitude and phase
        amp = dinput.zeros_like()
        amp.error = np.zeros((nlat, nlon))
        ph = dinput.zeros_like()
        ph.error = np.zeros((nlat, nlon))
        # calculating amplitude and phase of spatial field
        amp.data[indy,indx],ph.data[indy,indx] = gravtk.time_series.amplitude(
            out.data[indy,indx,j], out.data[indy,indx,j+1]
        )
        # convert phase from -180:180 to 0:360
        ii,jj = np.nonzero((ph.data < 0) & np.logical_not(ph.mask))
        ph.data[ii,jj] += 360.0
        # Amplitude Error
        comp1=out.error[indy,indx,j]*out.data[indy,indx,j]/amp.data[indy,indx]
        comp2=out.error[indy,indx,j+1]*out.data[indy,indx,j+1]/amp.data[indy,indx]
        amp.error[indy,indx] = np.sqrt(comp1**2 + comp2**2)
        # Phase Error (degrees)
        comp1=out.error[indy,indx,j]*out.data[indy,indx,j+1]/(amp.data[indy,indx]**2)
        comp2=out.error[indy,indx,j+1]*out.data[indy,indx,j]/(amp.data[indy,indx]**2)
        ph.error[indy,indx] = (180.0/np.pi)*np.sqrt(comp1**2 + comp2**2)

        # output file names for amplitude and phase
        f2 = (ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            flag, '_AMPL', START, END, suffix)
        f3 = (ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            flag, '_PHASE', START, END, suffix)
        file2 = OUTPUT_DIRECTORY.joinpath(output_format.format(*f2))
        file3 = OUTPUT_DIRECTORY.joinpath(output_format.format(*f3))
        # full attributes
        AMP_UNITS = 'centimeters'
        PH_UNITS = 'degrees'
        LONGNAME = 'Equivalent_Water_Thickness'
        AMP_TITLE = f'Sea_Level_Fingerprint_{amp_title[flag]}'
        PH_TITLE = f'Sea_Level_Fingerprint_{ph_title[flag]}'
        # Output seasonal amplitude and phase to files
        output_data(amp, FILENAME=file2, DATAFORM=DATAFORM,
            UNITS=AMP_UNITS, LONGNAME=LONGNAME, TITLE=AMP_TITLE,
            CONF=CONF, VERBOSE=VERBOSE, MODE=MODE)
        output_data(ph, FILENAME=file3, DATAFORM=DATAFORM,
            UNITS=PH_UNITS, LONGNAME='Phase', TITLE=PH_TITLE,
            CONF=CONF, VERBOSE=VERBOSE, MODE=MODE)
        # add output files to list object
        output_files.append(file2)
        output_files.append(file3)

    # Output fit significance
    signif_longname = {}
    signif_longname['SSE'] = 'Sum of Squares Error'
    signif_longname['AIC'] = 'Akaike information criterion'
    signif_longname['BIC'] = 'Bayesian information criterion'
    signif_longname['R2Adj'] = 'Adjusted Coefficient of Determination'
    # for each fit significance term
    for key,fs in FS.items():
        # output file names for fit significance
        signif_str = f'{key}_'
        f4 = (ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            signif_str, coef_str[ORDER], START, END, suffix)
        file4 = OUTPUT_DIRECTORY.joinpath(output_format.format(*f4))
        # full attributes
        LONGNAME = signif_longname[key]
        # output fit significance to file
        output_data(fs, FILENAME=file4, DATAFORM=DATAFORM,
            UNITS=key, LONGNAME=LONGNAME, TITLE=nu,
            VERBOSE=VERBOSE, MODE=MODE)
        # add output files to list object
        output_files.append(file4)

    # return the list of output files
    return output_files

# PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, DATAFORM=None, UNITS=None,
    LONGNAME=None, TITLE=None, CONF=0, VERBOSE=0, MODE=0o775):
    # field mapping for output regression data
    field_mapping = {}
    field_mapping['lat'] = 'lat'
    field_mapping['lon'] = 'lon'
    field_mapping['data'] = 'data'
    # add data error to output field mapping
    if hasattr(data, 'error'):
        field_mapping['error'] = 'error'
    # output attributes
    attributes = dict(ROOT={})
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    attributes['data'] = {}
    attributes['data']['description'] = 'Model_fit'
    attributes['data']['long_name'] = LONGNAME
    attributes['data']['units'] = UNITS
    attributes['error'] = {}
    attributes['error']['description'] = 'Uncertainty_in_model_fit'
    attributes['error']['long_name'] = LONGNAME
    attributes['error']['units'] = UNITS
    attributes['error']['confidence'] = 100*CONF
    # output global attributes
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # write to output file
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        data.to_ascii(FILENAME, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        data.to_netCDF4(FILENAME, date=False, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attributes,
            title=TITLE, reference=REFERENCE)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        data.to_HDF5(FILENAME, date=False, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attributes,
            title=TITLE, reference=REFERENCE)
    # change the permissions mode of the output file
    FILENAME.chmod(mode=MODE)

# PURPOSE: print a file log for the sea level regression calculation
def output_log_file(input_arguments, output_files):
    # format: sea_level_regress_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'sea_level_regress_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the sea level regression calculation
def output_error_log_file(input_arguments):
    # format: failed_sea_level_regress_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'failed_sea_level_regress_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads in sea level grid files and calculates the
            trends at each grid point following an input regression model
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int,
        help='Starting GRACE/GRACE_FO month for time series regression')
    parser.add_argument('--end','-E',
        type=int,
        help='Ending GRACE/GRACE_FO month for time series regression')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # sea level fingerprint parameters
    parser.add_argument('--iteration','-I',
        type=int, default=1,
        help='Sea level fingerprint iteration')
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # regression parameters
    # 0: mean
    # 1: trend
    # 2: acceleration
    parser.add_argument('--order',
        type=int, default=2,
        help='Regression fit polynomial order')
    # regression fit cyclical terms
    parser.add_argument('--cycles',
        type=float, default=[0.5,1.0,161.0/365.25], nargs='+',
        help='Regression fit cyclical terms')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # sea_level_regress_run_2002-04-01_PID-00000.log
    # failed_sea_level_regress_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run sea_level_regress algorithm with parameters
        output_files = sea_level_regress(
            args.center,
            args.release,
            args.product,
            args.lmax,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            DATAFORM=args.format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            ORDER=args.order,
            CYCLES=args.cycles,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
