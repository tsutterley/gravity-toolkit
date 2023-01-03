#!/usr/bin/env python
u"""
calc_sensitivity_kernel.py
Written by Tyler Sutterley (01/2023)

Calculates spatial sensitivity kernels through a least-squares mascon procedure

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for mascon files
    --lmin X: minimum spherical harmonic degree
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input and output data format
        ascii
        netCDF4
        HDF5
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --mascon-file X: index file of mascons spherical harmonics
    --redistribute-mascons: redistribute mascon mass over the ocean
    --fit-method X: method for fitting sensitivity kernel to harmonics
        1: mass coefficients
        2: geoid coefficients
    -s, --spatial: Output spatial grid file for each mascon
    -S X, --spacing X: spatial resolution of output data (dlon,dlat)
    -I X, --interval X: Output grid interval
        1: global
        2: centered global
        3: non-global
    -B X, --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    I Velicogna, T C Sutterley and M R van den Broeke. "Regional acceleration
        in ice mass loss from Greenland and Antarctica using GRACE
        time-variable gravity data". Geophysical Research Letters,
        41(22):8130-8137, 2014. https://doi.org/10.1002/2014GL061052

    T Jacob, J Wahr, W Pfeffer, and S C Swenson "Recent contributions of
        glaciers and ice caps to sea level rise". Nature, 482, 514-518 (2012).
        https://doi.org/10.1038/nature10847

    V M Tiwari, J Wahr, S and Swenson, "Dwindling groundwater resources in
        northern India, from satellite gravity observations",
        Geophysical Research Letters, 36(18), L18401, (2009).
        https://doi.org/10.1029/2009GL039401

UPDATE HISTORY:
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 07/2022: create mask for output gridded variables
        made creating the spatial outputs optional to improve compute time
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: simplified file imports using wrappers in harmonics
        added path to default land-sea mask for mass redistribution
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 10/2019: changing Y/N flags to True/False
    Updated 10/2018: verify integers for python3 compatibility
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2018: added extrapolation of load love numbers if LMAX > 696
    Updated 09/2017: use a different land-sea mask for calculating ocean_Ylms
        use rcond=-1 in numpy least-squares algorithm
    Updated 05/2016: using __future__ print function
    Updated 02/2016: direct calculation of number of harmonics n_harm
        use getopt parameters to set number of PROCESSES to run in parallel,
            whether or not to output a log file, added new help module
    Updated 11/2015: create unique log filenames
    Updated 07/2015: added output of the sensitivity kernel Ylms in addition
        to the spatial fields (rather than just the spatial fields)
        will output logs with parameters and output_files
        added multiprocessing error handling with traceback
    Updated 05/2015: added parameter MMAX for LMAX != MMAX
        added portion to redistribute mascon mass uniformly over the ocean
    Updated 10/2014: distributed computing with the multiprocessing module
        added INTERVAL parameter for (degree spacing)/2
        input/output file type (ascii, netCDF4, HDF5)
    Updated 05/2014: added import functions
    Updated 02/2014: updated comments and added os.path.joins for connecting
        directories and files (generalizing code)
        some general updates to the program code
    Updated 08/2013: general updates to inputting data
    Updated 03/2012: edited to use new gen_stokes time-series option
    Written 02/2012
"""
from __future__ import print_function, division

import sys
import os
import re
import time
import logging
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate a regional time-series through a least
# squares mascon process
def calc_sensitivity_kernel(LMAX, RAD,
    LMIN=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MASCON_FILE=None,
    REDISTRIBUTE_MASCONS=False,
    FIT_METHOD=0,
    LANDMASK=None,
    SPATIAL=False,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]
    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)

    # Create output Directory if not currently existing
    if (not os.access(OUTPUT_DIRECTORY,os.F_OK)):
        os.mkdir(OUTPUT_DIRECTORY)

    # list object of output files for file logs (full path)
    output_files = []

    # read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
    # Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    # Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    # input/output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_MASCONS:
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX, MMAX=MMAX,
            LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        # not distributing uniformly over ocean
        ocean_str = ''

    # input mascon spherical harmonic datafiles
    with open(MASCON_FILE, mode='r', encoding='utf8') as f:
        mascon_files = [l for l in f.read().splitlines() if parser.match(l)]
    # number of mascons
    n_mas = len(mascon_files)
    # spatial area of the mascon
    total_area = np.zeros((n_mas))
    # name of each mascon
    mascon_name = []
    # for each valid file in the index (iterate over mascons)
    mascon_list = []
    for k,fi in enumerate(mascon_files):
        # read mascon spherical harmonics
        Ylms = gravtk.harmonics().from_file(os.path.expanduser(fi),
            format=DATAFORM, date=False)
        # Calculating the total mass of each mascon (1 cmwe uniform)
        total_area[k] = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
        # distribute mascon mass uniformly over the ocean
        if REDISTRIBUTE_MASCONS:
            # calculate ratio between total mascon mass and
            # a uniformly distributed cm of water over the ocean
            ratio = Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
            # for each spherical harmonic
            for m in range(0,MMAX+1):# MMAX+1 to include MMAX
                for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                    # remove ratio*ocean Ylms from mascon Ylms
                    # note: x -= y is equivalent to x = x - y
                    Ylms.clm[l,m] -= ratio*ocean_Ylms.clm[l,m]
                    Ylms.slm[l,m] -= ratio*ocean_Ylms.slm[l,m]
        # truncate mascon spherical harmonics to d/o LMAX/MMAX and add to list
        mascon_list.append(Ylms.truncate(lmax=LMAX, mmax=MMAX))
        # mascon base is the file without directory or suffix
        mascon_base = os.path.basename(mascon_files[k])
        mascon_base = os.path.splitext(mascon_base)[0]
        # if lower case, will capitalize
        mascon_base = mascon_base.upper()
        # if mascon name contains degree and order info, remove
        mascon_name.append(mascon_base.replace(f'_L{LMAX:d}', ''))
    # create single harmonics object from list
    mascon_Ylms = gravtk.harmonics().from_list(mascon_list, date=False)
    # clear mascon list variable
    del mascon_list

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # Initialing harmonics for least squares fitting
    # mascon kernel
    M_lm = np.zeros((n_harm,n_mas))
    # mascon kernel converted to output unit
    MA_lm = np.zeros((n_harm,n_mas))
    # sensitivity kernel
    A_lm = np.zeros((n_harm,n_mas))
    # Initializing conversion factors
    # factor for converting to smoothed coefficients of mass
    fact = np.zeros((n_harm))
    # factor for converting back into geoid coefficients
    fact_inv = np.zeros((n_harm))
    # smoothing factor
    wt_lm = np.zeros((n_harm))

    # ii is a counter variable for building the mascon column array
    ii = 0
    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    # Calculating factor to convert geoid spherical harmonic coefficients
    # to coefficients of mass (Wahr, 1998)
    coeff = rho_e*rad_e/3.0
    coeff_inv = 0.75/(np.pi*rho_e*rad_e**3)
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # copy cosine and sin harmonics
        mascon_harm = getattr(mascon_Ylms, csharm)
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # Mascon Spherical Harmonics
                M_lm[ii,:] = np.copy(mascon_harm[l,m,:])
                # degree dependent factor to convert to mass
                fact[ii] = (2.0*l + 1.0)/(1.0 + kl[l])
                # degree dependent factor to convert from mass
                fact_inv[ii] = coeff_inv*(1.0 + kl[l])/(2.0*l+1.0)
                # degree dependent smoothing
                wt_lm[ii] = np.copy(wt[l])
                # add 1 to counter
                ii += 1

    # Converting mascon coefficients to fit method
    if (FIT_METHOD == 1):
        # Fitting Sensitivity Kernel as mass coefficients
        # converting M_lm to mass coefficients of the kernel
        for i in range(n_harm):
            MA_lm[i,:] = M_lm[i,:]*wt_lm[i]*fact[i]
        fit_factor = wt_lm*fact
        inv_fit_factor = np.copy(fact_inv)
    else:
        # Fitting Sensitivity Kernel as geoid coefficients
        for i in range(n_harm):
            MA_lm[:,:] = M_lm[i,:]*wt_lm[i]
        fit_factor = wt_lm*np.ones((n_harm))
        inv_fit_factor = np.ones((n_harm))

    # Fitting the sensitivity kernel from the input kernel
    for i in range(n_harm):
        # setting kern_i equal to 1 for d/o
        kern_i = np.zeros((n_harm))
        # converting to mass coefficients if specified
        kern_i[i] = 1.0*fit_factor[i]
        # spherical harmonics solution for the
        # mascon sensitivity kernels
        # Least Squares Solutions: Inv(X'.X).(X'.Y)
        kern_lm = np.linalg.lstsq(MA_lm, kern_i, rcond=-1)[0]
        for k in range(n_mas):
            A_lm[i,k] = kern_lm[k]*total_area[k]
    # free up larger variables
    del M_lm, MA_lm, wt_lm, fact, fact_inv, fit_factor

    # reshaping harmonics of sensitivity kernel to LMAX+1,MMAX+1
    # calculating the spatial sensitivity kernel of each mascon
    # kernel calculated as outlined in Tiwari (2009) and Jacobs (2012)
    # Initializing output sensitivity kernel (both spatial and Ylms)
    kern_Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    kern_Ylms.clm = np.zeros((LMAX+1,MMAX+1,n_mas))
    kern_Ylms.slm = np.zeros((LMAX+1,MMAX+1,n_mas))
    kern_Ylms.time = np.copy(total_area)
    # counter variable for deconstructing the mascon column arrays
    ii = 0
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # inv_fit_factor: normalize from mass harmonics
                temp = getattr(kern_Ylms, csharm)
                temp[l,m,:] = inv_fit_factor[ii]*A_lm[ii,:]
                # add 1 to counter
                ii += 1
    # free up larger variables
    del A_lm, inv_fit_factor

    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # for each mascon
    for k in range(n_mas):
        # get harmonics for mascon
        Ylms = kern_Ylms.index(k, date=False)
        # output sensitivity kernel to file
        args = (mascon_name[k],ocean_str,LMAX,order_str,gw_str,suffix)
        FILE1 = '{0}_SKERNEL_CLM{1}_L{2:d}{3}{4}.{5}'.format(*args)
        Ylms.to_file(os.path.join(OUTPUT_DIRECTORY, FILE1),
            format=DATAFORM, date=False, **attributes)
        # change the permissions mode
        os.chmod(os.path.join(OUTPUT_DIRECTORY,FILE1),MODE)
        # add output files to list object
        output_files.append(os.path.join(OUTPUT_DIRECTORY,FILE1))

    # attributes for output files
    attributes = {}
    attributes['units'] = 'unitless'
    attributes['longname'] = 'Sensitivity_Kernel'
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # if outputting spatial grids
    if SPATIAL:
        # Output spatial data object
        grid = gravtk.spatial()
        # Output Degree Spacing
        dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
        # Output Degree Interval
        if (INTERVAL == 1):
            # (-180:180,90:-90)
            n_lon = np.int64((360.0/dlon)+1.0)
            n_lat = np.int64((180.0/dlat)+1.0)
            grid.lon = -180 + dlon*np.arange(0,n_lon)
            grid.lat = 90.0 - dlat*np.arange(0,n_lat)
        elif (INTERVAL == 2):
            # (Degree spacing)/2
            grid.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
            grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
            n_lon = len(grid.lon)
            n_lat = len(grid.lat)
        elif (INTERVAL == 3):
            # non-global grid set with BOUNDS parameter
            minlon,maxlon,minlat,maxlat = BOUNDS.copy()
            grid.lon = np.arange(minlon+dlon/2.0,maxlon+dlon/2.0,dlon)
            grid.lat = np.arange(maxlat-dlat/2.0,minlat-dlat/2.0,-dlat)
            nlon = len(grid.lon)
            nlat = len(grid.lat)

        # Computing plms for converting to spatial domain
        theta = (90.0-grid.lat)*np.pi/180.0
        PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

        # for each mascon
        for k in range(n_mas):
            # get harmonics for mascon
            Ylms = kern_Ylms.index(k, date=False)
            # convert spherical harmonics to output spatial grid
            grid.data = gravtk.harmonic_summation(Ylms.clm, Ylms.slm,
                grid.lon, grid.lat, LMAX=LMAX, MMAX=MMAX, PLM=PLM).T
            grid.mask = np.zeros_like(grid.data, dtype=bool)
            # output sensitivity kernel to file
            args = (mascon_name[k],ocean_str,LMAX,order_str,gw_str,suffix)
            FILE2 = '{0}_SKERNEL{1}_L{2:d}{3}{4}.{5}'.format(*args)
            grid.to_file(os.path.join(OUTPUT_DIRECTORY,FILE2),
                format=DATAFORM, date=False, **attributes)
            # change the permissions mode
            os.chmod(os.path.join(OUTPUT_DIRECTORY,FILE2),MODE)
            # add output files to list object
            output_files.append(os.path.join(OUTPUT_DIRECTORY,FILE2))

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the mascon sensitivity kernel analysis
def output_log_file(input_arguments, output_files):
    # format: calc_skernel_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_skernel_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
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

# PURPOSE: print a error file log for the mascon sensitivity kernel analysis
def output_error_log_file(input_arguments):
    # format: calc_skernel_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_skernel_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
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
        description="""Calculates spatial sensitivity kernels through a
            least-squares mascon procedure
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for mascon files')
    # minimum spherical harmonic degree
    parser.add_argument('--lmin',
        type=int, default=1,
        help='Minimum spherical harmonic degree')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for auxiliary files')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # 1: mass coefficients
    # 2: geoid coefficients
    parser.add_argument('--fit-method',
        type=int, default=1, choices=(1,2),
        help='Method for fitting sensitivity kernel to harmonics')
    # land-sea mask for redistributing mascon mass
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing mascon mass')
    # output spatial grid
    parser.add_argument('--spatial','-s',
        default=False, action='store_true',
        help='Output spatial grid file for each mascon')
    # output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds','-B',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    # Output log file for each job in forms
    # calc_skernel_run_2002-04-01_PID-00000.log
    # calc_skernel_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run calc_sensitivity_kernel algorithm with parameters
        output_files = calc_sensitivity_kernel(
            args.lmax,
            args.radius,
            LMIN=args.lmin,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            MASCON_FILE=args.mascon_file,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            FIT_METHOD=args.fit_method,
            LANDMASK=args.mask,
            SPATIAL=args.spatial,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as e:
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
