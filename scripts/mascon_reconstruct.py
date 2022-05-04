#!/usr/bin/env python
u"""
mascon_reconstruct.py
Written by Tyler Sutterley (04/2022)

Calculates the equivalent spherical harmonics from a mascon time series

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for mascon files
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
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
    --atm-correction: Apply atmospheric jump correction coefficients
    --mask X: Land-sea mask for redistributing mascon mass
    --mascon-file X: index file of mascons spherical harmonics
    --redistribute-mascons: redistribute mascon mass over the ocean
    --reconstruct-file X: reconstructed mascon time series file
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        added path to default land-sea mask for mass redistribution
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 06/2020: using harmonics class for spherical harmonic operations
    Updated 04/2020: updates to reading load love numbers
        reading land-sea mask as a parameter
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: use a different land-sea mask for calculating ocean_Ylms
    Updated 03/2018: include a flag if using atmospheric ECMWF "jump" files
    Updated 02/2017: added parameter to redistribute mascon mass over the ocean
    Updated 06/2016: using __future__ print function
    Updated 02/2016: use getopt parameters to set number of PROCESSES to run
    Updated 06/2015: output harmonics for each mascon
    Updated 05/2015: added parameter MMAX for LMAX != MMAX.
        minor update to have parameters converted in function.
        added ascii and HDF5 output option
    Updated 01/2015: added error handling for multiprocessing threads
    Updated 10/2014: Distribute computing with multiprocessing module
    Updated 09/2014: Converted to function with main args
    Updated 05/2014
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import numpy as np
import traceback

import gravity_toolkit.utilities as utilities
from gravity_toolkit.read_GIA_model import read_GIA_model
from gravity_toolkit.read_love_numbers import load_love_numbers
from gravity_toolkit.ocean_stokes import ocean_stokes
from gravity_toolkit.harmonics import harmonics
from gravity_toolkit.units import units

#-- PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: tilde-compress a file path string
def tilde_compress(file_path):
    return file_path.replace(os.path.expanduser('~'),'~')

#-- PURPOSE: Reconstruct spherical harmonic fields from the mascon
#-- time series calculated in calc_mascon
def mascon_reconstruct(DSET, LMAX, RAD,
    START=None,
    END=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    DATAFORM=None,
    MASCON_FILE=None,
    REDISTRIBUTE_MASCONS=False,
    RECONSTRUCT_FILE=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    #-- for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else '_{0}'.format(DSET)
    #-- atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    #-- Gaussian smoothing string for radius RAD
    gw_str = '_r{0:0.0f}km'.format(RAD) if (RAD != 0) else ''
    #-- input GIA spherical harmonic datafiles
    GIA_Ylms_rate = read_GIA_model(GIA_FILE,GIA=GIA,LMAX=LMAX,MMAX=MMAX)
    gia_str = '_{0}'.format(GIA_Ylms_rate['title']) if GIA else ''
    #-- output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = 'M{0:d}'.format(MMAX) if (MMAX != LMAX) else ''
    #-- filter grace coefficients flag
    ds_str = '_FL' if DESTRIPE else ''
    #-- output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    #-- file parser for reading index files
    #-- removes commented lines (can comment out files in the index)
    #-- removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)

    #-- create initial reconstruct index for calc_mascon.py
    fid = open(RECONSTRUCT_FILE,'w')
    #-- output file format
    file_format = '{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}_{9:03d}-{10:03d}.{11}'

    #-- read load love numbers
    hl,kl,ll = load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)
    #-- Earth Parameters
    factors = units(lmax=LMAX).harmonic(hl,kl,ll)
    #-- Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    #-- Average Radius of the Earth [cm]
    rad_e = factors.rad_e
    #-- Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_MASCONS:
        #-- read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = ocean_stokes(LANDMASK,LMAX,MMAX=MMAX,LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        #-- not distributing uniformly over ocean
        ocean_str = ''

    #-- input mascon spherical harmonic datafiles
    with open(MASCON_FILE, mode='r', encoding='utf8') as f:
        mascon_files = [l for l in f.read().splitlines() if parser.match(l)]
    for k,fi in enumerate(mascon_files):
        #-- read mascon spherical harmonics
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            Ylms=harmonics().from_ascii(os.path.expanduser(fi),date=False)
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            Ylms=harmonics().from_netCDF4(os.path.expanduser(fi),date=False)
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            Ylms=harmonics().from_HDF5(os.path.expanduser(fi),date=False)
        #-- Calculating the total mass of each mascon (1 cmwe uniform)
        total_area = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
        #-- distribute mascon mass uniformly over the ocean
        if REDISTRIBUTE_MASCONS:
            #-- calculate ratio between total mascon mass and
            #-- a uniformly distributed cm of water over the ocean
            ratio = Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
            #-- for each spherical harmonic
            for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
                for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
                    #-- remove ratio*ocean Ylms from mascon Ylms
                    #-- note: x -= y is equivalent to x = x - y
                    Ylms.clm[l,m] -= ratio*ocean_Ylms.clm[l,m]
                    Ylms.slm[l,m] -= ratio*ocean_Ylms.slm[l,m]
        #-- truncate mascon spherical harmonics to d/o LMAX/MMAX
        Ylms = Ylms.truncate(lmax=LMAX, mmax=MMAX)
        #-- mascon base is the file without directory or suffix
        mascon_base = os.path.basename(fi)
        mascon_base = os.path.splitext(mascon_base)[0]
        #-- if lower case, will capitalize
        mascon_base = mascon_base.upper()
        #-- if mascon name contains degree and order info, remove
        mascon_name = mascon_base.replace('_L{0:d}'.format(LMAX),'')

        #-- input filename format (for both LMAX==MMAX and LMAX != MMAX cases):
        #-- mascon name, GRACE dataset, GIA model, LMAX, (MMAX,)
        #-- Gaussian smoothing, filter flag, remove reconstructed fields flag
        #-- output GRACE error file
        args = (mascon_name,dset_str,gia_str.upper(),atm_str,ocean_str,
            LMAX,order_str,gw_str,ds_str)
        file_input = '{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}.txt'.format(*args)
        mascon_data_input=np.loadtxt(os.path.join(OUTPUT_DIRECTORY,file_input))

        #-- convert mascon time-series from Gt to cmwe
        mascon_sigma = 1e15*mascon_data_input[:,2]/total_area
        #-- mascon time-series Ylms
        mascon_Ylms = Ylms.scale(mascon_sigma)
        mascon_Ylms.time = mascon_data_input[:,1].copy()
        mascon_Ylms.month = mascon_data_input[:,0].astype(np.int64)

        #-- output to file: no ascii option
        args = (mascon_name,dset_str,gia_str.upper(),atm_str,ocean_str,
            LMAX,order_str,gw_str,ds_str,START,END,suffix[DATAFORM])
        FILE = file_format.format(*args)
        #-- output harmonics to file
        if (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            mascon_Ylms.to_netCDF4(os.path.join(OUTPUT_DIRECTORY,FILE))
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            mascon_Ylms.to_HDF5(os.path.join(OUTPUT_DIRECTORY,FILE))
        #-- print file name to index
        print(tilde_compress(os.path.join(OUTPUT_DIRECTORY,FILE)),file=fid)
        #-- change the permissions mode
        os.chmod(os.path.join(OUTPUT_DIRECTORY,FILE),MODE)
    #-- close the reconstruct index
    fid.close()
    #-- change the permissions mode of the index file
    os.chmod(RECONSTRUCT_FILE,MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
            description="""Calculates the equivalent spherical
            harmonics from a mascon time series
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for mascon files')
    #-- GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    #-- start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    #-- different treatments of the load Love numbers
    #-- 0: Han and Wahr (1995) values from PREM
    #-- 1: Gegout (2005) values from PREM
    #-- 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    #-- option for setting reference frame for gravitational load love number
    #-- reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    #-- Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    #-- Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    #-- GIA model type list
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
    #-- GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    #-- full path to GIA file
    parser.add_argument('--gia-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GIA file to read')
    #-- use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    #-- input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for auxiliary files')
    #-- mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    #-- mascon reconstruct parameters
    parser.add_argument('--reconstruct-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Reconstructed mascon time series file')
    #-- land-sea mask for redistributing mascon mass
    lsmask = utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing mascon mass')
    #-- print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    #-- try to run the analysis with listed parameters
    try:
        info(args)
        #-- run mascon_reconstruct algorithm with parameters
        mascon_reconstruct(
            args.product,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            DATAFORM=args.format,
            MASCON_FILE=args.mascon_file,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            RECONSTRUCT_FILE=args.reconstruct_file,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as e:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        logging.critical('process id {0:d} failed'.format(os.getpid()))
        logging.error(traceback.format_exc())

#-- run main program
if __name__ == '__main__':
    main()
