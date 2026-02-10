#!/usr/bin/env python
u"""
make_sea_level_mascon_shells.py
Written by Tyler Sutterley (05/2023)

Computes and outputs spherical harmonics for an input set of mascon files
Creates a shell script for running sea level variation code

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for mascon files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    -b X, --body X: Treatment of the body tide Love number
        0: Wahr (1981) and Wahr (1985) values from PREM
        1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    -f X, --fluid X: Treatment of the fluid Love number
        0: Han and Wahr (1989) fluid love number
        1: Munk and MacDonald (1960) secular love number
        2: Munk and MacDonald (1960) fluid love number
        3: Lambeck (1980) fluid love number
    --polar-feedback: Include polar feedback
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
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --mascon-file X: index file of mascons spherical harmonics
    --coordinate-file X: File with spatial coordinates of mascon centers
    -H X, --header X:  Number of header lines to skip in coordinate file
    --mascon-type X: input load type (DISC, POINT or CAP)
    --redistribute-mascons: redistribute mascon mass over the ocean
    -I X, --iteration X: Sea level fingerprint iteration
    -e X, --expansion X: Spherical harmonic expansion for sea level fingerprints
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    legendre_polynomials.py: Computes fully normalized Legendre polynomials
    gen_disc_load.py: Computes geoid Stokes coefficients for a disc load
    gen_point_load.py: Computes geoid Stokes coefficients for point masses
    gen_spherical_cap.py: Computes geoid Stokes coefficients for a spherical cap
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    units.py: class for converting spherical harmonic data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: output an index listing the output harmonic files
        add option to use mascons from point masses
    Updated 01/2021: harmonics object output from disc and cap generators
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: updates to reading load love numbers
        using harmonics class for operating and outputting spherical harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: check which python version is running
    Updated 12/2018: added parallel processing with multiprocessing
    Updated 11/2018: can vary the land-sea mask
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2018: added option --mode to set the output file permissions
        simplified love number extrapolation if LMAX is greater than 696
    Updated 02/2018: added parameter EXPANSION to specify the SLF truncation
    Updated 08/2017: added flag to incorporate polar feedback
    Written 08/2017
"""
from __future__ import print_function

import sys
import os
import re
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

# PURPOSE: calculate harmonics for a set of mascons from a coordinate index
def make_sea_level_mascon_shells(PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    BODY_TIDE_LOVE=0,
    FLUID_LOVE=0,
    REFERENCE=None,
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    DATAFORM=None,
    MASCON_FILE=None,
    COORDINATE_FILE=None,
    HEADER=0,
    MASCON_TYPE=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    POLAR_FEEDBACK=False,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    VERBOSE=0,
    MODE=0o775):

    # recursively create output directory if not currently existing
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # number of GRACE/GRACE-FO months
    nmon = len(sorted(set(np.arange(START,END+1)) - set(MISSING)))

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # Gaussian smoothing string for radius RAD (if RAD == 0: none)
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # used filtered (destriped) coefficients for stripe effects
    ds_str = '_FL' if DESTRIPE else ''
    # input GIA stokes coefficients to get titles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    gia_str = f'_{GIA_Ylms_rate.title.upper()}' if GIA else ''
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''

    # sea level equation program
    child_program = 'run_sea_level_equation.py'
    # spherical harmonic expansion flag
    expansion_flag = f' --lmax {EXPANSION:d}'
    # pole tide feedback
    polar_flag = ' --polar-feedback' if POLAR_FEEDBACK else ''
    # Love number flags
    love_flag = f' --love {LOVE_NUMBERS:d}'
    body_flag = f' --body {BODY_TIDE_LOVE:d}'
    fluid_flag = f' --fluid {FLUID_LOVE:d}'
    # Set gravitational load love number of degree 1 to 0.027
    reference_flag = f' --reference {REFERENCE}'
    # input and output data have date information
    date_flag = ' --date'
    # land-sea mask to use (if not default 0.5x0.5)
    if LANDMASK:
        mask_flag = f' --mask {LANDMASK}'
    else:
        mask_flag = ''
    # number of iterations in sea level program
    # limit the program iterations at 6
    iter_flag = ' -I 6'
    # input and output format flag
    format_flag = f' --format {DATAFORM}'
    # verbosity flag
    verbosity_flag = ' --verbose' if VERBOSE else ''

    # read coordinate file
    COORDINATE_FILE = pathlib.Path(COORDINATE_FILE).expanduser().absolute()
    coord = np.loadtxt(COORDINATE_FILE, skiprows=HEADER)
    # column 1: cap number
    # column 2: longitude of center point
    # column 3: latitude of center point
    num = coord[:,0].astype(np.int64)
    lon = coord[:,1]
    lat = coord[:,2]
    n_crd = len(num)

    # read load love numbers
    LOVE = gravtk.load_love_numbers(EXPANSION, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # input mascon spherical harmonic datafiles
    # Read mascon index file and get contents
    MASCON_FILE = pathlib.Path(MASCON_FILE).expanduser().absolute()
    with MASCON_FILE.open(mode='r', encoding='utf8') as f:
        mascon_files = f.read().splitlines()
    # number of mascons
    n_mas = len(mascon_files)
    if (n_crd != n_mas):
        errmsg = f'Mismatching number of mascons ({n_crd:d},{n_mas:d})'
        logging.critical(errmsg)
    # spatial area of the mascon
    total_area = np.zeros((n_mas))
    # name of each mascon
    mascon_name = []
    # allocate for mascon Ylms expanded up to degree of sea level fingerprints
    mascon_Ylms = gravtk.harmonics(lmax=EXPANSION, mmax=EXPANSION)
    mascon_Ylms.clm = np.zeros((EXPANSION+1, EXPANSION+1, nmon))
    mascon_Ylms.slm = np.zeros((EXPANSION+1, EXPANSION+1, nmon))
    # for each mascon
    # calculate the spherical harmonic coefs
    # write spherical harmonics to file
    for i,fi in enumerate(mascon_files):
        # stem is the mascon file without directory or suffix
        # if lower case: will capitalize
        # if mascon name contains degree and order info: scrub from string
        mascon_file = pathlib.Path(fi).expanduser().absolute()
        stem = re.sub(r'_L(\d+)(M\d+)?', r'', mascon_file.stem.upper())
        mascon_name.append(stem)
        # input filename format (for both LMAX==MMAX and LMAX != MMAX cases):
        # mascon name, GRACE dataset, GIA model, LMAX, (MMAX,)
        # Gaussian smoothing, filter flag, remove reconstructed fields flag
        # output GRACE error file
        file_input='{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}.txt'.format(mascon_name[i],
            dset_str,gia_str,atm_str,ocean_str,LMAX,order_str,gw_str,ds_str)
        dinput = np.loadtxt(OUTPUT_DIRECTORY.joinpath(file_input))
        mascon_Ylms.month = dinput[:,0].astype(np.int64)
        mascon_Ylms.time = dinput[:,1].copy()
        total_area[i] = dinput[0,4].copy()
        # Calculate spherical harmonic coefficients for given type
        # truncate spherical harmonics at degree EXPANSION
        if (MASCON_TYPE == 'DISC'):
            # Calculate for a disc load using 1 Gt
            Ylms = gravtk.gen_disc_load(1.0,lon[i],lat[i],total_area[i],
                LMAX=EXPANSION,LOVE=LOVE)
        elif (MASCON_TYPE == 'POINT'):
            # Calculate for a point load using 1 Gt
            Ylms = gravtk.gen_point_load(np.array(1.0),lon[i],lat[i],
                LMAX=EXPANSION,UNITS=2,LOVE=LOVE)
        elif (MASCON_TYPE == 'CAP'):
            # Calculate for a spherical cap using 1 Gt
            Ylms = gravtk.gen_spherical_cap(1.0,lon[i],lat[i],LMAX=EXPANSION,
                AREA=total_area[i]*1e10,UNITS=2,LOVE=LOVE)
        # calculate total coefficients for each date
        for t in range(nmon):
            mascon_Ylms.clm[:,:,t] += dinput[t,2]*Ylms.clm[:,:]
            mascon_Ylms.slm[:,:,t] += dinput[t,2]*Ylms.slm[:,:]

    # list of output files
    output_files = []
    # create sea level shell script and index file
    args = (MASCON_TYPE, ITERATION, dset_str, gia_str, ocean_str, EXPANSION, START, END)
    f1 = '{0}_ITERATION_{1}_INDEX{2}{3}{4}_L{5:d}_{6:03d}-{7:03}.sh'.format(*args)
    f2 = '{0}_ITERATION_{1}_INDEX{2}{3}{4}_CLM_L{5:d}_{6:03d}-{7:03}.txt'.format(*args)
    output_shell_script = OUTPUT_DIRECTORY.joinpath(f1)
    output_index_file = OUTPUT_DIRECTORY.joinpath(f2)
    fid1 = output_shell_script.open(mode='w', encoding='utf8')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # print the path to the shell script
    logging.info(str(output_shell_script))
    # output file format for input_distribution and output_slf
    file_format = '{0}_ITERATION_{1}{2}{3}{4}_L{5:d}_{6:03d}.{7}'
    # formatting string for each line in the shell script
    shell_format='{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11} --mode {12} {13} {14}'
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # output file for each date
    for t,grace_month in enumerate(mascon_Ylms.month):
        # output to file formatted for use in the sea level equation functions
        args = (MASCON_TYPE, ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            grace_month, suffix[DATAFORM])
        input_load = OUTPUT_DIRECTORY.joinpath(file_format.format(*args))
        # output spherical harmonic file in data format
        Ylms = mascon_Ylms.index(t)
        Ylms.to_file(input_load, format=DATAFORM, **attributes)
        # change the permissions mode of the output file
        input_load.chmod(mode=MODE)
        # print file name to index
        output_files.append(input_load)
        # print shell script commands
        args = ('SLF', ITERATION, dset_str, gia_str, ocean_str, EXPANSION,
            grace_month, suffix[DATAFORM])
        output_slf = OUTPUT_DIRECTORY.joinpath(file_format.format(*args))
        # shell script command
        args = (child_program, expansion_flag, polar_flag, love_flag,
            body_flag, fluid_flag, reference_flag, date_flag, mask_flag,
            iter_flag, format_flag, verbosity_flag, oct(MODE),
            Ylms.compressuser(input_load),
            gravtk.spatial().compressuser(output_slf))
        print(shell_format.format(*args), file=fid1)
        # print harmonics to index file
        print(Ylms.compressuser(input_load), file=fid2)
    # close the shell script and index files
    fid1.close()
    fid2.close()
    # change the permissions mode of the shell script and index files
    output_shell_script.chmod(mode=MODE)
    output_index_file.chmod(mode=MODE)

    # return list of output files
    return output_files

# PURPOSE: print a file log for the mascon harmonic calculation
def output_log_file(input_arguments, output_files):
    # format: mascon_disc_run_2002-04-01_PID-70335.log
    TYPE = arguments.mascon_type.lower()
    args = (TYPE,time.strftime('%Y-%m-%d',time.localtime()),os.getpid())
    LOGFILE = 'mascon_{0}_run_{1}_PID-{2:d}.log'.format(*args)
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

# PURPOSE: print a error file log for the mascon harmonic calculation
def output_error_log_file(input_arguments):
    # format: failed_mascon_disc_run_2002-04-01_PID-70335.log
    TYPE = arguments.mascon_type.lower()
    args = (TYPE,time.strftime('%Y-%m-%d',time.localtime()),os.getpid())
    LOGFILE = 'failed_mascon_{0}_run_{1}_PID-{2:d}.log'.format(*args)
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
        description="""Computes spherical harmonics for a set of mascon files.
            Creates a shell script for running sea level variation code.
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
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # different treatments of the body tide Love numbers of degree 2
    # 0: Wahr (1981) and Wahr (1985) values from PREM
    # 1: Farrell (1972) values from Gutenberg-Bullen oceanic mantle model
    parser.add_argument('--body','-b',
        type=int, default=0, choices=[0,1],
        help='Treatment of the body tide Love number')
    # different treatments of the fluid Love number of gravitational potential
    # 0: Han and Wahr (1989) fluid love number
    # 1: Munk and MacDonald (1960) secular love number
    # 2: Munk and MacDonald (1960) fluid love number
    # 3: Lambeck (1980) fluid love number
    parser.add_argument('--fluid','-f',
        type=int, default=0, choices=[0,1,2,3],
        help='Treatment of the fluid Love number')
    # option for polar feedback
    parser.add_argument('--polar-feedback',
        default=False, action='store_true',
        help='Include effects of polar feedback')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
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
    # use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=pathlib.Path,
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--coordinate-file',
        type=pathlib.Path,
        required=True,
        help='File with spatial coordinates of mascon centers')
    # number of header lines to skip in coordinate file
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines to skip in coordinate file')
    # input load type (DISC, POINT or CAP)
    parser.add_argument('--mascon-type','-T',
        type=str.upper, default='CAP', choices=['DISC','POINT','CAP'],
        help='Input load type')
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
    # land-sea mask for redistributing mascon mass and land water flux
    parser.add_argument('--mask',
        type=pathlib.Path,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # mascon_disc_run_2002-04-01_PID-00000.log
    # failed_mascon_disc_run_2002-04-01_PID-00000.log
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
        # run make_sea_level_mascon_shells algorithm with parameters
        output_files = make_sea_level_mascon_shells(
            args.center,
            args.release,
            args.product,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            BODY_TIDE_LOVE=args.body,
            FLUID_LOVE=args.fluid,
            REFERENCE=args.reference,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            DATAFORM=args.format,
            MASCON_FILE=args.mascon_file,
            COORDINATE_FILE=args.coordinate_file,
            HEADER=args.header,
            MASCON_TYPE=args.mascon_type,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            POLAR_FEEDBACK=args.polar_feedback,
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
