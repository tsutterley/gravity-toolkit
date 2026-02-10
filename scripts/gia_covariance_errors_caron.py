#!/usr/bin/env python
u"""
gia_covariance_errors_caron.py
Written by Yara Mohajerani (05/2019)
Updated by Tyler Sutterley (05/2023)

Calculate GIA errors based on full covariance matrix as given by
Caron et al [2019] https://doi.org/10.1002/2017GL076644

Error is calculated based on equation (3) in Wahr et al [2006]
https://doi.org/10.1029/2005GL025305

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for mascon files
    --lmin X: minimum spherical harmonic degree
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input data format for auxiliary files
        ascii
        netCDF4
        HDF5
    --mask X: Land-sea mask for redistributing mascon mass
    --mascon-file: index file of mascons spherical harmonics
    --redistribute-mascons: redistribute mascon mass over the ocean
    --fit-method: method for fitting sensitivity kernel to harmonics
        1: mass coefficients
        2: geoid coefficients
    -s X, --solver X: Least squares solver for sensitivity kernels
        inv: matrix inversion
        lstsq: least squares solution
        gelsy: complete orthogonal factorization
        gelss: singular value decomposition (SVD)
        gelsd: singular value decomposition (SVD) with divide and conquer method
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://scipy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 04/2023: add options for least-squares solver
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: added path to default land-sea mask for mass redistribution
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: add parser object for removing commented or empty lines
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
        use utilities to define path to load love numbers file
    Updated 04/2020: using the harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 10/2019: complete rewrite of program to be similar to calc_mascon.py
    Updated 06/2019: Instead of calculating errors individually for each cap
        add up kernel for given cap numbers and output 1 error
    Written 05/2019
"""
from __future__ import print_function, division

import sys
import os
import re
import logging
import pathlib
import argparse
import traceback
import numpy as np
import scipy.linalg
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: read GIA covariance matrix file and flatten to mascon form
def read_covariance_file(infile, LMAX, LMIN=None, MMAX=None):
    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)
    # read covariance matrix file
    infile = pathlib.Path(infile).expanduser().absolute()
    with infile.open(mode='r', encoding='utf8') as fid:
        # read the input file, split at lines and remove all commented lines
        contents = [i for i in fid.read().splitlines() if parser.match(i)]
    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_clm = (LMAX**2 - LMIN**2 + 3*LMAX - (LMAX-MMAX)**2 - (LMAX-MMAX))//2 + 1
    n_harm = int(LMAX**2 - LMIN**2 + 2*LMAX - (LMAX-MMAX)**2 - (LMAX-MMAX)) + 1
    # flattened covariance matrix
    cov = np.zeros((n_harm,n_harm))
    # for each line in the file
    for line in contents:
        # note clm orders are +ve and slm orders are -ve
        l1,m1,l2,m2,Ylms = line.split()
        l1,m1,l2,m2 = np.array([l1,m1,l2,m2],dtype=np.int64)
        # indice for filling flattened covariance matrix for lm harmonics
        if (m1 >= 0) and (l1 >= LMIN) and (l1 <= LMAX) and (m1 <= MMAX):
            # cosine harmonics
            i1 = m1 + ((l1-1)**2 - LMIN**2 + 3*(l1-1) - \
                np.max([l1-MMAX-1,0])**2 - np.max([l1-MMAX-1,0]))//2 + 1
        elif (l1 >= LMIN) and (l1 <= LMAX) and (np.abs(m1) <= MMAX):
            # sine harmonics
            i1 = n_clm + np.abs(m1) + ((l1-1)**2 - LMIN**2 + (l1-1) - \
                np.max([l1-MMAX-1,0])**2 - np.max([l1-MMAX-1,0]))//2
        else:
            i1 = None
        # indice for filling flattened covariance matrix for pq harmonics
        if (m2 >= 0) and (l2 >= LMIN) and (l2 <= LMAX) and (m2 <= MMAX):
            # cosine harmonics
            i2 = m2 + ((l2-1)**2 - LMIN**2 + 3*(l2-1) - \
                np.max([l2-MMAX-1,0])**2 - np.max([l2-MMAX-1,0]))//2 + 1
        elif (l2 >= LMIN) and (l2 <= LMAX) and (np.abs(m2) <= MMAX):
            # sine harmonics
            i2 = n_clm + np.abs(m2) + ((l2-1)**2 - LMIN**2 + (l2-1) - \
                np.max([l2-MMAX-1,0])**2 - np.max([l2-MMAX-1,0]))//2
        else:
            i2 = None
        # add data to flattened covariance matrix
        if (i1 and i2):
            cov[i1,i2] = np.float64(Ylms)
    # free up memory
    contents = None
    return cov

# calculate GIA error for given mascon configuration
def gia_covariance_errors_caron(base_dir, LMAX, RAD,
    LMIN=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    MASCON_FILE=None,
    MASCON_FORMAT=None,
    REDISTRIBUTE_MASCONS=False,
    FIT_METHOD=0,
    SOLVER=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # input directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)
    # list object of output files for file logs (full path)
    output_files = []

    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    # Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    # Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''

    # Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_MASCONS:
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX, MMAX=MMAX,
            LOVE=LOVE)
        ocean_str = '_OCN'
    else:
        # not distributing uniformly over ocean
        ocean_str = ''

    # input mascon spherical harmonic datafiles
    MASCON_FILE = pathlib.Path(MASCON_FILE).expanduser().absolute()
    with MASCON_FILE.open(mode='r', encoding='utf8') as f:
        mascon_files = [l for l in f.read().splitlines() if parser.match(l)]
    # number of mascons
    n_mas = len(mascon_files)
    # spatial area of the mascon
    area_tot = np.zeros((n_mas))
    # name of each mascon
    mascon_name = []
    # for each valid file in the index (iterate over mascons)
    mascon_list = []
    for k in range(n_mas):
        # read mascon spherical harmonics
        Ylms = gravtk.harmonics().from_file(mascon_files[k],
            format=MASCON_FORMAT, date=False)
        # Calculating the total mass of each mascon (1 cmwe uniform)
        area_tot[k] = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
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
        # stem is the mascon file without directory or suffix
        # if lower case: will capitalize
        # if mascon name contains degree and order info: scrub from string
        stem = re.sub(r'_L(\d+)(M\d+)?', r'', Ylms.filename.stem.upper())
        mascon_name.append(stem)
    # create single harmonics object from list
    mascon_Ylms = gravtk.harmonics().from_list(mascon_list, date=False)
    # clear mascon list variable
    del mascon_list

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # read Caron et al. (2018) GIA covariance matrix files
    # list containing files to read based on spherical harmonic degree range
    gia_files = []
    gia_files.append('covStokes_GIA_deg_2_to_60.txt')
    if (LMAX > 60):
        gia_files.append('covStokes_GIA_deg_61_to_75.txt')
    if (LMAX > 75):
        gia_files.append('covStokes_GIA_deg_76_to_83.txt')
    if (LMAX > 83):
        gia_files.append('covStokes_GIA_deg_84_to_89.txt')
    # flattened combined covariance matrix
    gia_cov = np.zeros((n_harm, n_harm))
    # for each GIA covariance file
    for fi in gia_files:
        gia_cov[:,:] += read_covariance_file(base_dir.joinpath(fi), LMAX,
            LMIN=LMIN, MMAX=MMAX)
    # GIA title string for covariance-derived errors
    gia_str = '_Caron_Error'

    # Initialing harmonics for least squares fitting
    # mascon kernel
    M_lm = np.zeros((n_harm, n_mas))
    # mascon kernel converted to output unit
    MA_lm = np.zeros((n_harm, n_mas))
    # sensitivity kernel
    A_lm = np.zeros((n_harm, n_mas))
    # Initializing conversion factors
    # factor for converting to smoothed coefficients of mass
    fact = np.zeros((n_harm))
    # smoothing factor
    wt_lm = np.zeros((n_harm))
    # total mascon error from covariance matrix
    M_err = np.zeros((n_mas))

    # ii is a counter variable for building the mascon column array
    ii = 0
    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    # Calculating factor to convert geoid spherical harmonic coefficients
    # to coefficients of mass (Wahr, 1998)
    coeff = rho_e*rad_e/3.0
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # copy cosine and sin harmonics
        mascon_harm = getattr(mascon_Ylms, csharm)
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mn = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mn+1):
                # Mascon Spherical Harmonics
                M_lm[ii,:] = np.copy(mascon_harm[l,m,:])
                # degree dependent factor to convert to mass
                fact[ii] = (2.0*l + 1.0)/(1.0 + LOVE.kl[l])
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
    elif (FIT_METHOD == 2):
        # Fitting Sensitivity Kernel as geoid coefficients
        for i in range(n_harm):
            MA_lm[:,:] = M_lm[i,:]*wt_lm[i]
        fit_factor = wt_lm*np.ones((n_harm))

    # Fitting the sensitivity kernel from the input kernel
    for i in range(n_harm):
        # setting kern_i equal to 1 for d/o
        kern_i = np.zeros((n_harm))
        # converting to mass coefficients if specified
        kern_i[i] = 1.0*fit_factor[i]
        # spherical harmonics solution for the
        # mascon sensitivity kernels
        if (SOLVER == 'inv'):
            kern_lm = np.dot(np.linalg.inv(MA_lm), kern_i)
        elif (SOLVER == 'lstsq'):
            kern_lm = np.linalg.lstsq(MA_lm, kern_i, rcond=-1)[0]
        elif SOLVER in ('gelsd', 'gelsy', 'gelss'):
            kern_lm, res, rnk, s = scipy.linalg.lstsq(MA_lm, kern_i,
                lapack_driver=SOLVER)
        # calculate the sensitivity kernel for each mascon
        for k in range(n_mas):
            A_lm[i,k] = kern_lm[k]*area_tot[k]

    # calculate total error for each kernel
    # for each spherical harmonic lm (order is [C00...Clm,S11...Slm])
    for i in range(n_harm):
        # for each spherical harmonic pq (order is [C00...Cpq,S11...Spq])
        for j in range(n_harm):
            # add to total mascon errors
            M_err[:] += (A_lm[i,:]*A_lm[j,:]*gia_cov[i,j])

    # for each mascon
    for k in range(n_mas):
        # output filename format (for both LMAX==MMAX and LMAX != MMAX cases):
        # mascon name, GIA model, LMAX, (MMAX,), Gaussian smoothing, filter
        fargs = (mascon_name[k],gia_str.upper(),ocean_str,LMAX,order_str,gw_str)
        file_format = '{0}{1}{2}_L{3:d}{4}{5}.txt'
        output_file = OUTPUT_DIRECTORY.joinpath(file_format.format(*fargs))
        # take sqrt and convert from from g to gigatonnes
        args = (np.sqrt(M_err[k])/1e15, area_tot[k]/1e10)
        with output_file.open(mode='w', encoding='utf8') as fid1:
            print('{0:16.10f} {1:16.10f}'.format(*args), file=fid1)
        # change the permissions mode
        output_file.chmod(mode=MODE)
        # add output files to list object
        output_files.append(output_file)

    # return the list of output files
    return output_files

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculate GIA errors based on full covariance matrix
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
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
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
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
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=pathlib.Path,
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--mascon-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for mascon files')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # 1: mass coefficients
    # 2: geoid coefficients
    parser.add_argument('--fit-method',
        type=int, default=1, choices=(1,2),
        help='Method for fitting sensitivity kernel to harmonics')
    # least squares solver
    choices = ('inv','lstsq','gelsd', 'gelsy', 'gelss')
    parser.add_argument('--solver','-s',
        type=str, default='lstsq', choices=choices,
        help='Least squares solver for sensitivity kernel solutions')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass')
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
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run gia_covariance_errors_caron algorithm with parameters
        output_files = gia_covariance_errors_caron(
            args.directory,
            args.lmax,
            args.radius,
            LMIN=args.lmin,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            MASCON_FILE=args.mascon_file,
            MASCON_FORMAT=args.mascon_format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            FIT_METHOD=args.fit_method,
            SOLVER=args.solver,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
