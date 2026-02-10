#!/usr/bin/env python
u"""
combine_sea_level_data.py
Written by Tyler Sutterley (05/2023)
Combines the sea level fingerprint data with the input load harmonics

INPUTS:
    infile: index file with spherical harmonic data

COMMAND LINE OPTIONS:
    -O X, --output-directory X: output directory for file
    -P X, --file-prefix: prefix string for input and output files
    --mask X: Land-sea mask for sea level fingerprints
    -F X, --format X: Input and output data format
    -l X, --lmax X: Maximum spherical harmonic degree and order
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
    -V, --verbose: verbose output of processing run
    -M X, --mode X: permissions mode of the output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: netCDF4: Python interface to the netCDF C library
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://h5py.org)

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gen_stokes.py: Computes geoid Stokes coefficients for an input grid
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
        iterate over harmonics object
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: added path to default land-sea mask from sea level equation
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py
    Updated 12/2020: added more love number options
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: updates to reading load love numbers
        using harmonics class to read and write spherical harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 11/2018: can vary the land-sea mask for ascii files
    Written 09/2018
"""
from __future__ import print_function

import sys
import os
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

# PURPOSE: combine the sea level data with the input dataset
def combine_sea_level_data(index_file,
    LMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    LANDMASK=None,
    DIRECTORY=None,
    FILE_PREFIX=None,
    MODE=0o775):

    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]
    # input and output file format
    file_format = '{0}{1}L{2:d}_{3:03d}.{4}'

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    # degree spacing and grid dimensions
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape
    # longitude and colatitude in radians
    th = (90.0 - np.squeeze(landsea.lat))*np.pi/180.0

    # Calculating Legendre Polynomials using Holmes and Featherstone relation
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))
    # read load Love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # input spherical harmonic datafile index
    data_Ylms = gravtk.harmonics().from_index(index_file, format=DATAFORM)
    # truncate to degree and order
    data_Ylms.truncate(lmax=LMAX,mmax=LMAX)

    # index file listing all output spherical harmonic files
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    output_index_file = DIRECTORY.joinpath('index.txt')
    fid = output_index_file.open(mode='w', encoding='utf8')
    # print the path to the index file
    logging.info(str(output_index_file))

    # for each grace month and input file
    for Ylms in data_Ylms:
        # read sea level file
        SLF = file_format.format(FILE_PREFIX,'',LMAX,Ylms.month,suffix)
        input_file = DIRECTORY.joinpath(SLF)
        if (DATAFORM == 'ascii'):
            dinput = gravtk.spatial(spacing=[dlon,dlat],
                nlon=nlon, nlat=nlat).from_ascii(input_file, date=False)
        elif (DATAFORM == 'netCDF4'):
            dinput = gravtk.spatial().from_netCDF4(input_file, date=False)
        elif (DATAFORM == 'HDF5'):
            dinput = gravtk.spatial().from_HDF5(input_file, date=False)

        # Converting sea level field into spherical harmonics
        slf_Ylms = gravtk.gen_stokes(dinput.data, dinput.lon, dinput.lat,
            UNITS=1, LMIN=0, LMAX=LMAX, PLM=PLM, LOVE=LOVE)
        # add sea level harmonics to input harmonics
        Ylms.add(slf_Ylms)

        # attributes for output files
        attributes = {}
        attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
        # output combined sea level harmonics
        fargs = (FILE_PREFIX,'CLM_',LMAX,Ylms.month,suffix)
        output_file = DIRECTORY.joinpath(file_format.format(*fargs))
        Ylms.to_file(output_file, format=DATAFORM, **attributes)
        # change output file permissions mode to MODE
        output_file.chmod(mode=MODE)
        # add to output index file
        print(Ylms.compressuser(output_file), file=fid)
    # close the index file
    fid.close()
    # change the permissions mode of the output index file
    output_index_file.chmod(mode=MODE)

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
    parser.add_argument('infile',
        type=pathlib.Path,
        help='Input index file with spherical harmonic data files')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for input and output files')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # land-sea mask
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for sea level fingerprints')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
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

    # create logger for verbosity level
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run calc_sensitivity_kernel algorithm with parameters
        combine_sea_level_data(args.infile,
            LMAX=args.lmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            LANDMASK=args.mask,
            DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
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
