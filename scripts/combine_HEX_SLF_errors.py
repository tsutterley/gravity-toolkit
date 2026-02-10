#!/usr/bin/env python
u"""
combine_HEX_SLF_errors.py (04/2024)
calculates the estimated SLF uncertainty for spherical cap mascon regions

before running this program, you need to run:
    make_sea_level_error_shells.py with the estimated mascon satellite errors
    run_sea_level_equation.py with the Monte Carlo simulations
    sea_level_error.py with the resultant sea level fingerprints
    least_squares_mascons.py with the harmonics created in the error program

COMMAND LINE OPTIONS:
    -D X, --directory X: Working directory for reanalysis data
    -O X, --output-directory X: output directory for mascon files
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    --redistribute-mascons: redistribute mascon mass over the ocean
    --redistribute-removed: redistribute removed mass fields over the ocean
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

UPDATE HISTORY:
    Updated 04/2024: increase precision of output mass and thickness data
    Updated 02/2024: verify shape of data from input mascon files
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: output filenames similar to other combine programs
    Updated 10/2020: use argparse to set command line parameters
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: add three ice shelf caps for Antarctica
    Updated 08/2019: can set the base data directory as a command line option
        adjusted output filename to fit a more common format
    Written 04/2018
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

# Spherical Cap Parameters
RAD_CAP = 1.5

# Regions
region = []
region.extend(['AAp', 'ApB', 'BC', 'CCp', 'CpD', 'DDp', 'DpE', 'EEp', 'EpFp',
    'FpG', 'GH', 'HHp', 'HpI', 'IIpp', 'IppJ', 'JJpp', 'JppK','KKp', 'KpA',
    'INTERIOR', 'AIS', 'EAIS', 'WAIS', 'APIS', 'QML', 'NoQML', 'CpDc','CpDi',
    'DDpc','DDpi','GH2','HHp2','JJpp2','GH3','TM','PIG','THSPK'])
region.extend(['NW','NN','NE','SW','SE','GIS','CDE','CBI','ICL','SVB','ALK',
    'DEN','PNW','PAT','FJL','SZEM','NZEM','RUS','ARC','GIC'])

# Cap numbers for each basin
cap = {}
# the contents of the cap variable will be referenced to the basin name
# the cap numbers are the global numbers for the spherical cap RADius
cap['AAp'] = np.array([25,26,39,40,41,53,54])
cap['ApB'] = np.array([55,56,67,68,69,70,83,97])
cap['BC'] = np.array([80,81,82,93,94,95,96,107,108,109,121,122])
cap['CCp'] = np.array([110,123,124,136,137,138,150,151])
cap['CpD'] = np.array([135,148,149,162,163,164,175,176,177,189,190])
cap['DDp'] = np.array([174,186,187,188,200,201,202])
cap['DpE'] = np.array([158,159,160,161,172,173,185,199])
cap['EEp'] = np.array([106,117,118,119,120,131,132,133,134,145,146,147])
cap['EpFp'] = np.array([102,103,115,116,128,129,130,142,143])
cap['FpG'] = np.array([127,141])
cap['GH'] = np.array([86,100,101,113,114])
cap['HHp'] = np.array([72,85,99])
cap['HpI'] = np.array([45,59])
cap['IIpp'] = np.array([4,17,18,31,32])
cap['IppJ'] = np.array([46])
cap['JJpp'] = np.array([60,73,74,75,87,88,89])
cap['JppK'] = np.array([49,50,51,52,62,63,64,65,66,76,77,78,79,90,91,92,104,105])
cap['KKp'] = np.array([23,36,37])
cap['KpA'] = np.array([24,38])
cap['SHELVES'] = np.array([61,144,157])
# Grouping Basins into larger regions (AIS, EAIS, WAIS and APEN)
cap['AIS'] = np.concatenate((cap['AAp'], cap['ApB'], cap['BC'], cap['CCp'], \
    cap['CpD'], cap['DDp'], cap['DpE'], cap['EEp'], cap['EpFp'], cap['FpG'], \
    cap['GH'], cap['HHp'], cap['HpI'], cap['IIpp'], cap['IppJ'], cap['JJpp'], \
    cap['JppK'], cap['KKp'], cap['KpA']),axis=0) #, cap['SHELVES']
cap['EAIS'] = np.concatenate((cap['AAp'], cap['ApB'], cap['BC'], cap['CCp'], \
    cap['CpD'], cap['DDp'], cap['DpE'], cap['EEp'], cap['JppK'], cap['KKp'], \
    cap['KpA']),axis=0)
cap['WAIS'] = np.concatenate((cap['EpFp'], cap['FpG'], cap['GH'], cap['HHp'], cap['JJpp']),axis=0)
cap['APIS'] = np.concatenate((cap['HpI'], cap['IIpp'], cap['IppJ']),axis=0)
cap['INTERIOR'] = np.array([38,51,52,53,65,66,67,68,78,79,80,81,91,92,93,94,\
    104,105,106,107,108,118,119,120,121,122,123,133,134,135,136,147,148,149,\
    160,161,162,174,175])
cap['QML'] = np.concatenate((cap['KpA'],cap['AAp'],cap['ApB']),axis=0)
cap['NoQML'] = np.concatenate((cap['BC'], cap['CCp'], \
    cap['CpD'], cap['DDp'], cap['DpE'], cap['EEp'], cap['JppK'], cap['KKp']),axis=0)
cap['ISLAND'] = np.array([62])
cap['CpDc'] = np.array([163,164,176,177,189,190])
#cap['CpDc'] = np.array([163,164,177])
cap['CpDi'] = np.array([135,148,149,162,175])
cap['TM'] = np.array([135,148,149,162,163,176,177,190])
cap['DDpc'] = np.array([186,187,188,200,201,202])
#cap['DDpi'] = np.array([174])
cap['DDpi'] = np.concatenate((cap['DDp'],cap['DpE']),axis=0)
# different version of west ant
cap['GH2'] = np.array([86,87,99,100,101,113,114])
cap['HHp2'] = np.array([72,85])
cap['JJpp2'] = np.array([60,73,74,75,88,89])
cap['GH3'] = np.array([86,99,100,101,113,114])
# Amundsen Sea Embayment regions
cap['PIG'] = np.array([86,87,99])
cap['THSPK'] = np.array([100,101,113,114])

# Greenland Mascons
cap['NW'] = np.array([308,309,312,313,316,317])
cap['NN'] = np.array([301,302,303,304,305,306])
cap['NE'] = np.array([307,310,311,314,315])
cap['SE'] = np.array([318,319,322,323,325,327])
cap['SW'] = np.array([320,321,324,326])
cap['GIS'] = np.concatenate((cap['NW'],cap['NN'],cap['NE'],cap['SW'],cap['SE']),axis=0)
# Canadian Archipelago
cap['CDE'] = np.array([328,329,330,331,332,333])
cap['CBI'] = np.array([334,335,336,337,338])
# Iceland, Svalbard, Franz Josef Land, Svernaya Zemlya and Novaya Zemlya
cap['ICL'] = np.array([339])
cap['SVB'] = np.array([340])
cap['FJL'] = np.array([341])
cap['SZEM'] = np.array([342])
cap['NZEM'] = np.array([343,344])
# Alaska (Denali and Pacific Northwest)
cap['ALK'] = np.arange(345,368)
cap['DEN'] = np.array([345,346,348,349,350,351,352,353,354,355,356,367])
cap['PNW'] = np.array([347,357,358,359,360,361,362,363,364,365,366])
# Patagonia
cap['PAT'] = np.arange(400,407)

# All Arctic and all GIC (with patagonia)
cap['ARC'] = np.concatenate((cap['GIS'],cap['CDE'],cap['CBI'],cap['ICL'],
    cap['SVB'],cap['ALK'],cap['FJL'],cap['SZEM'],cap['NZEM']),axis=0)
cap['GIC'] = np.concatenate((cap['CDE'],cap['CBI'],cap['ICL'],cap['SVB'],
    cap['ALK'],cap['FJL'],cap['SZEM'],cap['NZEM'],cap['PAT']),axis=0)
# Russian Arctic (Franz Josef Land, Svernaya Zemlya and Novaya Zemlya)
cap['RUS'] = np.concatenate((cap['FJL'],cap['SZEM'],cap['NZEM']),axis=0)
# All caps (Greenland, Glaciers and Ice Caps, Antarctica)
cap['ALL'] = np.concatenate((cap['GIS'],cap['GIC'],cap['AIS']),axis=0)

# PURPOSE: keep track of threads
def info(args):
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

def combine_HEX_SLF_errors(LMAX, RAD,
    MMAX=None,
    DESTRIPE=False,
    REDISTRIBUTE_MASCONS=False,
    REDISTRIBUTE_REMOVED=False,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # Gaussian smoothing string for radius RAD (if 0: no flag)
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # input/output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # filtered (destriped) GRACE coefficients flag
    ds_str = '_FL' if DESTRIPE else ''
    # Read Ocean function and convert to Ylms for redistribution
    ocean_str = '_OCN' if (REDISTRIBUTE_MASCONS | REDISTRIBUTE_REMOVED) else ''

    # read each mascon file and store in dictionary with k of the cap number
    SLF_data = {}
    SLF_reg = {}
    # allocate for total area of region
    area_reg = {}
    # Read each spherical cap
    for k in cap['ALL']:
        # GRACE files
        a1=(RAD_CAP,k,LMAX,order_str,gw_str,ds_str,ocean_str)
        f1='MC_SPH_CAP_RAD{0:0.1f}_{1:d}_L{2:d}{3}{4}{5}{6}.txt'.format(*a1)
        # read cap time-series and verify shape
        SLF_data[k] = np.loadtxt(OUTPUT_DIRECTORY.joinpath(f1), ndmin=2)

    # monte carlo information
    RUNS = len(SLF_data[k][:,0])

    # for each region
    for i in region:
        # regional time-series
        SLF_reg[i] = np.zeros((RUNS))
        # calculate total area to convert to cm w.e.
        area_reg[i] = 0.0
        # sum mascons in region
        for k in cap[i]:
            # sum of SLF uncertainty data
            SLF_reg[i] += SLF_data[k][:,0]**2
            # add mascon area to total area (cm^2)
            area_reg[i] += 1e10*SLF_data[k][0,1]

        # output data files
        args = (i,RAD_CAP,ocean_str,LMAX,order_str,gw_str,ds_str)
        FILE1='MC_{0}_SPH_CAP_RAD{1:0.1f}{2}_L{3:d}{4}{5}{6}.txt'.format(*args)
        FILE2='MC_{0}_SPH_CAP_RAD{1:0.1f}{2}_L{3:d}{4}{5}{6}_cmwe.txt'.format(*args)
        # open files for writing region time-series
        mass_file = OUTPUT_DIRECTORY.joinpath(FILE1)
        cmwe_file = OUTPUT_DIRECTORY.joinpath(FILE2)
        mass_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
        fid1 = mass_file.open(mode='w', encoding='utf8')
        fid2 = cmwe_file.open(mode='w', encoding='utf8')
        # calculate the RMS of the SLF for the region in mass and cm thickness
        SLF_RMS_reg = np.sqrt(SLF_reg[i])
        SLF_RMS_thick = 1e15*np.sqrt(SLF_reg[i])/area_reg[i]
        for n in range(RUNS):
            print(f'{n:05d} {SLF_RMS_reg[n]:14.6f}', file=fid1)
            print(f'{n:05d} {SLF_RMS_thick[n]:14.6f}', file=fid2)
        # close the output file
        fid1.close()
        fid2.close()
        # change the permissions mode
        mass_file.chmod(mode=MODE)
        cmwe_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates the estimated SLF uncertainty for
            spherical cap mascon regions
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    parser.add_argument('--redistribute-removed',
        default=False, action='store_true',
        help='Redistribute removed mass fields over the ocean')
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
        # run combine_HEX_SLF_errors algorithm with parameters
        combine_HEX_SLF_errors(
            args.lmax,
            args.radius,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
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
