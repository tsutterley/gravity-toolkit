#!/usr/bin/env python
u"""
combine_HEX_TWC_errors.py (04/2024)
calculates the time-series of terrestrial mass leakage for spherical cap mascons
    using outputs from least_squares_mascon_timeseries.py

before running this program, you need to run:
    least_squares_mascon_timeseries.py for each model of GLDAS data

INPUTS:
    GLDAS land surface model
        CLM: Common Land Model (CLM)
        CLSM: Catchment Land Surface Model (CLSM)
        MOS: Mosaic model
        NOAH: Noah model
        VIC: Variable Infiltration Capacity (VIC) model

COMMAND LINE OPTIONS:
    -D X, --directory X: Working directory for GLDAS data
    -O X, --output-directory X: output directory for mascon files
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    --version X: GLDAS model version
    --spacing X: spatial resolution of GLDAS model
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
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
    Updated 05/2023: fix output subdirectory when changing permissions mode
        use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: output filenames similar to other combine programs
    Updated 02/2020: calculate RMS using version 2.1 of NOAH, CLSM and VIC
        remove 2003-2014 time series mean
    Updated 10/2019: calculate RMS from 4 models (NOAH, CLM, MOS and VIC)
        changing Y/N flags to True/False
    Updated 09/2019: add three ice shelf caps for Antarctica
    Updated 08/2019: adjusted output filename to fit a more common format
    Updated 03/2018: simplified flag setting algorithms to single lines
    Written 03/2018
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
    'DDpc','DDpi','GH2','HHp2','JJpp2','GH3','TM'])
region.extend(['NW','NN','NE','SW','SE','GIS','CDE','CBI','ICL','SVB','ALK',
    'DEN','PNW','PAT','FJL','SZEM','NZEM','RUS','ARC','GIC'])

# Cap numbers for each basin
cap = {}
# the contents of the cap variable will be referenced to the basin name
# the cap numbers are the global numbers for the spherical cap radius
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
cap['NZEM'] = np.array([343])
# Alaska (Denali and Pacific Northwest)
cap['ALK'] = np.arange(344,364)
cap['DEN'] = np.array([344,347,348,349,350,351,352,353,354,355])
cap['PNW'] = np.array([345,346,356,357,358,359,360,361,362,363])
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

def combine_HEX_TWC_errors(base_dir, MODEL, LMAX, RAD,
    START=None,
    END=None,
    VERSION=None,
    SPACING=None,
    MMAX=None,
    DESTRIPE=False,
    REDISTRIBUTE_MASCONS=False,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # Gaussian smoothing string for radius RAD
    gw_str = f'_r{RAD:0.0f}km' if (RAD > 0) else ''
    # input/output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # set flag if using filtered spherical harmonic coefficients
    ds_str = '_FL' if DESTRIPE else ''
    # using coefficients that were redistributed over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''
    # models to run
    if (VERSION == '1'):
        # GLDAS Version 1
        V1,V2 = ('_V1','')
    else:
        # GLDAS Version 2.1
        V1,V2 = ('','V{0}_'.format(VERSION))

    # subdirectory with output RMS data
    a = ('HEX_GLDAS_TWC_',V2,LMAX,START,END)
    sd = '{0}{1}RMS_SPH_CAP_MSCNS_L{2:d}_{3:03d}-{4:03d}'.format(*a)

    # read each mascon file and store in dictionary with k of the cap number
    TWC = {}
    TWC_mass = {}
    for M in MODEL:
        TWC[M] = {}
        TWC_mass[M] = {}
    # terrestrial water content anomalies over regions
    TWC_reg = {}
    # allocate for total area of region
    area_reg = {}

    # input files from least_squares_mascon_timeseries.py
    ff='GLDAS_{0}{1}_TWC_SPH_CAP_RAD{2:0.1f}_{3}_L{4:d}{5}{6}{7}{8}.txt'
    for m in MODEL:
        # subdirectory with data
        a=(m,SPACING,V1,LMAX,START,END)
        s='GLDAS_{0}{1}{2}_TWC_SPH_CAP_MSCNS_L{3:d}_{4:03d}-{5:03d}'.format(*a)
        # Read each spherical cap
        for k in cap['ALL']:
            # read cap time-series for mascon and verify shape
            a=(m,SPACING,RAD_CAP,k,LMAX,order_str,gw_str,ds_str,ocean_str)
            TWC[m][k] = np.loadtxt(base_dir.joinpath(s,ff.format(*a)), ndmin=2)

    # GRACE/GRACE-FO months
    mon = np.arange(START,END+1)
    TWC_months = TWC[m][k][:,0].astype(np.int64)
    missing = sorted(set(np.arange(START,END+1))-set(TWC_months))
    n_mon = len(mon)
    # GRACE/GRACE-FO dates
    calendar_year = 2002 + (mon-1)//12
    calendar_month = np.mod(mon-1,12) + 1
    tdec = gravtk.time.convert_calendar_decimal(calendar_year,calendar_month)
    # remove mean of 2003--2014
    m0314, = np.nonzero((mon >= 13) & (mon <= 156))

    # for each region
    for i in region:
        # regional time-series
        for m in MODEL:
            TWC_mass[m][i] = np.zeros((n_mon))
        # calculate total area to convert to cm w.e.
        area_reg[i] = 0.0
        # sum mascons in region
        for k in cap[i]:
            for m in MODEL:
                TWC_months = TWC[M][k][:,0].astype(np.int64)
                ind = np.ravel([np.flatnonzero(mon == m) for m in TWC_months])
                TWC_mass[m][i][ind] += TWC[m][k][:,2] - TWC[m][k][m0314,2].mean()
            # add mascon area to total area (cm^2)
            area_reg[i] += 1e10*TWC[m][k][0,3]

        # calculate mean between all hydrological models
        TWC_mean = np.zeros((n_mon))
        for m in MODEL:
            TWC_mean += TWC_mass[m][i]
        # convert from total to mean
        TWC_mean /= len(MODEL)
        # calculate variance off mean
        TWC_variance = np.zeros((n_mon))
        for m in MODEL:
            TWC_variance += (TWC_mass[m][i]-TWC_mean)**2
        # calculate RMS of mean differences
        TWC_reg[i] = np.sqrt(TWC_variance/(len(MODEL)-1.0))
        # replace invalid values with nan
        ind = np.ravel([np.flatnonzero(mon == m) for m in missing])
        TWC_reg[i][ind] = np.nan

        # output data files
        a = ('GLDAS_TWC_RMS_',i,RAD_CAP,ocean_str,LMAX,order_str,gw_str,ds_str)
        FILE1 = '{0}{1}_SPH_CAP_RAD{2:0.1f}{3}_L{4:d}{5}{6}.txt'.format(*a)
        FILE2 = '{0}{1}_SPH_CAP_RAD{2:0.1f}{3}_L{4:d}{5}{6}_cmwe.txt'.format(*a)
        # open files for writing region time-series
        mass_file = OUTPUT_DIRECTORY.joinpath(FILE1)
        cmwe_file = OUTPUT_DIRECTORY.joinpath(FILE2)
        mass_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
        fid1 = mass_file.open(mode='w', encoding='utf8')
        fid2 = cmwe_file.open(mode='w', encoding='utf8')
        # output regional mascon time-series
        for t in range(n_mon):
            mass_error = TWC_reg[i][t]
            thick_error = 1e15*TWC_reg[i][t]/area_reg[i]
            print(f'{mon[t]:03d} {tdec[t]:12.4f} {mass_error:14.6f}', file=fid1)
            print(f'{mon[t]:03d} {tdec[t]:12.4f} {thick_error:14.6f}', file=fid2)
        # close the output file
        fid1.close()
        fid2.close()
        # change the permissions mode
        mass_file.chmod(mode=MODE)
        cmwe_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""calculates the time-series of terrestrial water storage
            mass leakage for spherical cap mascons
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['CLM','CLSM','MOS','NOAH','VIC']
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['CLSM','NOAH','VIC'], choices=choices,
        help='GLDAS land surface model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    # GLDAS model version
    parser.add_argument('--version',
        type=str, default='2.1',
        help='GLDAS model version')
    # model spatial resolution
    # 10: 1.0 degrees latitude/longitude
    # 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of GLDAS model')
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
        # run combine_HEX_TWC_errors algorithm with parameters
        combine_HEX_TWC_errors(
            args.directory,
            args.model,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            VERSION=args.version,
            SPACING=args.spacing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
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
