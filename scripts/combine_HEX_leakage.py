#!/usr/bin/env python
u"""
combine_HEX_leakage.py (02/2024)
calculates the leakage component of ice mass change for spherical cap mascons

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for mascon files
    -P X, --file-prefix X: prefix string for mascon files
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
    --coordinate-file X: file with spatial coordinates of mascon centers
    -H X, --header X: Number of header lines to skip in coordinate file
    -r X, --region X: Major regions to run
        ANT: Antarctic major regions
        GRN: Greenland major regions
        GIC: Glacier and Ice Cap major regions
        ALL: Grouped major regions
    -t X, --type X: Leakage mascon type (Hole, Island)
    --redistribute-mascons: redistribute mascon mass over the ocean
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

UPDATE HISTORY:
    Updated 02/2024: verify shape of data from input mascon files
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: switch from parameter files to argparse arguments
    Updated 05/2021: use argparse to set parameters
        define int/float precision to prevent deprecation warning
    Updated 09/2019: add three ice shelf caps for Antarctica
    Updated 05/2018: output fraction of sum(leak)/sum(mass) as well
        added option type to set if running island or hole mascon configurations
        added option glaciers to include more regions (for SMB-derived leakage)
    Updated 03/2018: simplified flag setting algorithms to single lines
    Written 08/2017
"""
from __future__ import print_function

import sys
import os
import logging
import pathlib
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# Spherical Cap Parameters
RAD_CAP = 1.5

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

# groupings of major regions
major_regions = {}
major_regions['ANT'] = ['AAp','ApB','BC','CCp','CpD','DDp',
    'DpE','EEp','EpFp','FpG','GH','HHp','HpI','IIpp','IppJ','JJpp',
    'JppK','KKp','KpA','INTERIOR','AIS','EAIS','WAIS','APIS','QML',
    'NoQML','CpDc','CpDi','DDpc','DDpi','GH2','HHp2','JJpp2','GH3',
    'TM','PIG','THSPK']
major_regions['GRN'] = ['NW','NN','NE','SW','SE','GIS']
major_regions['GIC'] = ['CDE','CBI','ICL','SVB','FJL','SZEM','NZEM',
    'RUS','ALK','DEN','PNW','PAT']
major_regions['ALL'] = ['ALL','ARC','GIC']

# PURPOSE: keep track of threads
def info(args):
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

def combine_HEX_leakage(LMAX, RAD,
    MMAX=None,
    DESTRIPE=False,
    COORDINATE_FILE=None,
    HEADER=0,
    REDISTRIBUTE_MASCONS=False,
    REGION=None,
    TYPE=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
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

    # offset and scale for each mascon configuration (Hole or Island)
    if (TYPE == 'Hole'):
        offset = 0.0
        scale = 1.0
    elif (TYPE == 'Island'):
        offset = 1.0
        scale = -1.0

    # check if coordinate file exists
    COORDINATE_FILE = pathlib.Path(COORDINATE_FILE).expanduser().absolute()
    if not COORDINATE_FILE.exists():
        raise FileNotFoundError(str(COORDINATE_FILE))

    # input coordinate file (center points of each spherical cap)
    # read coordinate file for lat/lon of spherical cap centers
    coord = np.loadtxt(COORDINATE_FILE, skiprows=HEADER)
    # column 1: cap number
    num = coord[:,0].astype(int)

    # read each mascon file and store in dictionary with k of the cap number
    leakage_data = {}
    # Regional Mascon components
    leakage = {}
    leakage_reg = {}
    input_reg = {}
    # Read each spherical cap
    for k in cap['ALL']:
        # input files
        args = (FILE_PREFIX,RAD_CAP,k,LMAX,order_str,gw_str,ds_str,ocean_str)
        f1 = '{0}SPH_CAP_RAD{1:0.1f}_{2:d}_L{3:d}{4}{5}{6}{7}.txt'.format(*args)
        # read cap for mascon and verify shape
        leakage_data[k] = np.loadtxt(OUTPUT_DIRECTORY.joinpath(f1), ndmin=2)

    # major regions to run
    region = [region for reg in REGION for region in major_regions[reg]]

    # output data file with all regions
    args = (FILE_PREFIX,RAD_CAP,LMAX,order_str,gw_str,ds_str,ocean_str)
    f2 = '{0}SPH_CAP_RAD{1:0.1f}_L{2:d}{3}{4}{5}{6}.txt'.format(*args)
    # open files for writing region leakage mass
    mass_file = OUTPUT_DIRECTORY.joinpath(f2)
    mass_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
    fid = mass_file.open(mode='w', encoding='utf8')
    # for each region
    for i in region:
        # regional time-series
        leakage[i] = 0.0
        leakage_reg[i] = 0.0
        input_reg[i] = 0.0
        # sum mascons in region
        for k in cap[i]:
            # find cap order
            ii, = np.nonzero(num == k)
            if (leakage_data[k][ii,2] != 0):
                leakage[i] += (leakage_data[k][ii,0]/leakage_data[k][ii,2])**2
            leakage_reg[i] += leakage_data[k][ii,0]**2
            input_reg[i] += leakage_data[k][ii,2]**2
        # print to file
        leakage_percent, = offset + scale*np.sqrt(leakage[i]/len(cap[i]))
        leakage_fraction, = offset + scale*np.sqrt(leakage_reg[i]/input_reg[i])
        print(f'{i} {leakage_percent:12.8f} {leakage_fraction:12.8f}', file=fid)
    # close the output file
    fid.close()
    mass_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""calculates the leakage component of ice mass change
            errors for spherical cap mascons
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for mascon files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for mascon files')
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
    # mascon coordinate file and parameters
    parser.add_argument('--coordinate-file',
        type=pathlib.Path,
        required=True,
        help='File with spatial coordinates of mascon centers')
    # number of header lines to skip in coordinate file
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines to skip in coordinate file')
    # major regions to run
    choices = ['ANT','GRN','GIC','ALL']
    parser.add_argument('--region','-r',
        type=lambda x: str(x).upper(),
        metavar='REGION', nargs='+',
        default=['ANT','GRN'], choices=choices,
        help='Major regions to run')
    # leakage mascon type
    parser.add_argument('--type','-t',
        type=lambda x: str(x).title(),
        metavar='TYPE', required=True,
        choices=['Hole','Island'],
        help='Leakage mascon type')
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
        # run combine_HEX_leakage algorithm with parameters
        combine_HEX_leakage(
            args.lmax,
            args.radius,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            COORDINATE_FILE=args.coordinate_file,
            HEADER=args.header,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            REGION=args.region,
            TYPE=args.type,
            OUTPUT_DIRECTORY=args.output_directory,
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