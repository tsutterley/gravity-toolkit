#!/usr/bin/env python
u"""
calc_SLR_RMS.py (05/2023)
Reads low-degree zonal harmonics from Satellite Laser Ranging (SLR)
and estimates the uncertainty as the RMS off the monthly mean field

CALLING SEQUENCE:
    python calc_SLR_RMS.py --start 4 --end 232 --format netCDF4

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory with geocenter files
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    --missing X: Missing GRACE months in time series
    -F X, --format X: Output data format
        ascii
        netCDF4
        HDF5
    -V, --verbose: Verbose output of run
    -M X, --mode X: Permission mode of output files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 01/2023: refactored satellite laser ranging read functions
    Updated 12/2022: single implicit import of gravity toolkit
        use GFZ GRACE fields to include uncertainties additional harmonics
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Written 11/2021
"""
from __future__ import print_function, division

import sys
import logging
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: read SLR low-degree coefficients and calculate RMS
def calc_SLR_RMS(base_dir, START_MON, END_MON, MISSING, DATAFORM=None,
    MODE=0o775):

    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # GRACE/GRACE-FO mission gap
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    missing = sorted(set(MISSING) | set(GAP))
    # CSR 5x5 monthly harmonics
    SLR_file = base_dir.joinpath('CSR_Monthly_5x5_Gravity_Harmonics.txt')
    CSR55 = gravtk.read_SLR_harmonics(SLR_file, HEADER=True)
    # converting from MJD into month, day and year to calculate GRACE month
    YY,MM,*_ = gravtk.time.convert_julian(CSR55['MJD'] + 2400000.5,
        format='tuple')
    CSR55['month'] = gravtk.time.calendar_to_grace(YY,month=MM)
    CSR55['month'] = gravtk.time.adjust_months(CSR55['month'])
    # CSR TN-11 coefficients
    # SLR_file = base_dir.joinpath('TN-11_C20_SLR.txt')
    SLR_file = base_dir.joinpath('C20_RL06.txt')
    CSR_C20 = gravtk.SLR.C20(SLR_file)
    # GSFC TN-14 coefficients
    SLR_file = base_dir.joinpath('TN-14_C30_C20_GSFC_SLR.txt')
    GSFC_C20 = gravtk.SLR.C20(SLR_file)
    GSFC_C30 = gravtk.SLR.C30(SLR_file)
    SLR_file = base_dir.joinpath('gsfc_slr_5x5c61s61.txt')
    GSFC_CS21 = gravtk.SLR.CS2(SLR_file, ORDER=1, DATE=CSR55['time'])
    GSFC_CS22 = gravtk.SLR.CS2(SLR_file, ORDER=2, DATE=CSR55['time'])
    GSFC_C40 = gravtk.SLR.C40(SLR_file, DATE=CSR55['time'])
    GSFC_C50 = gravtk.SLR.C50(SLR_file, DATE=CSR55['time'])
    # GFZ GravIS coefficients
    SLR_file = base_dir.joinpath('GFZ_RL06_C20_SLR.dat')
    GFZ_C20 = gravtk.SLR.C20(SLR_file)
    SLR_file = base_dir.joinpath('GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat')
    GFZ_C30 = gravtk.SLR.C30(SLR_file)
    GFZ_CS21 = gravtk.SLR.CS2(SLR_file)
    # GFZ GRACE coefficients to d/o 5
    GFZ55 = gravtk.grace_input_months(base_dir, 'GFZ', 'RL06', 'GSM', 5,
        START_MON, END_MON, missing, None, None)

    # calculate common months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    common_months = np.copy(month)
    for v in [CSR55,GSFC_C30,GFZ_C30]:
        common_months = sorted(set(common_months) & set(v['month']))
    # number of common months
    nt = len(common_months)

    # calculate mean harmonics from SLR fields
    mean_Ylms = gravtk.harmonics().zeros(lmax=5, mmax=5, nt=nt)
    mean_Ylms.month[:] = np.copy(common_months)
    # calculate variance off mean harmonics
    variance_Ylms = gravtk.harmonics().zeros(lmax=5, mmax=5, nt=nt)
    variance_Ylms.month[:] = np.copy(common_months)
    # for each set of parameters
    parameters = [('clm',2,0),('clm',3,0),
                  ('clm',4,0),('clm',5,0),
                  ('clm',2,1),('slm',2,1),
                  ('clm',2,2),('slm',2,2)]
    # parameters for degree and order
    for a, (cs, l, m) in enumerate(parameters):
        if (a == 0):
            # C20
            centers = [CSR_C20,GSFC_C20,GFZ_C20]
        elif (a == 1):
            # C30
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            centers = [CSR55,GSFC_C30,GFZ_C30]
        elif (a == 2):
            # C40
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GFZ55['data'] = GFZ55[cs][l,m,:].copy()
            centers = [CSR55,GSFC_C40,GFZ55]
        elif (a == 3):
            # C50
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GFZ55['data'] = GFZ55[cs][l,m,:].copy()
            centers = [CSR55,GSFC_C50,GFZ55]
        elif (a == 4):
            # C21
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GSFC_CS21['data'] = GSFC_CS21['C2m'].copy()
            GFZ_CS21['data'] = GFZ_CS21['C2m'].copy()
            centers = [CSR55,GSFC_CS21,GFZ_CS21]
        elif (a == 5):
            # S21
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GSFC_CS21['data'] = GSFC_CS21['S2m'].copy()
            GFZ_CS21['data'] = GFZ_CS21['S2m'].copy()
            centers = [CSR55,GFZ_CS21,GFZ_CS21]
        elif (a == 6):
            # C22
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GSFC_CS22['data'] = GSFC_CS22['C2m'].copy()
            GFZ55['data'] = GFZ55[cs][l,m,:].copy()
            centers = [CSR55,GSFC_CS22,GFZ55]
        elif (a == 7):
            # S22
            CSR55['data'] = CSR55[cs][l,m,:].copy()
            GSFC_CS22['data'] = GSFC_CS22['S2m'].copy()
            GFZ55['data'] = GFZ55[cs][l,m,:].copy()
            centers = [CSR55,GSFC_CS22,GFZ55]

        # calculate the mean field for degree and order
        for i,v in enumerate(centers):
            ii = [i for i,m in enumerate(v['month']) if m in common_months]
            tmp = v['data'][ii]-v['data'][ii].mean()
            mean_Ylms.time[:] = v['time'][ii].copy()
            if (cs == 'clm'):
                mean_Ylms.clm[l,m,:] += tmp
            elif (cs == 'slm'):
                mean_Ylms.slm[l,m,:] += tmp

        # calculate variance off mean harmonics for degree and order
        for i,v in enumerate(centers):
            ii = [i for i,m in enumerate(v['month']) if m in common_months]
            tmp = v['data'][ii]-v['data'][ii].mean()
            variance_Ylms.time[:] = v['time'][ii].copy()
            if (cs == 'clm'):
                variance_Ylms.clm[l,m,:] += (tmp - mean_Ylms.clm[l,m,:]/3.0)**2
            elif (cs == 'slm'):
                variance_Ylms.slm[l,m,:] += (tmp - mean_Ylms.slm[l,m,:]/3.0)**2

    # calculate mean of harmonics
    mean_Ylms = mean_Ylms.scale(1.0/3.0)
    variance_Ylms = variance_Ylms.scale(1.0/3.0).power(0.5)
    # calculate RMS
    RMS_Ylms = gravtk.harmonics().zeros(lmax=5, mmax=5)
    RMS_Ylms.time = np.mean(variance_Ylms.time)
    RMS_Ylms.month = len(common_months)
    for Ylms in variance_Ylms:
        RMS_Ylms.add(Ylms.power(2.0))
    # convert from variance to RMS
    RMS_Ylms = RMS_Ylms.scale(1.0/nt).power(0.5)
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # output RMS harmonics to file
    args = (START_MON, END_MON, suffix[DATAFORM])
    FILE = base_dir.joinpath('SLR_RMS_{0:03d}-{1:03d}.{2}'.format(*args))
    RMS_Ylms.to_file(FILE, format=DATAFORM, **attributes)
    # change the permissions mode
    FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads low-degree zonal harmonics from
            Satellite Laser Ranging (SLR) and estimates theuncertainty
            as the RMS off the monthly mean field
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month for time series')
    parser.add_argument('--end','-E',
        type=int, default=231,
        help='Ending GRACE/GRACE-FO month for time series')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,172,
        177,178,182,200,201]
    parser.add_argument('--missing',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months in time series')
    # output data format
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=('ascii','netCDF4','HDF5'),
        help='Output data format')
    # print information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of output files')
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

    # run program with parameters
    calc_SLR_RMS(args.directory, args.start, args.end, args.missing,
        DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
