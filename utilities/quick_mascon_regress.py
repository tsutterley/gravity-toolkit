#!/usr/bin/env python
u"""
quick_mascon_regress.py
Written by Tyler Sutterley (06/2023)
Creates a regression summary file for a mascon time series file

COMMAND LINE OPTIONS:
    -h, --help: Lists the command line options
    -H X, --header X: Number of rows of header text to skip
    --order X: regression fit polynomial order
    --breakpoint X: breakpoint GRACE month for piecewise regression
    --cycles X: regression fit cyclical terms
    -U X, --units X: Units of input data
    -S, --stream: Stream regression summary to standard output

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

PROGRAM DEPENDENCIES:
    time_series.regress.py: calculates trend coefficients using least-squares
    time_series.amplitude.py: calculates the amplitude and phase of a harmonic

UPDATE HISTORY:
    Updated 06/2023: can choose different tidal aliasing periods
    Updated 05/2023: split S2 tidal aliasing terms into GRACE and GRACE-FO eras
        allow fit to be piecewise using a known breakpoint GRACE/GRACE-FO month
        use fit module for getting tidal aliasing terms
        use pathlib to define and operate on paths
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Written 10/2021
"""
from __future__ import print_function

import sys
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: Creates regression summary files for mascon files
def run_regress(input_file,
        header=0,
        order=None,
        breakpoint=None,
        cycles=None,
        units=None,
        stream=False
    ):
    """
    Creates a regression summary file for a mascon time series file

    Arguments
    ---------
    input_file: mascon time series file

    Keyword arguments
    -----------------
    header: Number of rows of header text to skip
    order: regression fit polynomial order
    breakpoint: breakpoint GRACE/GRACE-FO month for piecewise regression
    cycles regression fit cyclical terms
    units: Units of input data
    stream: Stream regression summary to standard output
    """

    # read data
    input_file = pathlib.Path(input_file).expanduser().absolute()
    dinput = np.loadtxt(input_file, skiprows=header)

    # fitting with either piecewise or polynomial regression
    if breakpoint is not None:
        # Setting output parameters for piecewise fit
        breakpoint_index, = np.nonzero(dinput[:,0] == breakpoint)
        cycle_index = 3
        coef_str = ['x0', 'px1', 'px1']
        unit_suffix = ['', ' yr^-1', ' yr^-1']
        # output regression filename
        output_file = f'{input_file.stem}_px1_{breakpoint:d}_SUMMARY.txt'
    elif order is not None:
        # Setting output parameters for each fit type
        cycle_index = 1 + order
        coef_str = ['x{0:d}'.format(o) for o in range(order+1)]
        unit_suffix = [' yr^{0:d}'.format(-o) if o else '' for o in range(order+1)]
        # output regression filename
        output_file = f'{input_file.stem}_x{order:d}_SUMMARY.txt'
    else:
        raise ValueError('No regression fit specified')

    # stream output or print to file
    if stream:
        fid = sys.stdout
    else:
        # full path to output regression summary file
        regress_file = input_file.with_name(output_file)
        # open the output regression summary file
        fid = regress_file.open(mode='w', encoding='utf8')

    # amplitude string for cyclical components
    amp_str = []
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
    terms = []
    term_index = []
    for i,c in enumerate(cycles):
        # check if fitting with semi-annual or annual terms
        if (c == 0.5):
            coef_str.extend(['SS','SC'])
            amp_str.append('SEMI')
            unit_suffix.extend(['',''])
        elif (c == 1.0):
            coef_str.extend(['AS','AC'])
            amp_str.append('ANN')
            unit_suffix.extend(['',''])
        # check if fitting with tidal aliasing terms
        for t,period in tidal_aliasing.items():
            if np.isclose(c, (period/365.25)):
                # terms for tidal aliasing during GRACE and GRACE-FO periods
                terms.extend(gravtk.time_series.aliasing_terms(dinput[:,1],
                    period=period))
                # labels for tidal aliasing during GRACE period
                coef_str.extend([f'{t}SGRC', f'{t}CGRC'])
                amp_str.append(f'{t}GRC')
                unit_suffix.extend(['',''])
                # labels for tidal aliasing during GRACE-FO period
                coef_str.extend([f'{t}SGFO', f'{t}CGFO'])
                amp_str.append(f'{t}GFO')
                unit_suffix.extend(['',''])
                # index to remove the original tidal aliasing term
                term_index.append(i)
    # remove the original tidal aliasing terms
    cycles = np.delete(cycles, term_index)

    # calculate regression
    if breakpoint is not None:
        fit = gravtk.time_series.piecewise(dinput[:,1], dinput[:,2],
            BREAKPOINT=breakpoint_index, CYCLES=cycles, TERMS=terms)
    elif order is not None:
        fit = gravtk.time_series.regress(dinput[:,1], dinput[:,2],
            ORDER=order, CYCLES=cycles, TERMS=terms)
    # Fitting seasonal components
    ncycles = 2*len(cycles) + len(terms)

    # Print output to regression summary file
    # Summary filename
    print('{0}: {1}'.format('Input File', input_file.name), file=fid)
    # Regression Formula with Correlation Structure if Applicable
    print('Regression Formula: {0}\n'.format('+'.join(coef_str)), file=fid)
    # Value, Error and Statistical Significance
    args = ('Coef.','Estimate','Std. Error','95% Conf.','Units')
    fid.write('{0:5}\t{1:12}\t{2:12}\t{3:12}\t{4:10}\n'.format(*args))
    print(64*'-',file=fid)
    for i,b in enumerate(fit['beta']):
        args=(coef_str[i],b,fit['std_err'][i],fit['error'][i],units+unit_suffix[i])
        fid.write('{0:5}\t{1:12.4f}\t{2:12.4f}\t{3:12.4f}\t{4:10}\n'.format(*args))
    # allocate for amplitudes and phases of cyclical components
    amp,ph = ({},{})
    for comp in ['beta','std_err','error']:
        amp[comp] = np.zeros((ncycles//2))
        ph[comp] = np.zeros((ncycles//2))
    # calculate amplitudes and phases of cyclical components
    for i, flag in enumerate(amp_str):
        # indice pointing to the cyclical components
        j = cycle_index + 2*i
        amp['beta'][i],ph['beta'][i] = gravtk.time_series.amplitude(
            fit['beta'][j], fit['beta'][j+1]
        )
        # convert phase from -180:180 to 0:360
        if (ph['beta'][i] < 0):
            ph['beta'][i] += 360.0
        # calculate standard error and 95% confidences
        for err in ['std_err','error']:
            # Amplitude Errors
            comp1 = fit[err][j]*fit['beta'][j]/amp['beta'][i]
            comp2 = fit[err][j+1]*fit['beta'][j+1]/amp['beta'][i]
            amp[err][i] = np.sqrt(comp1**2 + comp2**2)
            # Phase Error (degrees)
            comp1 = fit[err][j]*fit['beta'][j+1]/(amp['beta'][i]**2)
            comp2 = fit[err][j+1]*fit['beta'][j]/(amp['beta'][i]**2)
            ph[err][i] = (180.0/np.pi)*np.sqrt(comp1**2 + comp2**2)

    # Amplitude, Error and Statistical Significance
    args = ('Ampl.','Estimate','Std. Error','95% Conf.','Units')
    fid.write('\n{0:5}\t{1:12}\t{2:12}\t{3:12}\t{4:10}\n'.format(*args))
    print(64*'-',file=fid)
    for i,b in enumerate(amp['beta']):
        args=(amp_str[i],b,amp['std_err'][i],amp['error'][i],units)
        fid.write('{0:5}\t{1:12.4f}\t{2:12.4f}\t{3:12.4f}\t{4:10}\n'.format(*args))

    # Phase, Error and Statistical Significance
    args = ('Phase','Estimate','Std. Error','95% Conf.','Units')
    fid.write('\n{0:5}\t{1:12}\t{2:12}\t{3:12}\t{4:10}\n'.format(*args))
    print(64*'-',file=fid)
    for i,b in enumerate(ph['beta']):
        args=(amp_str[i],b,ph['std_err'][i],ph['error'][i],'Degree')
        fid.write('{0:5}\t{1:12.4f}\t{2:12.4f}\t{3:12.4f}\t{4:10}\n'.format(*args))

    # Fit Significance Criteria
    print('\nFit Criterion', file=fid)
    print(64*'-',file=fid)
    fid.write('{0}: {1:d}\n'.format('DOF', fit['DOF']))
    for fitstat in ['AIC','BIC','LOGLIK','MSE','NRMSE','R2','R2Adj']:
        fid.write('{0}: {1:f}\n'.format(fitstat, fit[fitstat]))
    # close the output file
    fid.write('\n') if stream else fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates a regression summary file for a mascon
            time series file
        """
    )
    group = parser.add_mutually_exclusive_group(required=True)
    # command line parameters
    parser.add_argument('file',
        type=pathlib.Path, nargs='+',
        help='Mascon data files')
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of rows of header text to skip')
    # regression parameters
    # 0: mean
    # 1: trend
    # 2: acceleration
    group.add_argument('--order',
        type=int,
        help='Regression fit polynomial order')
    # breakpoint month for piecewise regression
    group.add_argument('--breakpoint',
        type=int,
        help='Breakpoint GRACE/GRACE-FO month for piecewise regression')
    # regression fit cyclical terms
    parser.add_argument('--cycles',
        type=float, default=[0.5,1.0,161.0/365.25], nargs='+',
        help='Regression fit cyclical terms')
    parser.add_argument('--units','-U',
        type=str, default='Gt',
        help='Units of input data')
    # stream to sys.stdout
    parser.add_argument('--stream','-S',
        default=False, action='store_true',
        help='Stream regression summary to standard output')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run plot program for each input file
    for i,f in enumerate(args.file):
        run_regress(f, header=args.header, order=args.order,
                    breakpoint=args.breakpoint, cycles=args.cycles,
                    units=args.units, stream=args.stream)

# run main program
if __name__ == '__main__':
    main()
