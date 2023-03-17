#!/usr/bin/env python
u"""
quick_mascon_regress.py
Written by Tyler Sutterley (01/2023)
Creates a regression summary file for a mascon time series file

COMMAND LINE OPTIONS:
    -h, --help: Lists the command line options
    -H X, --header X: Number of rows of header text to skip
    --order X: regression fit polynomial order
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
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 05/2022: use argparse descriptions within documentation
    Written 10/2021
"""
from __future__ import print_function

import sys
import os
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: Creates regression summary files for mascon files
def run_regress(input_file,header=0,order=None,cycles=None,
    units=None,stream=False):
    """
    Creates a regression summary file for a mascon time series file

    Arguments
    ---------
    input_file: mascon time series file

    Keyword arguments
    -----------------
    header: Number of rows of header text to skip
    order: regression fit polynomial order
    cycles regression fit cyclical terms
    units: Units of input data
    stream: Stream regression summary to standard output
    """

    # Setting output parameters for each fit type
    order_str = 'x{0:d}'.format(order)
    coef_str = ['x{0:d}'.format(o) for o in range(order+1)]
    unit_suffix = [' yr^{0:d}'.format(-o) if o else '' for o in range(order+1)]
    # filename strings for cyclical terms
    cyclic_str = {}
    cyclic_str['SEMI'] = ['SS','SC']
    cyclic_str['ANN'] = ['AS','AC']
    cyclic_str['S2'] = ['S2S','S2C']
    # unit longnames for cyclical terms
    cyclic_longname = {}
    cyclic_longname['SEMI'] = ['Semi-Annual Sine', 'Semi-Annual Cosine']
    cyclic_longname['ANN'] = ['Annual Sine', 'Annual Cosine']
    cyclic_longname['S2'] = ['S2 Tidal Alias Sine', 'S2 Tidal Alias Cosine']
    amp_str = []
    for i,c in enumerate(cycles):
        if (c == 0.5):
            flag = 'SEMI'
        elif (c == 1.0):
            flag = 'ANN'
        elif (c == (161.0/365.25)):
            flag = 'S2'
        coef_str.extend(cyclic_str[flag])
        unit_suffix.extend(['',''])
        amp_str.append(flag)
    # Fitting seasonal components
    ncomp = len(coef_str)
    ncycles = 2*len(cycles)

    # stream output or print to file
    if stream:
        fid = sys.stdout
    else:
        # input file directory and basename
        ddir = os.path.dirname(input_file)
        fileBasename,_ = os.path.splitext(os.path.basename(input_file))
        # open the output regression summary file
        regress_file = '{0}_x{1}_SUMMARY.txt'.format(fileBasename,order_str[-1])
        fid = open(os.path.join(ddir, regress_file), mode='w', encoding='utf8')

    # read data
    dinput = np.loadtxt(os.path.expanduser(input_file), skiprows=header)
    # calculate regression
    fit = gravtk.time_series.regress(dinput[:,1], dinput[:,2],
        ORDER=order, CYCLES=cycles)

    # Print output to regression summary file
    # Summary filename
    print('{0}: {1}'.format('Input File',os.path.basename(input_file)), file=fid)
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
    for i,flag in enumerate(amp_str):
        # indice pointing to the cyclical components
        j = 1 + order + 2*i
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
    # command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Mascon data files')
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of rows of header text to skip')
    # regression parameters
    # 0: mean
    # 1: trend
    # 2: acceleration
    parser.add_argument('--order',
        type=int, default=2,
        help='Regression fit polynomial order')
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
                    cycles=args.cycles, units=args.units,
                    stream=args.stream)

# run main program
if __name__ == '__main__':
    main()
