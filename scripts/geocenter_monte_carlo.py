#!/usr/bin/env python
u"""
geocenter_monte_carlo.py
Written by Tyler Sutterley (12/2022)

CALLING SEQUENCE:
    python geocenter_monte_carlo.py --start 4 --end 237

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory with geocenter files
    -c X, --center X: GRACE/GRACE-FO data processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series

UPDATE HISTORY:
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: adjust minimum x limit based on starting GRACE month
    Written 11/2021
"""
from __future__ import print_function

import os
import argparse
import numpy as np
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.offsetbox import AnchoredText
import gravity_toolkit as gravtk

# rebuilt the matplotlib fonts and set parameters
matplotlib.font_manager._load_fontmanager()
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams['mathtext.default'] = 'regular'

# PURPOSE: plots the GRACE/GRACE-FO geocenter time series
def geocenter_monte_carlo(grace_dir,PROC,DREL,START_MON,END_MON,MISSING):
    # GRACE months
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    months = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    nmon = len(months)
    # labels for Release-6
    model_str = 'OMCT' if DREL in ('RL04','RL05') else 'MPIOM'
    # GIA and processing labels
    input_flag = 'SLF'
    gia_str = '_AW13_ice6g_GA'
    delta_str = '_monte_carlo'
    ds_str = '_FL'

    # degree one coefficient labels
    fig_labels = ['C11','S11','C10']
    axes_labels = dict(C10='c)',C11='a)',S11='b)')
    ylabels = dict(C10='z',C11='x',S11='y')

    # 3 row plot (C10, C11 and S11)
    ax = {}
    fig,(ax[0],ax[1],ax[2])=plt.subplots(num=1,ncols=3,sharey=True,figsize=(9,4))

    # read geocenter file for processing center and model
    fargs = (PROC,DREL,model_str,input_flag,gia_str,delta_str,ds_str)
    grace_file = '{0}_{1}_{2}_{3}{4}{5}{6}.nc'.format(*fargs)
    DEG1 = gravtk.geocenter().from_netCDF4(os.path.join(grace_dir,grace_file))
    # setting Load Love Number (kl) to 0.021 to match Swenson et al. (2008)
    DEG1.to_cartesian(kl=0.021)
    # number of monte carlo runs
    _,nruns = np.shape(DEG1.C10)

    # plot each coefficient
    for j,key in enumerate(fig_labels):
        # create a time series with nans for missing months
        tdec = np.full((nmon),np.nan,dtype=np.float64)
        data = np.full((nmon,nruns),np.nan,dtype=np.float64)
        val = getattr(DEG1, ylabels[key].upper())
        for i,m in enumerate(months):
            valid = np.count_nonzero(DEG1.month == m)
            if valid:
                mm, = np.nonzero(DEG1.month == m)
                tdec[i] = DEG1.time[mm]
                data[i,:] = val[mm,:]

        # show solutions for each iteration
        plot_colors = iter(cm.rainbow(np.linspace(0,1,nruns)))
        # mean of all monte carlo solutions
        MEAN = np.mean(data, axis=1)
        nvalid = np.count_nonzero(np.isfinite(MEAN))
        # calculate variance off of the mean
        variance = np.zeros((nruns))
        max_var = 0.0
        for k in range(nruns):
            color_k = next(plot_colors)
            # plot all dates
            ax[j].plot(tdec, data[:,k], color=color_k)
            # variance off of the mean
            variance[k] = np.nansum((data[:,k] - MEAN)**2)/nvalid
            if (np.nanmax(np.abs(data[:,k] - MEAN)) > max_var):
                max_var = np.nanmax(np.abs(data[:,k] - MEAN))
        # add mean solution
        ax[j].plot(tdec, MEAN, color='k', lw=1)
        # calculate total RMS
        RMS = np.nansum(np.sqrt(variance))/nruns

        # add axis labels and adjust font sizes for axis ticks
        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,day=3,hour=12,minute=12)
        ax[j].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(DEG1.month == 186)
        kk, = np.flatnonzero(DEG1.month == 198)
        vs = ax[j].axvspan(DEG1.time[jj],DEG1.time[kk],
            color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (4,2)
        # axis label
        ax[j].set_title(ylabels[key], style='italic', fontsize=14)
        ax[j].add_artist(AnchoredText(axes_labels[key], pad=0.,
            prop=dict(size=16,weight='bold'), frameon=False, loc=2))
        lbl = '$\sigma$ = {0:0.2f} mm\nmax = {1:0.2f} mm'.format(RMS,max_var)
        ax[j].add_artist(AnchoredText(lbl, pad=0.,
            prop=dict(size=12), frameon=False, loc=3))
        ax[j].set_xlabel('Time [Yr]', fontsize=14)
        # set ticks
        xmin = 2002 + (START_MON + 1.0)//12.0
        xmax = 2002 + (END_MON + 1.0)/12.0
        major_ticks = np.arange(2005, xmax, 5)
        ax[j].xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(xmin, xmax, 1)) - set(major_ticks))
        ax[j].xaxis.set_ticks(minor_ticks, minor=True)
        ax[j].set_xlim(xmin, xmax)
        ax[j].set_ylim(-9.5,8.5)
        # axes tick adjustments
        ax[j].get_xaxis().set_tick_params(which='both', direction='in')
        ax[j].get_yaxis().set_tick_params(which='both', direction='in')
        for tick in ax[j].xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax[j].yaxis.get_major_ticks():
            tick.label.set_fontsize(14)

    # labels and set limits
    ax[0].set_ylabel(f'{PROC} Geocenter Variation [mm]', fontsize=14)
    # adjust locations of subplots
    fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.94,wspace=0.05)
    # save figure to file
    OUTPUT_FIGURE = f'SV19_{PROC}_{DREL}_monte_carlo.pdf'
    plt.savefig(os.path.join(grace_dir,OUTPUT_FIGURE), format='pdf', dpi=300)
    plt.clf()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Plots the GRACE/GRACE-FO geocenter time series for
            each iteration of a monte carlo solution
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str,
        default='RL06', choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month for time series')
    parser.add_argument('--end','-E',
        type=int, default=236,
        help='Ending GRACE/GRACE-FO month for time series')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,172,
        177,178,182,200,201]
    parser.add_argument('--missing','-M',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months in time series')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program with parameters
    geocenter_monte_carlo(args.directory, args.center, args.release,
        args.start, args.end, args.missing)

# run main program
if __name__ == '__main__':
    main()
