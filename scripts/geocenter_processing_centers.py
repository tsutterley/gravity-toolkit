#!/usr/bin/env python
u"""
geocenter_processing_centers.py
Written by Tyler Sutterley (05/2022)
Plots the GRACE/GRACE-FO geocenter time series for different
    GRACE/GRACE-FO processing centers

CALLING SEQUENCE:
    python geocenter_processing_centers.py --start 4 --end 216

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory with geocenter files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within documentation
    Updated 12/2021: adjust minimum x limit based on starting GRACE month
        make the list of processing centers an option
    Updated 11/2021: use gravity_toolkit geocenter class for operations
    Updated 05/2021: additionally plot GFZ with pole tide replaced with SLR
    Updated 04/2021: reload the matplotlib font manager
        use GRACE/GRACE-FO months to update the ticks
    Updated 02/2021: using argparse to set parameters
    Updated 04/2020: use units class for setting earth parameters
    Updated 02/2020: add minor ticks and adjust x axes
    Updated 11/2019: adjust axes and set directory to full path
    Updated 09/2019: for public release of time series to references page
"""
from __future__ import print_function

import os
import argparse
import numpy as np
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.offsetbox import AnchoredText,AnchoredOffsetbox,TextArea,VPacker
from gravity_toolkit.time import convert_calendar_decimal
from gravity_toolkit.geocenter import geocenter

#-- rebuilt the matplotlib fonts and set parameters
matplotlib.font_manager._load_fontmanager()
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams['mathtext.default'] = 'regular'

#-- PURPOSE: plots the GRACE/GRACE-FO geocenter time series
def geocenter_processing_centers(grace_dir,PROC,DREL,START_MON,END_MON,MISSING):
    #-- GRACE months
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    months = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    #-- labels for each scenario
    input_flags = ['','iter','SLF_iter','SLF_iter_wSLR21','SLF_iter_wSLR21_wSLR22']
    input_labels = ['Static','Iterated','Iterated SLF']
    #-- labels for Release-6
    model_str = 'OMCT' if DREL in ('RL04','RL05') else 'MPIOM'
    #-- degree one coefficient labels
    fig_labels = ['C11','S11','C10']
    axes_labels = dict(C10='c)',C11='a)',S11='b)')
    ylabels = dict(C10='z',C11='x',S11='y')

    #-- plot colors for each dataset
    plot_colors = dict(CSR='darkorange',GFZ='darkorchid',JPL='mediumseagreen')
    plot_colors['GFZwPT'] = 'dodgerblue'
    plot_colors['GFZ+CS21'] = 'darkorchid'
    plot_colors['GFZ+CS21+CS22'] = 'darkorchid'

    #-- 3 row plot (C10, C11 and S11)
    ax = {}
    fig,(ax[0],ax[1],ax[2])=plt.subplots(num=1,ncols=3,sharey=True,figsize=(9,4))
    #-- plot geocenter estimates for each processing center
    for k,pr in enumerate(PROC):
        #-- additionally plot GFZ with SLR replaced pole tide
        if pr in ('GFZwPT','GFZ+CS21'):
            fargs = ('GFZ',DREL,model_str,input_flags[3])
        elif (pr == 'GFZ+CS21+CS22'):
            fargs = ('GFZ',DREL,model_str,input_flags[4])
        else:
            fargs = (pr,DREL,model_str,input_flags[2])
        #-- read geocenter file for processing center and model
        grace_file = '{0}_{1}_{2}_{3}.txt'.format(*fargs)
        DEG1 = geocenter().from_UCI(os.path.join(grace_dir,grace_file))
        #-- indices for mean months
        kk, = np.nonzero((DEG1.month >= START_MON) & (DEG1.month <= 176))
        DEG1.mean(apply=True, indices=kk)
        #-- setting Load Love Number (kl) to 0.021 to match Swenson et al. (2008)
        DEG1.to_cartesian(kl=0.021)
        #-- plot each coefficient
        for j,key in enumerate(fig_labels):
            #-- create a time series with nans for missing months
            tdec = np.full_like(months,np.nan,dtype=np.float64)
            data = np.full_like(months,np.nan,dtype=np.float64)
            val = getattr(DEG1, ylabels[key].upper())
            for i,m in enumerate(months):
                valid = np.count_nonzero(DEG1.month == m)
                if valid:
                    mm, = np.nonzero(DEG1.month == m)
                    tdec[i] = DEG1.time[mm]
                    data[i] = val[mm]
            #-- plot all dates
            ax[j].plot(tdec, data, color=plot_colors[pr], label=pr)

    #-- add axis labels and adjust font sizes for axis ticks
    for j,key in enumerate(fig_labels):
        #-- vertical line denoting the accelerometer shutoff
        acc = convert_calendar_decimal(2016,9,day=3,hour=12,minute=12)
        ax[j].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        #-- vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(DEG1.month == 186)
        kk, = np.flatnonzero(DEG1.month == 198)
        vs = ax[j].axvspan(DEG1.time[jj],DEG1.time[kk],
            color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (4,2)
        #-- axis label
        ax[j].set_title(ylabels[key], style='italic', fontsize=14)
        ax[j].add_artist(AnchoredText(axes_labels[key], pad=0.,
            prop=dict(size=16,weight='bold'), frameon=False, loc=2))
        ax[j].set_xlabel('Time [Yr]', fontsize=14)
        #-- set ticks
        xmin = 2002 + (START_MON + 1.0)//12.0
        xmax = 2002 + (END_MON + 1.0)/12.0
        major_ticks = np.arange(2005, xmax, 5)
        ax[j].xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(xmin, xmax, 1)) - set(major_ticks))
        ax[j].xaxis.set_ticks(minor_ticks, minor=True)
        ax[j].set_xlim(xmin, xmax)
        ax[j].set_ylim(-9.5,8.5)
        #-- axes tick adjustments
        ax[j].get_xaxis().set_tick_params(which='both', direction='in')
        ax[j].get_yaxis().set_tick_params(which='both', direction='in')
        for tick in ax[j].xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax[j].yaxis.get_major_ticks():
            tick.label.set_fontsize(14)

    #-- add legend
    lgd = ax[0].legend(loc=3,frameon=False)
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    for i,text in enumerate(lgd.get_texts()):
        text.set_weight('bold')
        text.set_color(plot_colors[text.get_text()])
    #-- labels and set limits
    ax[0].set_ylabel('Geocenter Variation [mm]', fontsize=14)
    #-- adjust locations of subplots
    fig.subplots_adjust(left=0.06,right=0.98,bottom=0.12,top=0.94,wspace=0.05)
    #-- save figure to file
    OUTPUT_FIGURE = 'SV19_{0}_centers.pdf'.format(DREL)
    plt.savefig(os.path.join(grace_dir,OUTPUT_FIGURE), format='pdf', dpi=300)
    plt.clf()

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Plots the GRACE/GRACE-FO geocenter time series for
            different GRACE/GRACE-FO processing centers
            """
    )
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Data processing center or satellite mission
    PROC = ['CSR','GFZ','GFZwPT','JPL']
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, nargs='+', default=PROC,
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str,
        default='RL06', choices=['RL04','RL05','RL06'],
        help='GRACE/GRACE-FO data release')
    #-- start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month for time series')
    parser.add_argument('--end','-E',
        type=int, default=230,
        help='Ending GRACE/GRACE-FO month for time series')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,172,
        177,178,182,200,201]
    parser.add_argument('--missing','-M',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months in time series')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- run program with parameters
    geocenter_processing_centers(args.directory, args.center, args.release,
        args.start, args.end, args.missing)

#-- run main program
if __name__ == '__main__':
    main()
