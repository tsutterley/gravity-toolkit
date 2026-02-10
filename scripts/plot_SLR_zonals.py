#!/usr/bin/env python
u"""
plot_SLR_zonals.py (05/2023)
Plots low-degree zonal harmonics from GRACE/GRACE-FO and SLR
Compares with a climatology calculated using GRACE months
    with both accelerometers operational

CALLING SEQUENCE:
    python plot_SLR_zonals.py --start 4 --end 216

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory with geocenter files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: place matplotlib import behind try/except statement
    Updated 01/2023: refactored satellite laser ranging read functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 07/2022: set plot tick formatter to not use offsets
    Updated 05/2022: use argparse descriptions within documentation
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 05/2021: use argparse to set parameters
    Forked 10/2019 from plot_slr_oblateness.py
    Updated 08/2019: updated for GRACE RL06
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 10/2017: update axis plot ticks to be inward
    Updated 04/2017: using __future__ print function.  minor code updates
    Updated 03/2017: output legend similar to matplotlib 1.0
    Written 05/2016
"""
from __future__ import print_function, division

import inspect
import pathlib
import warnings
import argparse
import numpy as np
import gravity_toolkit as gravtk

# attempt imports
try:
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['mathtext.default'] = 'regular'
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
    import matplotlib.offsetbox
except ModuleNotFoundError:
    warnings.warn("matplotlib not available", ImportWarning)

# current file path for the child programs
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# plot SLR low-degree zonal coefficients
def plot_SLR_zonals(base_dir, PROC, DREL, START_MON, END_MON, MISSING):
    # GRACE/GRACE-FO mission gap
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    missing = sorted(set(MISSING) | set(GAP))
    # GRACE/GRACE-FO monthly harmonics
    grace_Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, 'GSM', 5,
        START_MON, END_MON, missing, None, 0, MMAX=None, MODEL_DEG1=False,
        ATM=False, POLE_TIDE=False, SLR_C30='N')
    GRACE = dict(date=grace_Ylms['time'], month=grace_Ylms['month'])
    grace_Ylms = gravtk.harmonics().from_dict(grace_Ylms)
    grace_Ylms.mean(apply=True)
    # CSR 5x5 monthly harmonics
    SLR_file = base_dir.joinpath('CSR_Monthly_5x5_Gravity_Harmonics.txt')
    C55 = gravtk.read_SLR_harmonics(SLR_file, HEADER=True)
    # GSFC 5x5 weekly harmonics
    SLR_file = base_dir.joinpath('gsfc_slr_5x5c61s61.txt')
    G55 = gravtk.read_SLR_harmonics(SLR_file, HEADER=True)
    # converting from MJD into month, day and year to calculate GRACE month
    YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(C55['MJD']+2400000.5,
        format='tuple')
    C55['month'] = np.array(12.0*(YY - 2002.0)+MM, dtype=np.int64)

    # filtered coefficients from John Ries
    SLR_file = base_dir.joinpath('C30_LARES_filtered.txt')
    LARES_C30 = gravtk.SLR.C30(SLR_file)
    # CSR TN-11 coefficients
    SLR_file = base_dir.joinpath('TN-11_C20_SLR.txt')
    CSR_C20 = gravtk.SLR.C20(SLR_file)
    # GSFC TN-14 coefficients
    SLR_file = base_dir.joinpath('TN-14_C30_C20_GSFC_SLR.txt')
    GSFC_C20 = gravtk.SLR.C20(SLR_file)
    GSFC_C30 = gravtk.SLR.C30(SLR_file)
    # GFZ GravIS coefficients
    SLR_file = base_dir.joinpath('GFZ_RL06_C20_SLR.dat')
    GFZ_C20 = gravtk.SLR.C20(SLR_file)
    SLR_file = base_dir.joinpath('GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat')
    GFZ_C30 = gravtk.SLR.C30(SLR_file)

    # calculate common months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    common_months = np.copy(month)
    for v in [GRACE,C55,LARES_C30,GSFC_C30]:
        common_months = sorted(set(common_months) & set(v['month']))

    # read arrays of kl, hl, and ll Love Numbers
    love_numbers_file = gravtk.utilities.get_data_path(['data','love_numbers'])
    hl,kl,ll = gravtk.read_love_numbers(love_numbers_file, FORMAT='tuple')
    # Earth Parameters
    factors = gravtk.units(lmax=5).harmonic(hl,kl,ll)
    rho_e = factors.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e# Average Radius of the Earth [cm]

    # create output plot
    ax1 = {}
    ax2 = {}
    f1, (ax1[2],ax1[3],ax1[4],ax1[5]) = plt.subplots(num=1,nrows=4,
        sharex=True,figsize=(6,8))
    f2, (ax2[2],ax2[3],ax2[4],ax2[5]) = plt.subplots(num=2,nrows=4,
        sharex=True,figsize=(6,8))
    plot_colors = ['mediumseagreen','0.5','darkorchid','darkorange']
    plot_labels = [f'{PROC} GRACE/GRACE-FO',f'{PROC} GRACE Climatology',
        'GSFC SLR','GFZ GravIS']#'CSR TN-11 SLR'
    plot_zorder = [0,3,1,2]
    fig_text = {2:'a)',3:'b)',4:'c)',5:'d)'}
    for l,ax in ax1.items():
        # GRACE data for zonal harmonic
        GRACE['data'] = grace_Ylms.clm[l,0,:].copy()
        # CSR and GSFC 5x5 data for zonal harmonic
        C55['data'] = C55['clm'][l,0,:].copy()
        G55['data'] = G55['clm'][l,0,:].copy()
        # calculate GRACE climatology for "good months"
        CLIMATE = {}
        YY = np.arange(2002+1./24.,2024,1./12.).astype(np.int64)
        CLIMATE['month'] = 1 + np.arange(len(YY))
        MM = 1 + ((CLIMATE['month']-1) % 12)
        CLIMATE['date'] = gravtk.time.convert_calendar_decimal(YY, MM)
        ii, = np.nonzero((GRACE['month'] >= 13) & (GRACE['month'] <= 176))
        CLIMATE['data'] = regress_model(GRACE['date'][ii], GRACE['data'][ii],
            CLIMATE['date'], ORDER=2, CYCLES=[0.25,0.5,1.0,2.0,4.0,5.0,161.0/365.25],
            RELATIVE=GRACE['date'][0])
        # plot items
        if (l == 2):
            plot_items = [GRACE,CLIMATE,GSFC_C20,GFZ_C20]
        elif (l == 3):
            plot_items = [GRACE,CLIMATE,GSFC_C30,GFZ_C30]
        else:
            GSFC = gravtk.convert_weekly(G55['time'], G55['data'],
                DATE=GRACE['date'], NEIGHBORS=28)
            plot_items = [GRACE,CLIMATE,GSFC]
        for i,v in enumerate(plot_items):
            ii = [i for i,m in enumerate(v['month']) if m in common_months]
            # create a time series with nans for missing months
            tdec = np.full_like(month,np.nan,dtype=np.float64)
            cmwe = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(v['month'] == m)
                if valid:
                    mm, = np.nonzero(CLIMATE['month'] == m)
                    tdec[d] = CLIMATE['date'][mm]
                    mm, = np.nonzero(v['month'] == m)
                    # tdec[d] = v['date'][mm]
                    # cm w.e.: centimeters water equivalent
                    dfactor = factors.cmwe[l]
                    cmwe[d] = dfactor*(v['data'][mm]-v['data'][ii].mean())
            # plot all dates
            ax.plot(tdec, cmwe, zorder=plot_zorder[i],
                color=plot_colors[i], label=plot_labels[i])
            ax2[l].plot(tdec, cmwe, zorder=plot_zorder[i],
                color=plot_colors[i], label=plot_labels[i])

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        ax2[l].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(GRACE['month'] == 186)
        kk, = np.flatnonzero(GRACE['month'] == 198)
        # ax.axvline(GRACE['time'][jj],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # ax.axvline(GRACE['time'][kk],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        vs = ax.axvspan(GRACE['date'][jj],GRACE['date'][kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        vs = ax2[l].axvspan(GRACE['date'][jj],GRACE['date'][kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # set ticks
        xmax = 2002 + (END_MON + 1.0)/12.0
        # set limits for first plot
        major_ticks = np.arange(2014,xmax,2)
        minor_ticks = sorted(set(np.arange(2002, xmax, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        ax.set_xlim(2013,xmax)
        # set limits for second plot
        major_ticks = np.arange(2002,xmax,2)
        minor_ticks = sorted(set(np.arange(2002, xmax, 1)) - set(major_ticks))
        ax2[l].xaxis.set_ticks(major_ticks)
        ax2[l].xaxis.set_ticks(minor_ticks, minor=True)
        ax2[l].set_xlim(2001.75,xmax)
        # add labels
        ax.set_ylabel('C{0:d}0 [cm]'.format(l))
        ax2[l].set_ylabel('C{0:d}0 [cm]'.format(l))
        # adjust ticks
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax2[l].get_xaxis().set_tick_params(which='both', direction='in')
        ax2[l].get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax.xaxis.get_major_formatter().set_useOffset(False)
        ax2[l].xaxis.get_major_formatter().set_useOffset(False)
        # add labels to the first plot
        at = matplotlib.offsetbox.AnchoredText(fig_text[l],
            prop=dict(size=14,weight='bold'),
            pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        # add labels to the second plot
        at = matplotlib.offsetbox.AnchoredText(fig_text[l],
            prop=dict(size=14,weight='bold'),
            pad=0, frameon=False, loc=2)
        ax2[l].add_artist(at)

    # set y limits for each axis in the first plot
    ax1[2].set_ylim(-5.6,4.6)
    ax1[3].set_ylim(-2.6,2.1)
    ax1[4].set_ylim(-2.1,1.6)
    # ax1[4].yaxis.set_ticks([-1,0,1])
    ax1[5].set_ylim(-2.1,2.1)
    # set y limits for each axis in the second plot
    ax2[2].set_ylim(-5.1,4.6)
    ax2[3].set_ylim(-2.6,2.1)
    ax2[4].set_ylim(-2.1,2.1)
    # ax2[4].yaxis.set_ticks([-1,0,1])
    ax2[5].set_ylim(-2.1,3.6)

    # add label and legend
    ax1[5].set_xlabel('Time [Yr]')
    ax2[5].set_xlabel('Time [Yr]')
    # add legend to the first plot
    lgd = ax1[2].legend(loc=3,frameon=True,ncol=2)
    # set width, color and style of lines
    lgd.get_frame().set_boxstyle('square,pad=0.01')
    lgd.get_frame().set_edgecolor('white')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    # add legend to the second plot
    lgd = ax2[2].legend(loc=3,frameon=True,ncol=2)
    # set width, color and style of lines
    lgd.get_frame().set_boxstyle('square,pad=0.1')
    lgd.get_frame().set_edgecolor('white')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)

    # adjust the plot and save to file
    f1.subplots_adjust(left=0.08,right=0.97,bottom=0.05,top=0.99,hspace=0.06)
    figurefile = filepath.joinpath('fs02ad.pdf')
    f1.savefig(figurefile, format='pdf')
    # adjust the plot and save to file
    f2.subplots_adjust(left=0.08,right=0.97,bottom=0.05,top=0.99,hspace=0.06)
    figurefile = filepath.joinpath('fs02ad_all.pdf')
    f2.savefig(figurefile, format='pdf')
    # close everything
    plt.cla()
    plt.clf()
    plt.close()

    # create a separate plot just with C30 with the CSR LARES solutions
    l = 3
    # cm w.e. degree dependent factor for centimeters water equivalent
    dfactor = factors.cmwe[l]
    # GRACE data for zonal harmonic
    GRACE['data'] = grace_Ylms.clm[l,0,:].copy()
    # CSR 5x5 data for zonal harmonic
    C55['data'] = C55['clm'][l,0,:].copy()
    # calculate GRACE climatology for "good months"
    CLIMATE = {}
    YY = np.arange(2002+1./24.,2022,1./12.).astype(np.int64)
    CLIMATE['month'] = 1 + np.arange(len(YY))
    MM = 1 + ((CLIMATE['month']-1) % 12)
    CLIMATE['date'] = gravtk.time.convert_calendar_decimal(YY, MM)
    ii, = np.nonzero((GRACE['month'] >= 13) & (GRACE['month'] <= 176))
    CLIMATE['data'] = regress_model(GRACE['date'][ii], GRACE['data'][ii],
        CLIMATE['date'], ORDER=2, CYCLES=[0.25,0.5,1.0,2.0,4.0,5.0],
        RELATIVE=GRACE['date'][0])

    # create output plot
    fig, ax = plt.subplots(num=3, figsize=(6.5,4))
    plot_colors = ['mediumseagreen','0.5','darkorchid','darkorange']
    plot_labels = ['CSR GRACE/GRACE-FO','CSR GRACE Climatology',
        'GSFC TN-14 SLR','CSR LARES SLR']
    plot_zorder = [0,3,1,2]
    for i,v in enumerate([GRACE,CLIMATE,GSFC_C30,LARES_C30]):
        ii = [i for i,m in enumerate(v['month']) if m in common_months]
        # create a time series with nans for missing months
        tdec = np.full_like(month,np.nan,dtype=np.float64)
        cmwe = np.full_like(month,np.nan,dtype=np.float64)
        for d,m in enumerate(month):
            valid = np.count_nonzero(v['month'] == m)
            if valid:
                mm, = np.nonzero(CLIMATE['month'] == m)
                tdec[d] = CLIMATE['date'][mm]
                mm, = np.nonzero(v['month'] == m)
                # tdec[d] = v['date'][mm]
                # cm w.e.: centimeters water equivalent
                cmwe[d] = dfactor*(v['data'][mm]-v['data'][ii].mean())
        # plot all dates
        ax.plot(tdec, cmwe, zorder=plot_zorder[i],
            color=plot_colors[i], label=plot_labels[i])

    # vertical line denoting the accelerometer shutoff
    acc = gravtk.time.convert_calendar_decimal(2016,9,day=3,hour=12,minute=12)
    ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
    # vertical lines for end of the GRACE mission and start of GRACE-FO
    jj, = np.flatnonzero(GRACE['month'] == 186)
    kk, = np.flatnonzero(GRACE['month'] == 198)
    # ax.axvline(GRACE['time'][jj],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
    # ax.axvline(GRACE['time'][kk],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
    vs = ax.axvspan(GRACE['date'][jj],GRACE['date'][kk],color='0.5',ls='dashed',alpha=0.15)
    vs._dashes = (6,3)
    # set limits
    major_ticks = np.arange(2010,2022,2)
    minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
    ax.xaxis.set_ticks(minor_ticks, minor=True)
    ax.set_xlim(2010,2021.5)
    ax.set_ylim(-2.6,1.6)
    # add labels
    ax.set_xlabel('Time [Yr]')
    ax.set_ylabel('C{0:d}0 [cm]'.format(l))
    # adjust ticks
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    # set axis ticker
    ax.xaxis.get_major_formatter().set_useOffset(False)
    # add legend
    lgd = ax.legend(loc=3,frameon=True,ncol=2)
    # set width, color and style of lines
    lgd.get_frame().set_boxstyle('square,pad=0.1')
    lgd.get_frame().set_edgecolor('white')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    # adjust the plot and show
    fig.subplots_adjust(left=0.1,right=0.97,bottom=0.1,top=0.97)
    figurefile = filepath.joinpath('fs02b.pdf')
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

# PURPOSE: calculate a regression model for extrapolating values
def regress_model(t_in, d_in, t_out, ORDER=2, CYCLES=None, RELATIVE=None):
    # remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    t_out = np.squeeze(t_out)
    # check dimensions of output
    if (np.ndim(t_out) == 0):
        t_out = np.array([t_out])
    # CREATING DESIGN MATRIX FOR REGRESSION
    DMAT = []
    MMAT = []
    # add polynomial orders (0=constant, 1=linear, 2=quadratic)
    for o in range(ORDER+1):
        DMAT.append((t_in-RELATIVE)**o)
        MMAT.append((t_out-RELATIVE)**o)
    # add cyclical terms (0.5=semi-annual, 1=annual)
    for c in CYCLES:
        if (c == (161.0/365.25)):
            # terms for S2 tidal aliasing during GRACE and GRACE-FO periods
            DMAT.extend(gravtk.time_series.aliasing_terms(t_in))
            MMAT.extend(gravtk.time_series.aliasing_terms(t_out))
            # remove the original S2 tidal aliasing term from CYCLES list
            CYCLES.remove(c)
        else:
            DMAT.append(np.sin(2.0*np.pi*t_in/np.float64(c)))
            DMAT.append(np.cos(2.0*np.pi*t_in/np.float64(c)))
            MMAT.append(np.sin(2.0*np.pi*t_out/np.float64(c)))
            MMAT.append(np.cos(2.0*np.pi*t_out/np.float64(c)))
    # Calculating Least-Squares Coefficients
    # Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(np.transpose(DMAT), d_in, rcond=-1)[0]
    # return modeled time-series
    return np.dot(np.transpose(MMAT),beta_mat)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Plots low-degree zonal harmonics from
            GRACE/GRACE-FO and SLR
            """
    )
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-C',
        metavar='PROC', type=str,
        default=None, choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO data processing center')
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
        type=int, default=231,
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
    plot_SLR_zonals(args.directory, args.center, args.release,
        args.start, args.end, args.missing)

# run main program
if __name__ == '__main__':
    main()
