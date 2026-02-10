#!/usr/bin/env python
u"""
plot_SLR_azimuthal.py (05/2023)
Plots degree-two order-one harmonics from GRACE/GRACE-FO and SLR
Compares with a climatology calculated using GRACE months
    with both accelerometers operational

CALLING SEQUENCE:
    python plot_SLR_azimuthal.py --start 4 --end 216

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
    Written 05/2021
"""
from __future__ import print_function, division

import inspect
import pathlib
import argparse
import warnings
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

# plot SLR azimuthal dependence coefficients
def plot_SLR_azimuthal(base_dir, PROC, DREL, START_MON, END_MON, MISSING):
    # GRACE/GRACE-FO mission gap
    GAP = [187,188,189,190,191,192,193,194,195,196,197]
    missing = sorted(set(MISSING) | set(GAP))
    # CSR GRACE/GRACE-FO monthly harmonics
    grace_Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, 'GSM', 5,
        START_MON, END_MON, missing, None, 0, MMAX=None, MODEL_DEG1=False,
        ATM=False, POLE_TIDE=False, SLR_C30='N')
    GRACE = dict(date=grace_Ylms['time'], month=grace_Ylms['month'])
    grace_Ylms = gravtk.harmonics().from_dict(grace_Ylms)
    grace_Ylms.mean(apply=True)

    # CSR coefficients
    SLR_file = base_dir.joinpath('C21_S21_RL06.txt')
    CSR_CS2 = gravtk.SLR.CS2(SLR_file)
    # GSFC coefficients
    SLR_file = base_dir.joinpath('gsfc_slr_5x5c61s61.txt',DATE=GRACE['date'])
    GSFC_CS2 = gravtk.SLR.CS2(SLR_file)
    # GFZ GravIS coefficients
    SLR_file = base_dir.joinpath('GRAVIS-2B_GFZOP_GRACE+SLR_LOW_DEGREES_0002.dat')
    GFZ_CS2 = gravtk.SLR.CS2(SLR_file)

    # calculate common months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    common_months = np.copy(month)
    for v in [GRACE,CSR_CS2,GFZ_CS2]:
        common_months = sorted(set(common_months) & set(v['month']))

    # read arrays of kl, hl, and ll Love Numbers
    love_numbers_file = gravtk.utilities.get_data_path(['data','love_numbers'])
    hl,kl,ll = gravtk.read_love_numbers(love_numbers_file, FORMAT='tuple')
    # Earth Parameters
    factors = gravtk.units(lmax=5).harmonic(hl,kl,ll)
    rho_e = factors.rho_e# Average Density of the Earth [g/cm^3]
    rad_e = factors.rad_e# Average Radius of the Earth [cm]

    # create a separate plot with azimuthal dependence files
    l = 2
    # cm w.e. degree dependent factor for centimeters water equivalent
    dfactor = factors.cmwe[l]
    # create output plot
    ax1 = {}
    ax2 = {}
    f1, (ax1['C21'],ax1['S21']) = plt.subplots(num=1,nrows=2,
        sharex=True,figsize=(6,4))
    f2, (ax2['C21'],ax2['S21']) = plt.subplots(num=2,nrows=2,
        sharex=True,figsize=(6,4))
    plot_colors = ['mediumseagreen','0.5','darkorchid','dodgerblue','darkorange']
    plot_labels = ['CSR GRACE/GRACE-FO','CSR GRACE Climatology',
        'CSR SLR','GSFC SLR','GFZ GravIS']
    plot_zorder = [0,3,1,2]
    fig_text = {'C21':'a)','S21':'b)'}
    plot_ylabel = {'C21':'C$\mathregular{_{21}}$','S21':'S$\mathregular{_{21}}$'}
    CSR = dict(time=CSR_CS2['time'],month=CSR_CS2['month'])
    GSFC = dict(time=GSFC_CS2['time'],month=GSFC_CS2['month'])
    GFZ = dict(time=GFZ_CS2['time'],month=GFZ_CS2['month'])
    for cs,ax in ax1.items():
        # GRACE data for zonal harmonic
        if (cs == 'C21'):
            GRACE['data'] = grace_Ylms.clm[2,1,:].copy()
            # SLR azimuthal dependence data
            CSR['data'] = CSR_CS2['C2m'].copy()
            GSFC['data'] = GSFC_CS2['C2m'].copy()
            GFZ['data'] = GFZ_CS2['C2m'].copy()
        elif (cs == 'S21'):
            GRACE['data'] = grace_Ylms.slm[2,1,:].copy()
            # SLR azimuthal dependence data
            CSR['data'] = CSR_CS2['S2m'].copy()
            GSFC['data'] = GSFC_CS2['S2m'].copy()
            GFZ['data'] = GFZ_CS2['S2m'].copy()
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
        # plot items
        plot_items = [GRACE,CLIMATE,CSR,GSFC,GFZ]
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
                    cmwe[d] = dfactor*(v['data'][mm]-v['data'][ii].mean())
            # plot all dates
            ax.plot(tdec, cmwe, zorder=plot_zorder[i],
                color=plot_colors[i], label=plot_labels[i])
            ax2[cs].plot(tdec, cmwe, zorder=plot_zorder[i],
                color=plot_colors[i], label=plot_labels[i])

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        ax2[cs].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(GRACE['month'] == 186)
        kk, = np.flatnonzero(GRACE['month'] == 198)
        # ax.axvline(GRACE['time'][jj],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # ax.axvline(GRACE['time'][kk],color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        vs = ax.axvspan(GRACE['date'][jj],GRACE['date'][kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        vs = ax2[cs].axvspan(GRACE['date'][jj],GRACE['date'][kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # set ticks
        xmax = 2002 + (END_MON + 1.0)/12.0
        # set limits for first plot
        major_ticks = np.arange(2010,xmax,2)
        minor_ticks = sorted(set(np.arange(2002, xmax, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        ax.set_xlim(2010,xmax)
        # set limits for second plot
        major_ticks = np.arange(2002,xmax,2)
        minor_ticks = sorted(set(np.arange(2002, xmax, 1)) - set(major_ticks))
        ax2[cs].xaxis.set_ticks(major_ticks)
        ax2[cs].xaxis.set_ticks(minor_ticks, minor=True)
        ax2[cs].set_xlim(2001.75,xmax)
        # add labels
        ax.set_ylabel('{0} [cm]'.format(cs))
        ax2[cs].set_ylabel('{0} [cm]'.format(cs))
        # adjust ticks
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax2[cs].get_xaxis().set_tick_params(which='both', direction='in')
        ax2[cs].get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax.xaxis.get_major_formatter().set_useOffset(False)
        ax2[cs].xaxis.get_major_formatter().set_useOffset(False)
        # add labels to the first plot
        at = matplotlib.offsetbox.AnchoredText(fig_text[cs],
            prop=dict(size=14,weight='bold'),
            pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        # add labels to the second plot
        at = matplotlib.offsetbox.AnchoredText(fig_text[cs],
            prop=dict(size=14,weight='bold'),
            pad=0, frameon=False, loc=2)
        ax2[cs].add_artist(at)

    # set y limits for each axis in the first plot
    ax1['C21'].set_ylim(-1.6,1.1)
    ax1['S21'].set_ylim(-1.1,1.6)
    ax1['C21'].yaxis.set_ticks([-1,0,1])
    ax1['S21'].yaxis.set_ticks([-1,0,1])
    # set y limits for each axis in the second plot
    ax2['C21'].set_ylim(-2.1,2.1)
    ax2['S21'].set_ylim(-1.6,1.6)

    # add label and legend
    ax1['S21'].set_xlabel('Time [Yr]')
    ax2['S21'].set_xlabel('Time [Yr]')
    # add legend to the first plot
    lgd = ax1['C21'].legend(loc=1,frameon=True,ncol=2)
    # set width, color and style of lines
    lgd.get_frame().set_boxstyle('square,pad=0.01')
    lgd.get_frame().set_edgecolor('white')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    # add legend to the second plot
    lgd = ax2['C21'].legend(loc=1,frameon=True,ncol=2)
    # set width, color and style of lines
    lgd.get_frame().set_boxstyle('square,pad=0.1')
    lgd.get_frame().set_edgecolor('white')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)

    # adjust the plot and save to file
    f1.subplots_adjust(left=0.08,right=0.97,bottom=0.05,top=0.99,hspace=0.06)
    figurefile = filepath.joinpath('fs03ab.pdf')
    f1.savefig(figurefile, format='pdf')
    # adjust the plot and save to file
    f2.subplots_adjust(left=0.08,right=0.97,bottom=0.05,top=0.99,hspace=0.06)
    figurefile = filepath.joinpath('fs03ab_all.pdf')
    f2.savefig(figurefile, format='pdf')
    # close everything
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
        description="""Plots degree-two order-one harmonics
            from GRACE/GRACE-FO and SLR
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
    plot_SLR_azimuthal(args.directory, args.center, args.release,
        args.start, args.end, args.missing)

# run main program
if __name__ == '__main__':
    main()
