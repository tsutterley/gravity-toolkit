#!/usr/bin/env python
u"""
plot_mascon_SLF_timeseries.py
Written by Tyler Sutterley (05/2023)

CALLING SEQUENCE:
    python plot_mascon_SLF_timeseries.py --start 4 --end 230

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series
    -S X, --smb X: Plot Surface Mass Balance (SMB) time series

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: place imports behind try/except statements
        updated GEMB SMB time series to version 1.2
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: add JPL GEMB SMB time series for comparison
        make smb option an iterable list of models
    Updated 10/2022: add RACMO and GSFC-fdm time series for comparison
    Updated 07/2022: set plot tick formatter to not use offsets
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: added AW13 combination IJ05-R2 and ICE6G models
    Updated 05/2021: new TMB leakage directory path
        define int/float precision to prevent deprecation warning
    Updated 04/2021: using argparse to set parameters
    Updated 01/2021: updated for new open-source processing scheme
    Written 11/2019
"""
from __future__ import print_function

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
try:
    import scipy.stats
    import scipy.special
except ModuleNotFoundError:
    warnings.warn("scipy not available", ImportWarning)

# current file path for the child programs
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: plot mascon time series from Velicogna et al. (2014)
def plot_mascon_SLF_timeseries(base_dir,PROC,DREL,START_MON,END_MON,MISSING,SMB=False):
    # directory setup
    mascon_dir = base_dir.joinpath('GRACE','mascons')
    # GIA parameters
    GIA = {}
    gia_files = {}

    # Simpson (2009), Peltier (2015,2018), A (ICE6G), Caron (2018)
    GIA['N'] = ['SM09','ICE6G','ICE6G-D','AW13-ICE6G','Caron']
    # Ivins (2013), Whitehouse (2012), Peltier (2015,2018), A (ICE6G/IJ05), Caron (2018)
    GIA['S'] = ['IJ05-R2','W12a','ICE6G','ICE6G-D','AW13-ICE6G','ascii','Caron']
    # GIA files for each modeling group
    gia_files['AW13-ICE6G'] = []
    gia_files['AW13-ICE6G'].append(['GIA','AW13','ICE6G','stokes.ice6g_.1_10.'])
    gia_files['AW13-ICE6G'].append(['GIA','AW13','ICE6G','stokes.ice6g_1._10.'])
    gia_files['AW13-ICE6G'].append(['GIA','AW13','ICE6G','stokes.ice6g_.1_1.'])
    gia_files['AW13-ICE6G'].append(['GIA','AW13','ICE6G','stokes.ice6g_1._1.'])
    gia_files['AW13-ICE6G'].append(['GIA','AW13','ICE6G','stokes.ice6g_GA.txt'])
    gia_files['ascii'] = []
    gia_files['ascii'].append(['GIA','AW13','IJ05-R2','IJ05_R2_115_.2_1.5_ICE6G.txt'])
    gia_files['ascii'].append(['GIA','AW13','IJ05-R2','IJ05_R2_65_.2_1.5_ICE6G.txt'])
    gia_files['Caron'] = []
    gia_files['Caron'].append(['GIA','Caron','expStokes_GIA.txt'])
    gia_files['ICE6G-D'] = []
    gia_files['ICE6G-D'].append(['GIA','ICE6G','VersionD','Stokes_trend_High_Res.txt'])
    gia_files['ICE6G-D'].append(['GIA','ICE6G','VersionD','Stokes_trend_VM5a_O512.txt'])
    gia_files['ICE6G'] = []
    gia_files['ICE6G'].append(['GIA','ICE6G','VM5','Stokes_G_Rot_60_I6_A_VM5a'])
    gia_files['ICE6G'].append(['GIA','ICE6G','VM5','Stokes_G_Rot_60_I6_A_VM5b'])
    gia_files['IJ05-R2'] = []
    gia_files['IJ05-R2'].append(['GIA','IJ05-R2','Stokes.R2_115_.2_1.5_L120'])
    gia_files['IJ05-R2'].append(['GIA','IJ05-R2','Stokes.R2_115_.2_2._L120'])
    gia_files['IJ05-R2'].append(['GIA','IJ05-R2','Stokes.R2_115_.2_3.2_L120'])
    gia_files['IJ05-R2'].append(['GIA','IJ05-R2','Stokes.R2_115_.2_4._L120'])
    gia_files['IJ05-R2'].append(['GIA','IJ05-R2','Stokes.R2_65_.2_1.5_L120'])
    gia_files['SM09'] = []
    gia_files['SM09'].append(['GIA','SM09','grate_120p11.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_120p51.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_120p53.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_120p81.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_96p32.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_96p510.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_96p55.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_96p58.clm'])
    gia_files['SM09'].append(['GIA','SM09','grate_96p85.clm'])
    gia_files['W12a'] = []
    gia_files['W12a'].append(['GIA','W12a','grate_B.clm'])
    gia_files['W12a'].append(['GIA','W12a','grate_L.clm'])
    gia_files['W12a'].append(['GIA','W12a','grate_U.clm'])

    # Setting output error alpha with confidence interval
    CONF = 0.95
    alpha = 1.0 - CONF

    # GRACE data release and months
    GAP = [187,188,189,190,191,192,193,194,195,196,197,]
    months = sorted(set(np.arange(START_MON,END_MON+1))-set(MISSING)-set(GAP))
    nmon = len(months)
    # start and end GRACE months for correction data
    OBP_START,OBP_END = (4,254)
    ATM_START,ATM_END = (4,251)
    GLDAS_START,GLDAS_END = (4,254)

    # atmospheric ECMWF "jump" flag and ocean redistribution flag
    atm_str = 'wATM_' if (DREL == 'RL05') else ''
    ocean_str = 'OCN_'
    # maximum degree and order
    LMAX = 60
    # gaussian smoothing radius
    RAD = 250
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # destripe string
    ds_str = ''
    # create a set of months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    # version flags
    VERSION = ['v0','']
    plot_title = ['Version 0','Version 1']
    plot_colors = ['darkorchid','mediumseagreen','darkorange','red','dodgerblue']

    # subdirectory and input file formats
    sd = 'HEX_{0}_{1}{2}_SPH_CAP_MSCNS{3}_L{4:d}_{5:03d}-{6:03d}'
    ff = '{0}_{1}{2}_SPH_CAP_{3}{4}L{5:d}{6}{7}.txt'

    # create figure axis for Greenland plots
    ax1 = {}
    ax2 = {}
    fig, ((ax1['GIS'], ax1['NN'], ax1['NE']), (ax1['NW'], ax1['SW'], ax1['SE'])) = \
        plt.subplots(num=2, nrows=2, ncols=3, sharex=False, sharey=False, figsize=(12, 7.5))
    ylimits = {'NW':[-1200,1100,200],'NN':[-450,450,100],'NE':[-180,320,40],\
        'SW':[-650,550,100],'SE':[-900,1200,200],'GIS':[-2800,3000,400]}
    fig_text = {'GIS':'a)','NN':'b)','NE':'c)','NW':'d)','SW':'e)','SE':'f)'}
    reg_text = {'GIS':'GIS','NN':'N','NE':'NE','NW':'NW','SW':'SW','SE':'SE'}
    h = 'N'
    SLF = '_SLF3'
    RACMO_START,RACMO_END = (4,239)
    GEMB_START,GEMB_END = (4,251)
    ypad = dict(GIS=4,NN=4,NE=4,NW=4,SW=4,SE=4)
    for reg,ax in ax1.items():
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str)
        # ATM_file = ff.format('ATM_Differences',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        ATM_file = ff.format('ATM_Differences',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        ATM_input = np.loadtxt(mascon_dir.joinpath(subdir,ATM_file))[:nmon,:]
        # read GLDAS terrestrial water RMS file
        subdir = sd.format('GLDAS','TWC_V2.1_RMS','','',LMAX,GLDAS_START,GLDAS_END)
        TWC_file = ff.format('GLDAS_TWC_RMS',reg,'','RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
        TWC_input = np.loadtxt(base_dir.joinpath('GLDAS',subdir,TWC_file))[:nmon,:]
        isvalid, = np.nonzero(np.isfinite(TWC_input[:,2]) &
            (TWC_input[:,0] >= START_MON) & (TWC_input[:,0] <= END_MON))
        TWC_RMS = np.sqrt(np.sum(TWC_input[isvalid,2]**2)/len(isvalid))
        # input estimated SLF monte carlo variance file and calculate RMS
        subdir = sd.format(PROC,DREL,'','_MC',LMAX,START_MON,END_MON)
        SLF_file = ff.format('MC',reg,'','RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
        SLF_input = np.loadtxt(mascon_dir.joinpath(subdir,SLF_file))
        SLF_RMS = np.sqrt(np.sum(SLF_input[:,1]**2)/len(SLF_input))
        # calculate mean of RMS over multiple reanalyses
        OBP_RMS = 0.0
        ATM_RMS = 0.0
        # for j in range(2):
        #     ivalid, = np.nonzero(np.isfinite(OBP_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(OBP_input[:,j+3]))
        #     OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(OBP_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(OBP_input[:,2]))
        OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,2]**2)/valid_count)
        # for j in range(4):
        #     ivalid, = np.nonzero(np.isfinite(ATM_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(ATM_input[:,j+3]))
        #     ATM_RMS += np.sqrt(np.sum(OBP_input[ATM_input,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(ATM_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(ATM_input[:,2]))
        ATM_RMS += np.sqrt(np.sum(ATM_input[ivalid,2]**2)/valid_count)
        # # divide by the number of reanalyses
        # OBP_RMS /= 2.0
        # ATM_RMS /= 4.0
        # GIA modeling group for region
        g = GIA[h][0]
        # number of rheologies to iterate
        nRheology = len(gia_files[g])
        # iterate through solutions
        mon = np.zeros((nmon),dtype=np.int64)
        tdec = np.zeros((nmon),dtype=np.float64)
        mass = np.zeros((nmon,nRheology),dtype=np.float64)
        satellite_error = np.zeros((nmon),dtype=np.float64)
        for i,F in enumerate(VERSION):
            # subdirectory
            subdir = sd.format(PROC,DREL,F,SLF,LMAX,START_MON,END_MON)
            # read each GIA model
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_Ylms = gravtk.read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
                gia_str = gia_Ylms['title']
                input_file = ff.format(gia_str,reg,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                dinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                mon[:] = dinput[:nmon,0].astype(np.int64)
                tdec[:] = dinput[:nmon,1]
                mass[:,k] = dinput[:nmon,2]
                satellite_error[:] += dinput[:nmon,3]**2
            # calculate mean GIA-corrected mass change (for all Earth rheologies)
            gia_corrected_mean = np.mean(mass, axis=1)
            # GRACE satellite error component, ocean leakage (ECCO-GAD),
            # atmosphere leakage (ATM-GAA) and GLDAS TWC
            grace_error = np.sqrt(np.sum(satellite_error/nRheology + SLF_RMS**2 +
                OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon)
            # calculate variance off of mean for calculating GIA uncertainty
            gia_corrected_variance = np.zeros((nmon))
            gia_corrected_minmax = np.zeros((nmon))
            # calculate GIA uncertainty as "worst-case" not RMS
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_corrected_variance+=np.abs(mass[:,k]-gia_corrected_mean)
            # calculate GIA uncertainty as "worst-case" min max error
            for t in range(nmon):
                gia_corrected_minmax[t]=np.abs(np.max(mass[t,:])-np.min(mass[t,:]))
            # calculate uncertainty in mean GIA
            gia_corrected_error = gia_corrected_variance/(np.float64(nRheology)-1.0)
            gia_corrected_error = gia_corrected_error*np.sign(tdec-tdec.mean())
            gia_corrected_minmax = gia_corrected_minmax*np.sign(tdec-tdec.mean())
            gia_error_rate = (gia_corrected_error[-1]-gia_corrected_error[0])/(tdec[-1]-tdec[0])
            gia_minmax_rate = (gia_corrected_minmax[-1]-gia_corrected_minmax[0])/(tdec[-1]-tdec[0])
            # calculate GIA errors at confidence interval
            # t.ppf parallels tinv in matlab
            tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nRheology-1.0) if g in ('IJ05-R2','SM09') else 2.0
            gia_corrected_conf = tstar*np.abs(gia_error_rate)
            # add to plot with colors and label
            plot_label = plot_title[i]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(mon == m)
                if valid:
                    mm, = np.nonzero(mon == m)
                    tnan[d] = tdec[mm]
                    mnan[d] = gia_corrected_mean[mm] - np.mean(gia_corrected_mean)
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[i], label=plot_label, zorder=3+i)
            # fill between monthly errors
            ax.fill_between(tnan, mnan-grace_error, y2=mnan+grace_error,
                color=plot_colors[i], alpha=0.35, zorder=1+i)

        # plot RACMO SMB time series for comparison
        if ('RACMO' in SMB):
            # read RACMO surface mass balance file
            subdir = sd.format('RACMO2.3p2','FGRN055_DS1km_v4.0_SMB_cumul','','',LMAX,RACMO_START,RACMO_END)
            RACMO_file = ff.format('RACMO2.3p2','FGRN055_DS1km_v4.0_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            RACMO = np.loadtxt(base_dir.joinpath('RACMO','SMB1km_v4.0',subdir,RACMO_file))[:nmon,:]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(RACMO[:,0] == m)
                if valid:
                    mm, = np.nonzero(RACMO[:,0] == m)
                    tnan[d] = RACMO[mm,1]
                    mnan[d] = RACMO[mm,2] - np.mean(RACMO[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[2], label='RACMO2.3p2', zorder=5)

        # plot GSFC-fdm SMB time series for comparison
        if ('GSFC-fdm' in SMB):
            # read MERRA-2 hybrid surface mass balance file
            subdir = 'HEX_GSFC_FDM_{0}_{1}_{2}_L{3:d}'.format('v1_2_1','gris','SMB_a',LMAX)
            MERRA2_file = ff.format('GSFC_FDM','v1_2_1_SMB_a_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            MERRA2 = np.loadtxt(base_dir.joinpath('MERRA2_hybrid','v1.2.1',subdir,MERRA2_file))
            # plot all dates
            ax.plot(MERRA2[:,1], MERRA2[:,2] - np.mean(MERRA2[:,2]),
                color=plot_colors[3], label='GSFC-fdm v1.2.1', zorder=5)

        # plot JPL-GEMB SMB time series for comparison
        if ('GEMB' in SMB):
            # read GEMB surface mass balance file
            subdir = sd.format('GEMB','v1_2_Greenland_SMB_cumul','','',LMAX,GEMB_START,GEMB_END)
            GEMB_file = ff.format('GEMB','v1_2_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            GEMB = np.loadtxt(base_dir.joinpath('GEMB','v1.2',subdir,GEMB_file))
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(GEMB[:,0] == m)
                if valid:
                    mm, = np.nonzero(GEMB[:,0] == m)
                    tnan[d] = GEMB[mm,1]
                    mnan[d] = GEMB[mm,2] - np.mean(GEMB[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[4], label='GEMB v1.2', zorder=5)

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax.axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax.axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax.axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[reg],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(reg_text[reg],
            prop=textprops, pad=0, frameon=False, loc=1)
        ax.add_artist(at)
        # set ticks
        major_ticks = np.arange(2002, 2024, 4)
        ax.xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        axlim = ax.set_xlim(2002, 2023.0)
        axlim = ax.set_ylim(ylimits[reg][0], ylimits[reg][1])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax.xaxis.get_major_formatter().set_useOffset(False)
        # add x and y labels
        ax.set_xlabel('Time [Yr]')
        ax.set_ylabel('Mass [Gt]', labelpad=ypad[reg])

    # add legend
    lgd = ax1['GIS'].legend(loc=3,frameon=False)
    # lgd = ax1['GIS'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
    # set width, color and style of lines
    # lgd.get_frame().set_boxstyle('square,pad=0.1')
    # lgd.get_frame().set_edgecolor('black')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    #     line.set_linewidth(0)
    for i,text in enumerate(lgd.get_texts()):
        text.set_color(plot_colors[i])
        text.set_weight('bold')

    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.0625,right=0.99,bottom=0.05,top=0.99,wspace=0.2,hspace=0.125)
    figurefile = filepath.joinpath('fig3af_{0}_{1}.pdf'.format(PROC,DREL))
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

    # create figure axis for Antarctic plots
    ax1 = {}
    ax2 = {}
    fig, ((ax1['AIS'], ax1['APIS'], ax1['GH3']), (ax1['QML'], ax1['CpDc'], ax1['DDpi'])) = \
       plt.subplots(num=2, nrows=2, ncols=3, sharex=False, sharey=False, figsize=(12, 7.5))
    ylimits = {'AIS':[-1700,1600,250],'APIS':[-350,350,50],'GH3':[-1400,1400,200],\
        'QML':[-750,850,100],'CpDc':[-400,300,50],'DDpi':[-200,250,50]}
    fig_text = {'AIS':'a)','APIS':'b)','GH3':'c)','QML':'d)','CpDc':'e)','DDpi':'f)'}
    reg_text = {'AIS':'AIS','APIS':'APIS','GH3':'ASE','QML':'QML','CpDc':'TMF','DDpi':'VW'}
    h = 'S'
    SLF = '_SLF3'
    RACMO_START,RACMO_END = (4,243)
    GEMB_START,GEMB_END = (4,251)
    ypad = dict(AIS=0,APIS=4,GH3=0,QML=4,CpDc=4,DDpi=4)
    for reg,ax in ax1.items():
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str)
        # ATM_file = ff.format('ATM_Differences',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        ATM_file = ff.format('ATM_Differences',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        ATM_input = np.loadtxt(mascon_dir.joinpath(subdir,ATM_file))[:nmon,:]
        # read GLDAS terrestrial water RMS file
        subdir = sd.format('GLDAS','TWC_V2.1_RMS','','',LMAX,GLDAS_START,GLDAS_END)
        TWC_file = ff.format('GLDAS_TWC_RMS',reg,'','RAD1.5_','',LMAX,gw_str,ds_str)
        TWC_input = np.loadtxt(base_dir.joinpath('GLDAS',subdir,TWC_file))[:nmon,:]
        isvalid, = np.nonzero(np.isfinite(TWC_input[:,2]) &
            (TWC_input[:,0] >= START_MON) & (TWC_input[:,0] <= END_MON))
        TWC_RMS = np.sqrt(np.sum(TWC_input[isvalid,2]**2)/len(isvalid))
        # input estimated SLF monte carlo variance file and calculate RMS
        subdir = sd.format(PROC,DREL,'','_MC',LMAX,START_MON,END_MON)
        SLF_file = ff.format('MC',reg,'','RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
        SLF_input = np.loadtxt(mascon_dir.joinpath(subdir,SLF_file))
        SLF_RMS = np.sqrt(np.sum(SLF_input[:,1]**2)/len(SLF_input))
        # calculate mean of RMS over multiple reanalyses
        OBP_RMS = 0.0
        ATM_RMS = 0.0
        # for j in range(2):
        #     ivalid, = np.nonzero(np.isfinite(OBP_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(OBP_input[:,j+3]))
        #     OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(OBP_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(OBP_input[:,2]))
        OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,2]**2)/valid_count)
        # for j in range(4):
        #     ivalid, = np.nonzero(np.isfinite(ATM_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(ATM_input[:,j+3]))
        #     ATM_RMS += np.sqrt(np.sum(OBP_input[ATM_input,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(ATM_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(ATM_input[:,2]))
        ATM_RMS += np.sqrt(np.sum(ATM_input[ivalid,2]**2)/valid_count)
        # # divide by the number of reanalyses
        # OBP_RMS /= 2.0
        # ATM_RMS /= 4.0
        # GIA modeling group for region
        g = GIA[h][0]
        # number of rheologies to iterate
        nRheology = len(gia_files[g])
        # iterate through solutions
        mon = np.zeros((nmon),dtype=np.int64)
        tdec = np.zeros((nmon),dtype=np.float64)
        mass = np.zeros((nmon,nRheology),dtype=np.float64)
        satellite_error = np.zeros((nmon),dtype=np.float64)
        for i,F in enumerate(VERSION):
            # subdirectory
            subdir = sd.format(PROC,DREL,F,SLF,LMAX,START_MON,END_MON)
            # read each GIA model
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_Ylms = gravtk.read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
                gia_str = gia_Ylms['title']
                input_file = ff.format(gia_str,reg,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                dinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                mon[:] = dinput[:nmon,0].astype(np.int64)
                tdec[:] = dinput[:nmon,1]
                mass[:,k] = dinput[:nmon,2]
                satellite_error[:] += dinput[:nmon,3]**2
            # calculate mean GIA-corrected mass change (for all Earth rheologies)
            gia_corrected_mean = np.mean(mass, axis=1)
            # GRACE satellite error component, ocean leakage (ECCO-GAD),
            # atmosphere leakage (ATM-GAA) and GLDAS TWC
            grace_error = np.sqrt(np.sum(satellite_error/nRheology + SLF_RMS**2 +
                OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon)
            # calculate variance off of mean for calculating GIA uncertainty
            gia_corrected_variance = np.zeros((nmon))
            gia_corrected_minmax = np.zeros((nmon))
            # calculate GIA uncertainty as "worst-case" not RMS
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_corrected_variance+=np.abs(mass[:,k]-gia_corrected_mean)
            # calculate GIA uncertainty as "worst-case" min max error
            for t in range(nmon):
                gia_corrected_minmax[t]=np.abs(np.max(mass[t,:])-np.min(mass[t,:]))
            # calculate uncertainty in mean GIA
            gia_corrected_error = gia_corrected_variance/(np.float64(nRheology)-1.0)
            gia_corrected_error = gia_corrected_error*np.sign(tdec-tdec.mean())
            gia_corrected_minmax = gia_corrected_minmax*np.sign(tdec-tdec.mean())
            gia_error_rate = (gia_corrected_error[-1]-gia_corrected_error[0])/(tdec[-1]-tdec[0])
            gia_minmax_rate = (gia_corrected_minmax[-1]-gia_corrected_minmax[0])/(tdec[-1]-tdec[0])
            # calculate GIA errors at confidence interval
            # t.ppf parallels tinv in matlab
            tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nRheology-1.0) if g in ('IJ05-R2','SM09') else 2.0
            gia_corrected_conf = tstar*np.abs(gia_error_rate)
            # add to plot with colors and label
            plot_label = plot_title[i]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(mon == m)
                if valid:
                    mm, = np.nonzero(mon == m)
                    tnan[d] = tdec[mm]
                    mnan[d] = gia_corrected_mean[mm] - np.mean(gia_corrected_mean)
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[i], label=plot_label, zorder=3+i)
            # fill between monthly errors
            ax.fill_between(tnan, mnan-grace_error, y2=mnan+grace_error,
                color=plot_colors[i], alpha=0.35, zorder=1+i)

        # plot RACMO SMB time series for comparison
        if ('RACMO' in SMB):
            # read RACMO surface mass balance file
            subdir = sd.format('RACMO2.3p2','ANT27_SMB_cumul','','',LMAX,RACMO_START,RACMO_END)
            RACMO_file = ff.format('RACMO2.3p2','ANT27_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            RACMO = np.loadtxt(base_dir.joinpath('RACMO','XANT27_1979-2022',subdir,RACMO_file))[:nmon,:]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(RACMO[:,0] == m)
                if valid:
                    mm, = np.nonzero(RACMO[:,0] == m)
                    tnan[d] = RACMO[mm,1]
                    mnan[d] = RACMO[mm,2] - np.mean(RACMO[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[2], label='RACMO2.3p2', zorder=5)

        # plot GSFC-fdm SMB time series for comparison
        if ('GSFC-fdm' in SMB):
            # read MERRA-2 hybrid surface mass balance file
            subdir = 'HEX_GSFC_FDM_{0}_{1}_{2}_L{3:d}'.format('v1_2_1','ais','SMB_a',LMAX)
            MERRA2_file = ff.format('GSFC_FDM','v1_2_1_SMB_a_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            MERRA2 = np.loadtxt(base_dir.joinpath('MERRA2_hybrid','v1.2.1',subdir,MERRA2_file))
            # plot all dates
            ax.plot(MERRA2[:,1], MERRA2[:,2] - np.mean(MERRA2[:,2]),
                color=plot_colors[3], label='GSFC-fdm v1.2.1', zorder=5)

        # plot JPL-GEMB SMB time series for comparison
        if ('GEMB' in SMB):
            # read GEMB surface mass balance file
            subdir = sd.format('GEMB','v1_2_Antarctica_SMB_cumul','','',LMAX,GEMB_START,GEMB_END)
            GEMB_file = ff.format('GEMB','v1_2_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            GEMB = np.loadtxt(base_dir.joinpath('GEMB','v1.2',subdir,GEMB_file))
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(GEMB[:,0] == m)
                if valid:
                    mm, = np.nonzero(GEMB[:,0] == m)
                    tnan[d] = GEMB[mm,1]
                    mnan[d] = GEMB[mm,2] - np.mean(GEMB[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[4], label='GEMB v1.2', zorder=5)

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax.axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax.axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax.axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[reg],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(reg_text[reg],
            prop=textprops, pad=0, frameon=False, loc=1)
        ax.add_artist(at)
        # set ticks
        major_ticks = np.arange(2002, 2024, 4)
        ax.xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        axlim = ax.set_xlim(2002, 2023.0)
        axlim = ax.set_ylim(ylimits[reg][0], ylimits[reg][1])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax.xaxis.get_major_formatter().set_useOffset(False)
        # add x and y labels
        ax.set_xlabel('Time [Yr]')
        ax.set_ylabel('Mass [Gt]', labelpad=ypad[reg])

    # add legend
    lgd = ax1['AIS'].legend(loc=3,frameon=False)
    # lgd = ax1['AIS'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
    # set width, color and style of lines
    # lgd.get_frame().set_boxstyle('square,pad=0.1')
    # lgd.get_frame().set_edgecolor('black')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
       line.set_linewidth(6)
    #     line.set_linewidth(0)
    for i,text in enumerate(lgd.get_texts()):
        text.set_color(plot_colors[i])
        text.set_weight('bold')

    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.0625,right=0.99,bottom=0.05,top=0.99,wspace=0.2,hspace=0.125)
    figurefile = filepath.joinpath('fig4af_{0}_{1}.pdf'.format(PROC,DREL))
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

    # create figure axis for Antarctic plots
    ax1 = {}
    ax2 = {}
    fig, (ax1['EAIS'], ax1['INTERIOR']) = plt.subplots(num=1,ncols=2,figsize=(8,3.75))
    ylimits = {'EAIS':[-1700,1600,250],'INTERIOR':[-650,650,100]}
    fig_text = {'EAIS':'a)','INTERIOR':'b)'}
    reg_text = {'EAIS':'EAIS','INTERIOR':'EAIS Interior'}
    h = 'S'
    SLF = '_SLF3'
    ypad = dict(EAIS=4,INTERIOR=4)
    for reg,ax in ax1.items():
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'','',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str)
        # ATM_file = ff.format('ATM_Differences',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        ATM_file = ff.format('ATM_Differences',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        ATM_input = np.loadtxt(mascon_dir.joinpath(subdir,ATM_file))[:nmon,:]
        # read GLDAS terrestrial water RMS file
        subdir = sd.format('GLDAS','TWC_V2.1_RMS','','',LMAX,GLDAS_START,GLDAS_END)
        TWC_file = ff.format('GLDAS_TWC_RMS',reg,'','RAD1.5_','',LMAX,gw_str,ds_str)
        TWC_input = np.loadtxt(base_dir.joinpath('GLDAS',subdir,TWC_file))[:nmon,:]
        isvalid, = np.nonzero(np.isfinite(TWC_input[:,2]) &
            (TWC_input[:,0] >= START_MON) & (TWC_input[:,0] <= END_MON))
        TWC_RMS = np.sqrt(np.sum(TWC_input[isvalid,2]**2)/len(isvalid))
        # input estimated SLF monte carlo variance file and calculate RMS
        subdir = sd.format(PROC,DREL,'','_MC',LMAX,START_MON,END_MON)
        SLF_file = ff.format('MC',reg,'','RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
        SLF_input = np.loadtxt(mascon_dir.joinpath(subdir,SLF_file))
        SLF_RMS = np.sqrt(np.sum(SLF_input[:,1]**2)/len(SLF_input))
        # calculate mean of RMS over multiple reanalyses
        OBP_RMS = 0.0
        ATM_RMS = 0.0
        # for j in range(2):
        #     ivalid, = np.nonzero(np.isfinite(OBP_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(OBP_input[:,j+3]))
        #     OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(OBP_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(OBP_input[:,2]))
        OBP_RMS += np.sqrt(np.sum(OBP_input[ivalid,2]**2)/valid_count)
        # for j in range(4):
        #     ivalid, = np.nonzero(np.isfinite(ATM_input[:,j+3]))
        #     valid_count = np.count_nonzero(np.isfinite(ATM_input[:,j+3]))
        #     ATM_RMS += np.sqrt(np.sum(OBP_input[ATM_input,j+3]**2)/valid_count)
        ivalid, = np.nonzero(np.isfinite(ATM_input[:,2]))
        valid_count = np.count_nonzero(np.isfinite(ATM_input[:,2]))
        ATM_RMS += np.sqrt(np.sum(ATM_input[ivalid,2]**2)/valid_count)
        # # divide by the number of reanalyses
        # OBP_RMS /= 2.0
        # ATM_RMS /= 4.0
        # GIA modeling group for region
        g = GIA[h][0]
        # number of rheologies to iterate
        nRheology = len(gia_files[g])
        # iterate through solutions
        mon = np.zeros((nmon),dtype=np.int64)
        tdec = np.zeros((nmon),dtype=np.float64)
        mass = np.zeros((nmon,nRheology),dtype=np.float64)
        satellite_error = np.zeros((nmon),dtype=np.float64)
        for i,F in enumerate(VERSION):
            # subdirectory
            subdir = sd.format(PROC,DREL,F,SLF,LMAX,START_MON,END_MON)
            # read each GIA model
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_Ylms = gravtk.read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
                gia_str = gia_Ylms['title']
                input_file = ff.format(gia_str,reg,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                dinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                mon[:] = dinput[:nmon,0].astype(np.int64)
                tdec[:] = dinput[:nmon,1]
                mass[:,k] = dinput[:nmon,2]
                satellite_error[:] += dinput[:nmon,3]**2
            # calculate mean GIA-corrected mass change (for all Earth rheologies)
            gia_corrected_mean = np.mean(mass, axis=1)
            # GRACE satellite error component, ocean leakage (ECCO-GAD),
            # atmosphere leakage (ATM-GAA) and GLDAS TWC
            grace_error = np.sqrt(np.sum(satellite_error/nRheology + SLF_RMS**2 +
                OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon)
            # calculate variance off of mean for calculating GIA uncertainty
            gia_corrected_variance = np.zeros((nmon))
            gia_corrected_minmax = np.zeros((nmon))
            # calculate GIA uncertainty as "worst-case" not RMS
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_corrected_variance+=np.abs(mass[:,k]-gia_corrected_mean)
            # calculate GIA uncertainty as "worst-case" min max error
            for t in range(nmon):
                gia_corrected_minmax[t]=np.abs(np.max(mass[t,:])-np.min(mass[t,:]))
            # calculate uncertainty in mean GIA
            gia_corrected_error = gia_corrected_variance/(np.float64(nRheology)-1.0)
            gia_corrected_error = gia_corrected_error*np.sign(tdec-tdec.mean())
            gia_corrected_minmax = gia_corrected_minmax*np.sign(tdec-tdec.mean())
            gia_error_rate = (gia_corrected_error[-1]-gia_corrected_error[0])/(tdec[-1]-tdec[0])
            gia_minmax_rate = (gia_corrected_minmax[-1]-gia_corrected_minmax[0])/(tdec[-1]-tdec[0])
            # calculate GIA errors at confidence interval
            # t.ppf parallels tinv in matlab
            tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nRheology-1.0) if g in ('IJ05-R2','SM09') else 2.0
            gia_corrected_conf = tstar*np.abs(gia_error_rate)
            # add to plot with colors and label
            plot_label = plot_title[i]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(mon == m)
                if valid:
                    mm, = np.nonzero(mon == m)
                    tnan[d] = tdec[mm]
                    mnan[d] = gia_corrected_mean[mm] - np.mean(gia_corrected_mean)
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[i], label=plot_label, zorder=3+i)
            # fill between monthly errors
            ax.fill_between(tnan, mnan-grace_error, y2=mnan+grace_error,
               color=plot_colors[i], alpha=0.35, zorder=1+i)

        # plot RACMO SMB time series for comparison
        if ('RACMO' in SMB):
            # read RACMO surface mass balance file
            subdir = sd.format('RACMO2.3p2','ANT27_SMB_cumul','','',LMAX,RACMO_START,RACMO_END)
            RACMO_file = ff.format('RACMO2.3p2','ANT27_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            RACMO = np.loadtxt(base_dir.joinpath('RACMO','XANT27_1979-2022',subdir,RACMO_file))[:nmon,:]
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(RACMO[:,0] == m)
                if valid:
                    mm, = np.nonzero(RACMO[:,0] == m)
                    tnan[d] = RACMO[mm,1]
                    mnan[d] = RACMO[mm,2] - np.mean(RACMO[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[2], label='RACMO2.3p2', zorder=5)

        # plot GSFC-fdm SMB time series for comparison
        if ('GSFC-fdm' in SMB):
            # read MERRA-2 hybrid surface mass balance file
            subdir = 'HEX_GSFC_FDM_{0}_{1}_{2}_L{3:d}'.format('v1_2_1','ais','SMB_a',LMAX)
            MERRA2_file = ff.format('GSFC_FDM','v1_2_1_SMB_a_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            MERRA2 = np.loadtxt(base_dir.joinpath('MERRA2_hybrid','v1.2.1',subdir,MERRA2_file))
            # plot all dates
            ax.plot(MERRA2[:,1], MERRA2[:,2] - np.mean(MERRA2[:,2]),
                color=plot_colors[3], label='GSFC-fdm v1.2.1', zorder=5)

        # plot JPL-GEMB SMB time series for comparison
        if ('GEMB' in SMB):
            # read GEMB surface mass balance file
            subdir = sd.format('GEMB','v1_2_Antarctica_SMB_cumul','','',LMAX,GEMB_START,GEMB_END)
            GEMB_file = ff.format('GEMB','v1_2_SMB_cumul_',reg,'RAD1.5_',ocean_str,LMAX,gw_str,ds_str)
            GEMB = np.loadtxt(base_dir.joinpath('GEMB','v1.2',subdir,GEMB_file))
            # create a time series with nans for missing months
            tnan = np.full_like(month,np.nan,dtype=np.float64)
            mnan = np.full_like(month,np.nan,dtype=np.float64)
            for d,m in enumerate(month):
                valid = np.count_nonzero(GEMB[:,0] == m)
                if valid:
                    mm, = np.nonzero(GEMB[:,0] == m)
                    tnan[d] = GEMB[mm,1]
                    mnan[d] = GEMB[mm,2] - np.mean(GEMB[:,2])
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[4], label='GEMB v1.2', zorder=5)

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax.axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax.axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax.axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[reg],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(reg_text[reg],
            prop=textprops, pad=0, frameon=False, loc=1)
        ax.add_artist(at)
        # set ticks
        major_ticks = np.arange(2002, 2024, 4)
        ax.xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        axlim = ax.set_xlim(2002, 2023.0)
        axlim = ax.set_ylim(ylimits[reg][0], ylimits[reg][1])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax.xaxis.get_major_formatter().set_useOffset(False)
        # add x and y labels
        ax.set_xlabel('Time [Yr]')
        ax.set_ylabel('Mass [Gt]', labelpad=ypad[reg])

    # add legend
    lgd = ax1['EAIS'].legend(loc=3,frameon=False)
    # lgd = ax1['EAIS'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
    # set width, color and style of lines
    # lgd.get_frame().set_boxstyle('square,pad=0.1')
    # lgd.get_frame().set_edgecolor('black')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
       line.set_linewidth(6)
    #     line.set_linewidth(0)
    for i,text in enumerate(lgd.get_texts()):
        text.set_color(plot_colors[i])
        text.set_weight('bold')

    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.08,right=0.99,bottom=0.1,top=0.99,wspace=0.2,hspace=0.125)
    figurefile = filepath.joinpath('fig4ab_{0}_{1}.pdf'.format(PROC,DREL))
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser()
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, default=None,
        choices=['CSR','GFZ','JPL'],
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # start and end GRACE/GRACE-FO months
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
    parser.add_argument('--smb','-s',
        type=str, default=[], choices=('RACMO','GSFC-fdm','GEMB'), nargs='+',
        help='Plot Surface Mass Balance (SMB) time series')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program for parameters
    plot_mascon_SLF_timeseries(args.directory,args.center,args.release,
        args.start,args.end,args.missing,SMB=args.smb)

# run main program
if __name__ == '__main__':
    main()
