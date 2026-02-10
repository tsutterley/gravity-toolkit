#!/usr/bin/env python
u"""
plot_mascon_SLF_combined.py
Written by Tyler Sutterley (05/2023)

CALLING SEQUENCE:
    python plot_mascon_SLF_combined.py --start 4 --end 230

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: place imports behind try/except statements
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
import gravity_toolkit.time
import gravity_toolkit.units
import gravity_toolkit.spatial
import gravity_toolkit.utilities
from gravity_toolkit.gen_harmonics import gen_harmonics
from gravity_toolkit.read_GIA_model import read_GIA_model

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

# PURPOSE: plot and calculate effects of SLF on mascon time series
def plot_mascon_SLF_combined(base_dir,PROC,DREL,START_MON,END_MON,MISSING):
    # directory setup
    mascon_dir = base_dir.joinpath('GRACE','mascons')
    # GIA and REGION
    GIA = {}
    gia_files = {}
    gia_mean_str = {}
    REGION = {}

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
    # mean title for each GIA string
    gia_mean_str['IJ05-R2'] = 'Ivins_R2_mean'
    gia_mean_str['W12a'] = 'W12a_mean'
    gia_mean_str['ICE6G'] = 'ICE6G_VM5_mean'
    gia_mean_str['ICE6G-D'] = 'ICE6G-D_mean'
    gia_mean_str['AW13-ICE6G'] = 'AW13_ice6g_mean'
    gia_mean_str['ascii'] = '_AW13_IJ05_ICE6G_mean'
    gia_mean_str['SM09'] = 'SM09_HUY2_mean'
    gia_mean_str['Caron'] = 'Caron_expt'
    # regions
    REGION['N'] = 'GIS'
    REGION['S'] = 'AIS'
    REMOVE = dict(S='AIS',N='ARC')

    # Land-Sea Mask (with Antarctica from Rignot et al., 2017)
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    LSMASK = gravity_toolkit.utilities.get_data_path(['data','landsea_hd.nc'])
    landsea = gravity_toolkit.spatial().from_netCDF4(LSMASK,
        date=False, varname='LSMASK')
    # create land function
    nth,nphi = landsea.shape
    land_function = np.zeros((nth,nphi),dtype=np.float64)
    # combine land and island levels for land function
    indx,indy = np.nonzero((landsea.data >= 1) & (landsea.data <= 3))
    land_function[indx,indy] = 1.0
    # calculate ocean function from land function
    ocean_function = 1.0 - land_function
    # convert ocean function into a series of spherical harmonics
    # (Note that LMAX=0 in the call to gen_harmonics)
    Ylms = gen_harmonics(ocean_function.T, landsea.lon, landsea.lat,
        LMAX=0, PLM=np.ones((1,1,nth)))
    # total area of ocean calculated by integrating the ocean function
    # Average Radius of the Earth [mm]
    rad_e = 10.0*gravity_toolkit.units().rad_e
    ocean_area = 4.0*np.pi*(rad_e**2)*Ylms.clm[0,0]

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

    # input leakage file
    leakage_file = mascon_dir.joinpath('HEX_TMB_LEAKAGE_SPH_CAP_MSCNS_L60',
        'HEX_TMB_LEAKAGE_HOLE_SPH_CAP_RAD1.5_L60_r250km_OCN.txt')
    regional_leakage = {}
    with open(leakage_file,'r') as f:
        file_contents = f.read().splitlines()
        for line in file_contents:
            line_contents = line.split()
            regional_leakage[line_contents[0]] = np.float64(line_contents[1])

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

    # subdirectory and input file formats
    sd = 'HEX_{0}_{1}_SPH_CAP_MSCNS{2}_L{3:d}_{4:03d}-{5:03d}'
    ff = '{0}_{1}{2}_SPH_CAP_{3}{4}L{5:d}{6}{7}.txt'

    # create figure axis
    ax1 = {}
    ax2 = {}
    ax3 = {}
    fig, ((ax1['GIS'],ax3['GIS']),(ax1['WAIS'],ax3['WAIS']),
        (ax1['EAIS'],ax3['EAIS']),(ax1['APIS'],ax3['APIS'])) = \
        plt.subplots(num=1,nrows=4,ncols=2,sharex=True,figsize=(9,10))

    # create a set of months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    # figure labels
    FLAG = ['HEX','HEX']
    plot_colors = ['darkorchid','mediumseagreen']
    plot_title = ['No Sea Level Correction','Sea Level Fingerprint']
    iter_label = ['No SL','SLF']
    fig_text = dict(GIS='a)', WAIS='b)', EAIS='c)', APIS='d)')
    ylimits = dict(GIS=[-5500,1000],WAIS=[-3000,600],EAIS=[-400,1500],APIS=[-550,200])
    hem = dict(GIS='N',WAIS='S',EAIS='S',APIS='S')
    ypad = dict(GIS=4,WAIS=4,EAIS=9.5,APIS=9.5)
    ysea = dict(GIS=[14,-2,-2],WAIS=[8,-1,-1],EAIS=[1,-4,-1],APIS=[1,0,-1])
    seapad = dict(GIS=4,WAIS=4,EAIS=4,APIS=9.5)
    for reg,ax in ax1.items():
        # hemisphere flag for region
        h = hem[reg]
        # set ice sheet only for AIS and GIS
        SLF = ['','_SLF3','_SLF6'] if (reg == 'GIS') else ['','_SLF3','_SLF5']
        # leakage fraction for regions
        leakage_fraction = regional_leakage[reg]
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        # ATM_file = ff.format('ATM_Differences',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        ATM_file = ff.format('ATM_Differences',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        ATM_input = np.loadtxt(mascon_dir.joinpath(subdir,ATM_file))[:nmon,:]
        # read GLDAS terrestrial water RMS file
        subdir = sd.format('GLDAS','TWC_V2.1_RMS','',LMAX,GLDAS_START,GLDAS_END)
        TWC_file = ff.format('GLDAS_TWC_RMS',reg,'','RAD1.5_','',LMAX,gw_str,ds_str)
        TWC_input = np.loadtxt(base_dir.joinpath('GLDAS',subdir,TWC_file))[:nmon,:]
        isvalid, = np.nonzero(np.isfinite(TWC_input[:,2]) &
            (TWC_input[:,0] >= START_MON) & (TWC_input[:,0] <= END_MON))
        TWC_RMS = np.sqrt(np.sum(TWC_input[isvalid,2]**2)/len(isvalid))
        # input estimated SLF monte carlo variance file and calculate RMS
        subdir = sd.format(PROC,DREL,'_MC',LMAX,START_MON,END_MON)
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
        for i,F in enumerate(FLAG):
            # subdirectory
            subdir = sd.format(PROC,DREL,SLF[i],LMAX,START_MON,END_MON)
            # read each GIA model
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_Ylms = read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
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
                    mnan[d] = gia_corrected_mean[mm]-gia_corrected_mean[0]
            # plot all dates
            ax.plot(tnan, mnan, color=plot_colors[i], label=plot_label, zorder=2)
            # fill between monthly errors
            ax.fill_between(tnan, mnan-grace_error, y2=mnan+grace_error,
                color=plot_colors[i], alpha=0.5, zorder=1)
            # converting gigatonnes to milligrams then to mm sea level
            mm_sealevel = -1e18*(gia_corrected_mean-gia_corrected_mean[0])/ocean_area
            mm_sealevel_error = 1e18*np.sqrt(grace_error**2 +
                (tstar*gia_corrected_error[-1])**2)/ocean_area

        # vertical line denoting the accelerometer shutoff
        acc = gravity_toolkit.time.convert_calendar_decimal(2016,9,
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
        lab = ax1[reg].set_ylabel('Mass [Gt]', labelpad=ypad[reg])
        # set ticks
        major_ticks = np.arange(2002, 2024, 4)
        ax.xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        axlim = ax.set_xlim(2002, 2022.25)
        axlim = ax.set_ylim(ylimits[reg])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker
        ax.xaxis.get_major_formatter().set_useOffset(False)
        # y labels on both sides of plot
        ax2[reg] = ax.twinx()
        # add plot to hidden mm sea level axis
        ax2[reg].plot(tdec, mm_sealevel, visible=False)
        seaticks = np.arange(ysea[reg][0],ysea[reg][1]+ysea[reg][2],ysea[reg][2])
        ax2[reg].set_yticks(seaticks)
        ax2[reg].set_ylim(-1e18*ylimits[reg][0]/ocean_area,-1e18*ylimits[reg][1]/ocean_area)
        ax2[reg].tick_params(axis='y',colors='black',which='both',direction='in')
        # formatted ticks on sea level axis
        ax2[reg].set_ylabel('Equivalent Sea Level [mm]',labelpad=seapad[reg],color='black')
        for tl in ax2[reg].get_yticklabels():
            tl.set_color('black')

    # add x label
    ax1['APIS'].set_xlabel('Time [Yr]')
    # add legend
    lgd = ax1['APIS'].legend(loc=3,frameon=False)
    # lgd = ax1['APIS'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
    # set width, color and style of lines
    # lgd.get_frame().set_boxstyle('square,pad=0.1')
    # lgd.get_frame().set_edgecolor('black')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    #     line.set_linewidth(0)
    # for i,text in enumerate(lgd.get_texts()):
    #     text.set_color(plot_colors[i])
    #     text.set_weight('bold')


    # flags for creating plots
    FLAG = ['HEX','HEX']
    SLF = ['','_SLF3']
    plot_colors = ['black','mediumseagreen']
    plot_titles = ['No SL Correction','Iteration 3']
    fig_text = dict(GIS='e)', WAIS='f)', EAIS='g)', APIS='h)')
    # y limits for plot
    ylimits = [-90,170,40]
    for reg,ax in ax3.items():
        # hemisphere flag for region
        h = hem[reg]
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        # ATM_file = ff.format('ATM_Differences',reg,'_3D','',ocean_str,LMAX,gw_str,ds_str)
        ATM_file = ff.format('ATM_Differences',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        ATM_input = np.loadtxt(mascon_dir.joinpath(subdir,ATM_file))[:nmon,:]
        # read GLDAS terrestrial water RMS file
        subdir = sd.format('GLDAS','TWC_V2.1_RMS','',LMAX,GLDAS_START,GLDAS_END)
        TWC_file = ff.format('GLDAS_TWC_RMS',reg,'','RAD1.5_','',LMAX,gw_str,ds_str)
        TWC_input = np.loadtxt(base_dir.joinpath('GLDAS',subdir,TWC_file))[:nmon,:]
        isvalid, = np.nonzero(np.isfinite(TWC_input[:,2]) &
            (TWC_input[:,0] >= START_MON) & (TWC_input[:,0] <= END_MON))
        TWC_RMS = np.sqrt(np.sum(TWC_input[isvalid,2]**2)/len(isvalid))
        # input estimated SLF monte carlo variance file and calculate RMS
        subdir = sd.format(PROC,DREL,'_MC',LMAX,START_MON,END_MON)
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
        mass_m1 = np.zeros((nmon),dtype=np.float64)
        satellite_error = np.zeros((nmon),dtype=np.float64)
        for i,F in enumerate(FLAG):
            # subdirectory
            subdir = sd.format(PROC,DREL,SLF[i],LMAX,START_MON,END_MON)
            # read each GIA model
            for k,GIA_FILE in enumerate(gia_files[g]):
                gia_Ylms = read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
                gia_str = gia_Ylms['title']
                input_file = ff.format(gia_str,reg,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                dinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                mon[:] = dinput[:nmon,0]
                tdec[:] = dinput[:nmon,1]
                mass[:,k] = dinput[:nmon,2]
                satellite_error[:] += dinput[:nmon,3]**2
            # calculate mean GIA-corrected mass change (for all Earth rheologies)
            gia_corrected_mean = np.mean(mass, axis=1)
            # GRACE satellite error component, ocean leakage (ECCO-GAD),
            # atmosphere leakage (ATM-GAA) and GLDAS TWC
            grace_error = np.sqrt(np.sum(satellite_error/nRheology + SLF_RMS**2 +
                OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon)
            if (i > 0):
                # add to plot with colors and label
                plot_label = u'{0} \u2013 {1}'.format(plot_titles[i],plot_titles[i-1])
                residual = gia_corrected_mean - gia_corrected_mean[0] - mass_m1
                # create a time series with nans for missing months
                tnan = np.full_like(month,np.nan,dtype=np.float64)
                rnan = np.full_like(month,np.nan,dtype=np.float64)
                for d,m in enumerate(month):
                    valid = np.count_nonzero(mon == m)
                    if valid:
                        mm, = np.nonzero(mon == m)
                        tnan[d] = tdec[mm]
                        rnan[d] = residual[mm]
                # plot all dates
                ax.plot(tnan, rnan, color=plot_colors[i], label=plot_label)
            # save mass for iteration
            mass_m1 = gia_corrected_mean - gia_corrected_mean[0]
            error_m1 = np.copy(grace_error)

        # vertical line denoting the accelerometer shutoff
        acc = gravity_toolkit.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax.axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax.axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax.axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax.axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add horizontal line at 0
        ax.axhline(0.0, color='black', ls='dashed', dashes=(11,5), lw=0.5)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[reg],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax.add_artist(at)
        ax.set_ylabel('Mass Difference [Gt]')
        # set ticks
        major_ticks = np.arange(2002, 2024, 4)
        ax.xaxis.set_ticks(major_ticks)
        minor_ticks = sorted(set(np.arange(2002, 2022, 1)) - set(major_ticks))
        ax.xaxis.set_ticks(minor_ticks, minor=True)
        axlim = ax.set_xlim(2002, 2022.25)
        data_ticks = np.arange(ylimits[1]-10,ylimits[0]-ylimits[2],-ylimits[2])
        ax.yaxis.set_ticks(data_ticks[::-1])
        axlim = ax.set_ylim(ylimits[:2])
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker
        ax.xaxis.get_major_formatter().set_useOffset(False)

    # add x labels
    ax3['APIS'].set_xlabel('Time [Yr]')
    # add legend
    lgd = ax3['APIS'].legend(loc=3,frameon=False)
    # lgd = ax3['APIS'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
    # set width, color and style of lines
    # lgd.get_frame().set_boxstyle('square,pad=0.1')
    # lgd.get_frame().set_edgecolor('black')
    lgd.get_frame().set_alpha(1.0)
    for line in lgd.get_lines():
        line.set_linewidth(6)
    #     line.set_linewidth(0)
    # for i,text in enumerate(lgd.get_texts()):
    #     text.set_color(plot_colors[i+1])
    #     text.set_weight('bold')

    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.08,right=0.99,bottom=0.04,top=0.99,
        hspace=0.08,wspace=0.28)
    figurefile = filepath.joinpath('fig5ah_{0}_{1}.pdf'.format(PROC,DREL))
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
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program for parameters
    plot_mascon_SLF_combined(args.directory,args.center,args.release,
        args.start,args.end,args.missing)

# run main program
if __name__ == '__main__':
    main()
