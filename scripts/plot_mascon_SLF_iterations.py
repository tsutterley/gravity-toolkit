#!/usr/bin/env python
u"""
plot_mascon_SLF_iterations.py
Written by Tyler Sutterley (05/2023)

CALLING SEQUENCE:
    python plot_mascon_SLF_iterations.py --start 4 --end 230

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Releases
    -S X, --start X: starting GRACE month for time series
    -E X, --end X: ending GRACE month for time series
    -M X, --missing X: Missing GRACE months in time series
    -s, --scenarios: Run with all sea level scenarios

UPDATE HISTORY:
    Updated 05/2023: split S2 tidal aliasing terms into GRACE and GRACE-FO eras
        edit monte carlo regression function to use orders and cycles keywords
        use fit module for getting tidal aliasing terms
        use pathlib to define and operate on paths
    Updated 03/2023: place imports behind try/except statements
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: simplify to not run piecewise regression schemes
    Updated 10/2022: simplify to only Greenland and Antarctic plots
    Updated 07/2022: set plot tick formatter to not use offsets
    Updated 05/2022: use argparse descriptions within documentation
    Updated 04/2022: added AW13 combination IJ05-R2 and ICE6G models
    Updated 05/2021: new TMB leakage directory path. option for SLF scenarios
        define int/float precision to prevent deprecation warning
    Updated 04/2021: using argparse to set parameters
    Updated 01/2021: updated for new open-source processing scheme
    Updated 11/2019: adjust zorder of plot lines and error shading
    Updated 10/2019: changing Y/N flags to True/False.  getopt to set parameters
        Include Caron et al. (2018) and A et al. ICE-6G model outputs
    Updated 08/2019: include GRACE-FO months. empty span over GRACE/GRACE-FO gap
        updated filenames to use a more common format
    Updated 01/2019: include 161-day S2 tidal aliasing in regressions
    Updated 05/2018: output estimated trend and error for each leakage component
        include estimated SLF error from monte carlo simulation of GRACE errors
        include SLF leakage effects for single ice sheet (AIS/GIS) configuration
        print the total equivalent sea level contribution for each region
    Updated 03/2018: estimate atmospheric dealiasing errors from reanalysis
        include estimated total water storage error.  trends for glacier regions
    Updated 01/2018: estimate oceanic dealiasing errors from JPL ECCO2 outputs
    Updated 10/2017: output trends to file.  calculate total sea level change
    Updated 09/2017: calculate with longer time series and method improvements
    Updated 08/2017: complete rewrite of program to create figures for paper
        iterate through different estimates from different GIA models
        include estimates using ICE-6G GIA. calculate ice sheet sea level change
        estimate sea level contribution on righthand plot axis of time series
        create a plot showing the effect of the iterations for supplement
    Written 08/2017
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

# PURPOSE: plot and calculate effects of SLF on mascon time series
def plot_mascon_SLF_iterations(base_dir,PROC,DREL,START_MON,END_MON,MISSING,SCENARIOS=False):
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
    gia_mean_str['ICE6G-D'] = 'ICE6G-D_High_Res'
    gia_mean_str['AW13-ICE6G'] = 'AW13_ice6g_mean'
    gia_mean_str['ascii'] = 'AW13_IJ05_ICE6G_mean'
    gia_mean_str['SM09'] = 'SM09_HUY2_mean'
    gia_mean_str['Caron'] = 'Caron_expt'
    # regions
    REGION['N'] = 'GIS'
    REGION['S'] = 'AIS'
    REMOVE = dict(S='AIS',N='ARC')

    # figure labels
    fig_text = dict(N='b)', S='a)', GIS='a)', WAIS='b)', EAIS='c)', APIS='d)')

    # Land-Sea Mask (with Antarctica from Rignot et al., 2017)
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF file for reading
    LSMASK = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    landsea = gravtk.spatial().from_netCDF4(LSMASK,
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
    Ylms = gravtk.gen_harmonics(ocean_function.T, landsea.lon, landsea.lat,
        LMAX=0, PLM=np.ones((1,1,nth)))
    # total area of ocean calculated by integrating the ocean function
    # Average Radius of the Earth [mm]
    rad_e = 10.0*gravtk.units().rad_e
    ocean_area = 4.0*np.pi*(rad_e**2)*Ylms.clm[0,0]

    # Setting output error alpha with confidence interval
    CONF = 0.95
    alpha = 1.0 - CONF
    # monte carlo regression runs
    RUNS = 20000

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

    # flags for calculating trends to populate tables
    FLAG = ['HEX','HEX']
    SLF = ['','_SLF3']
    iter_label = ['No SL','SLF']
    # if running with all scenarios
    if SCENARIOS:
        FLAG.extend(['HEX','HEX','HEX'])
        SLF.extend(['_SLF4','_SLF5','_SLF6'])
        iter_label.extend(['ISSLF','AISSLF','GISSLF'])
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

    # run for multiple regions to get trends for the table
    regions = []
    AIS_regions = ['AIS','WAIS','EAIS','APIS','INTERIOR','QML','IIpp','CpDc','DDpi','FpG','GH2','GH3']
    GIS_regions = ['GIS','NW','NN','NE','SW','SE']
    GIC_regions = ['CBI','CDE','ICL','SVB','FJL','SZEM','NZEM','ALK','DEN','PNW','PAT']
    HEM = ['S']*len(AIS_regions) + ['N']*len(GIS_regions) + ['N']*len(GIC_regions)
    remove = ['AIS']*len(AIS_regions) + ['ARC']*len(GIS_regions) + ['ARC']*len(GIC_regions)
    regions.extend(AIS_regions)
    regions.extend(GIS_regions)
    regions.extend(GIC_regions)
    # file with output table
    tablefile = filepath.joinpath('table_1_{0}_{1}_obp_atm_twc.txt'.format(PROC,DREL))
    fid = tablefile.open(mode='w', encoding='utf8')
    # fid = sys.stdout
    for h,reg,rem in zip(HEM,regions,remove):
        # leakage fraction for regions
        if reg in regional_leakage.keys():
            leakage_fraction = regional_leakage[reg]
        else:
            leakage_fraction = 0.0
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
        # read GIA uncertainty from Caron et al. (2018)
        subdir = sd.format(PROC,DREL,'_SLF3',LMAX,START_MON,END_MON)
        input_file = ff.format('Caron_Error',reg,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
        CARON_RMS = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
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
        # iterate through GIA solutions
        for j,g in enumerate(GIA[h]):
            print('{0} {1}'.format(reg, gia_mean_str[g]), file=fid)
            mon = np.zeros((nmon),dtype=np.int64)
            tdec = np.zeros((nmon),dtype=np.float64)
            # number of rheologies to iterate
            nRheology = len(gia_files[g])
            mass = np.zeros((nmon,nRheology),dtype=np.float64)
            complement = np.zeros((nmon),dtype=np.float64)
            satellite_error = np.zeros((nmon),dtype=np.float64)
            bx1 = np.zeros_like(FLAG,dtype=np.float64)
            ex1 = np.zeros_like(FLAG,dtype=np.float64)
            bx2 = np.zeros_like(FLAG,dtype=np.float64)
            ex2 = np.zeros_like(FLAG,dtype=np.float64)
            # iterate through solutions
            for i,F in enumerate(FLAG):
                # subdirectory
                subdir = sd.format(PROC,DREL,SLF[i],LMAX,START_MON,END_MON)
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
                    # read remove file to calculate complement
                    input_file = ff.format(gia_str,rem,'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                    rinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                    complement += (rinput[:nmon,2] - mass[:,k])
                # calculate mean GIA-corrected mass change (for all Earth rheologies)
                gia_corrected_mean = np.mean(mass, axis=1)
                # GRACE satellite error component with SLF uncertainty
                grace_error = np.sqrt(satellite_error/np.float64(nRheology))
                # calculate variance off of mean for calculating GIA uncertainty
                gia_corrected_variance = np.zeros((nmon))
                gia_corrected_minmax = np.zeros((nmon))
                # calculate complement time series of region
                complement = complement/np.float64(nRheology)
                # calculate GIA uncertainty as "worst-case" not RMS
                for k,GIA_FILE in enumerate(gia_files[g]):
                    gia_corrected_variance+=np.abs(mass[:,k]-gia_corrected_mean)
                # calculate GIA uncertainty as "worst-case" min max error
                for t in range(nmon):
                    gia_corrected_minmax[t]=np.abs(np.max(mass[t,:])-np.min(mass[t,:]))
                # if using Caron et al. (2018) use RMS of covariance errors
                if (g == 'Caron'):
                    # calculate uncertainty in mass drift
                    gia_corrected_error = CARON_RMS*(tdec-tdec.mean())
                    gia_corrected_minmax = CARON_RMS*(tdec-tdec.mean())
                    gia_error_rate = np.copy(CARON_RMS)
                    gia_minmax_rate = np.copy(CARON_RMS)
                    tstar = 1.0
                else:
                    # calculate uncertainty in mean GIA
                    gia_corrected_error = gia_corrected_variance/(np.float64(nRheology)-1.0)
                    gia_corrected_error = gia_corrected_error*np.sign(tdec-tdec.mean())
                    gia_corrected_minmax = gia_corrected_minmax*np.sign(tdec-tdec.mean())/2.0
                    gia_error_rate = (gia_corrected_error[-1]-gia_corrected_error[0])/(tdec[-1]-tdec[0])
                    gia_minmax_rate = (gia_corrected_minmax[-1]-gia_corrected_minmax[0])/(tdec[-1]-tdec[0])
                    # calculate GIA errors at confidence interval
                    # t.ppf parallels tinv in matlab
                    tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nRheology-1.0) if g in ('IJ05-R2','SM09') else 2.0
                # GIA errors at confidence interval
                gia_corrected_conf = tstar*np.abs(gia_error_rate)
                # build terms for S2 tidal aliasing
                TERMS = gravtk.time_series.aliasing_terms(tdec)
                # fit a linear model for trend and acceleration
                bfit = gravtk.time_series.regress(tdec, gia_corrected_mean,
                    ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
                ofit = gravtk.time_series.regress(tdec, gia_corrected_mean,
                    DATA_ERR=OBP_RMS, ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
                afit = gravtk.time_series.regress(tdec, gia_corrected_mean,
                    DATA_ERR=ATM_RMS, ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
                tfit = gravtk.time_series.regress(tdec, gia_corrected_mean,
                    DATA_ERR=TWC_RMS, ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
                cfit = gravtk.time_series.regress(tdec, complement,
                    ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
                # run a monte carlo regression for trend and acceleration
                combined_error = np.sqrt(grace_error**2 + SLF_RMS**2)
                mc = monte_carlo_regress(tdec, gia_corrected_mean, combined_error,
                    tstar*(gia_corrected_error-gia_corrected_error[0]),
                    leakage_fraction*complement, OBP_RMS, ATM_RMS, TWC_RMS,
                    ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS, RUNS=RUNS,
                    RMS=False, CONF=CONF)
                # print trend (x1) coefficients to file
                bx1[i],cx1 = bfit['beta'][1],cfit['beta'][1]
                ex1[i] = bfit['error'][1] + gia_corrected_conf + \
                    np.abs(leakage_fraction*cx1) + np.abs(ofit['error'][1]) + \
                    np.abs(afit['error'][1]) + np.abs(tfit['error'][1])
                # print acceleration (x2) coefficients to file
                bx2[i],cx2 = 2.0*bfit['beta'][2],2.0*cfit['beta'][2]
                ex2[i] = 2.0*bfit['error'][2] + np.abs(leakage_fraction*cx2) + \
                    np.abs(2.0*ofit['error'][2]) + np.abs(2.0*afit['error'][2]) + \
                    np.abs(2.0*tfit['error'][2])
                args = (iter_label[i],bx1[i],ex1[i],bx2[i],ex2[i])
                print('{0}\tx1={1:f}+/-{2:f}\tx2={3:f}+/-{4:f}'.format(*args),file=fid)
                # print monte carlo trend and acceleration coefficients to file
                a=(mc['beta'][1],mc['error'][1],2.*mc['beta'][2],2.*mc['error'][2])
                print('\tmx1={0:f}+/-{1:f}\tmx2={2:f}+/-{3:f}'.format(*a),file=fid)
                # converting gigatonnes to milligrams then to mm sea level
                mm_sealevel = -1e18*(gia_corrected_mean-gia_corrected_mean[0])/ocean_area
                cumulative_gia = tstar*(gia_corrected_error[-1]-gia_corrected_error[0])
                mm_sealevel_error = 1e18*np.sqrt(np.sum(satellite_error/nRheology +
                    SLF_RMS**2 + OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon +
                    cumulative_gia**2)/ocean_area
                args = (mm_sealevel[-1],mm_sealevel_error)
                print('\ttotal sea level: {0:f}+/-{1:f} mm'.format(*args),file=fid)
            # if running with all scenarios
            if SCENARIOS:
                # print percent differences
                dx1 = 100.0*np.abs((bx1[1]-bx1[0])/bx1[0])
                dx2 = 100.0*np.abs((bx2[1]-bx2[0])/bx2[0])
                print('(2-1)%\tx1={0:f}\t\t\tx2:{1:f}'.format(dx1,dx2),file=fid)
                dx1 = 100.0*np.abs((bx1[2]-bx1[1])/bx1[1])
                dx2 = 100.0*np.abs((bx2[2]-bx2[1])/bx2[1])
                print('(3-2)%\tx1={0:f}\t\t\tx2:{1:f}'.format(dx1,dx2),file=fid)
            # calculate uncertainties associated with geophysical corrections
            emc = monte_carlo_regress(tdec, np.ones((nmon)), grace_error,
                0.0, 0.0, 0.0, 0.0, 0.0, ORDER=2, CYCLES=[0.5,1.0],
                TERMS=TERMS, RUNS=RUNS, RMS=False, CONF=CONF)
            smc = monte_carlo_regress(tdec, np.ones((nmon)), SLF_RMS,
                0.0, 0.0, 0.0, 0.0, 0.0, ORDER=2, CYCLES=[0.5,1.0],
                TERMS=TERMS, RUNS=RUNS, RMS=False, CONF=CONF)
            lmc = monte_carlo_regress(tdec, np.ones((nmon)), 0.0, 0.0,
                leakage_fraction*complement, 0.0, 0.0, 0.0,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS, RUNS=RUNS,
                RMS=False, CONF=CONF)
            omc = monte_carlo_regress(tdec, np.ones((nmon)), 0.0, 0.0, 0.0,
                OBP_RMS, 0.0, 0.0, ORDER=2, CYCLES=[0.5,1.0],
                TERMS=TERMS, RUNS=RUNS, RMS=False, CONF=CONF)
            amc = monte_carlo_regress(tdec, np.ones((nmon)), 0.0, 0.0, 0.0,
                0.0, ATM_RMS, 0.0, ORDER=2, CYCLES=[0.5,1.0],
                TERMS=TERMS, RUNS=RUNS, RMS=False, CONF=CONF)
            tmc = monte_carlo_regress(tdec, np.ones((nmon)), 0.0, 0.0, 0.0,
                0.0, 0.0, TWC_RMS, ORDER=2, CYCLES=[0.5,1.0],
                TERMS=TERMS, RUNS=RUNS, RMS=False, CONF=CONF)
            # print uncertainties associated with geophysical corrections
            args = (emc['error'][1],2.0*emc['error'][2])
            print('grace={0:f}\tacc={1:f}'.format(*args),file=fid)
            args = (smc['error'][1],2.0*smc['error'][2])
            print('slf={0:f}\tacc={1:f}'.format(*args),file=fid)
            args = (omc['error'][1],2.0*omc['error'][2])
            print('obp={0:f}\tacc={1:f}'.format(*args),file=fid)
            args = (amc['error'][1],2.0*amc['error'][2])
            print('atm={0:f}\tacc={1:f}'.format(*args),file=fid)
            args = (tmc['error'][1],2.0*tmc['error'][2])
            print('tws={0:f}\tacc={1:f}'.format(*args),file=fid)
            args = (lmc['error'][1],2.0*lmc['error'][2])
            print('leak={0:f}\tacc={1:f}'.format(*args),file=fid)
            print('gia={0:f}\n'.format(gia_corrected_conf),file=fid)

    # flags for creating plots
    FLAG = ['HEX','HEX']
    mm_total = np.zeros((2))
    mm_error = np.zeros((2))
    plot_colors = ['darkorchid','mediumseagreen']
    plot_title = ['No Sea Level Correction','Sea Level Fingerprint']
    iter_label = ['No SL','SLF']
    # create a set of months
    month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))

    # create figure axis
    ax = {}
    ax1 = {}
    fig, (ax['S'],ax['N']) = plt.subplots(num=1,ncols=2,figsize=(8.5,4))
    # y limits for plot
    ylimits = [-5200,1200,500]
    for h,reg in REGION.items():
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
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str)
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
        complement = np.zeros((nmon),dtype=np.float64)
        bx1 = np.zeros_like(FLAG,dtype=np.float64)
        ex1 = np.zeros_like(FLAG,dtype=np.float64)
        bx2 = np.zeros_like(FLAG,dtype=np.float64)
        ex2 = np.zeros_like(FLAG,dtype=np.float64)
        for i,F in enumerate(FLAG):
            # subdirectory
            subdir = sd.format(PROC,DREL,SLF[i],LMAX,START_MON,END_MON)
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
                # read remove file to calculate complement
                input_file = ff.format(gia_str,REMOVE[h],'',atm_str,ocean_str,LMAX,gw_str,ds_str)
                rinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
                complement += (rinput[:nmon,2] - mass[:,k])
                # ax[h].plot(tdec,mass[:,k]-mass[0,k])
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
            # tstar = 2.0
            gia_corrected_conf = tstar*np.abs(gia_error_rate)
            # build terms for S2 tidal aliasing
            TERMS = gravtk.time_series.aliasing_terms(tdec)
            # fit a linear model for trend and acceleration
            bfit = gravtk.time_series.regress(tdec, gia_corrected_mean,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
            ofit = gravtk.time_series.regress(tdec, gia_corrected_mean, DATA_ERR=OBP_RMS,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
            afit = gravtk.time_series.regress(tdec, gia_corrected_mean, DATA_ERR=ATM_RMS,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
            tfit = gravtk.time_series.regress(tdec, gia_corrected_mean, DATA_ERR=TWC_RMS,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
            cfit = gravtk.time_series.regress(tdec, complement,
                ORDER=2, CYCLES=[0.5,1.0], TERMS=TERMS)
            bx1[i],cx1 = bfit['beta'][1], cfit['beta'][1]
            ex1[i] = bfit['error'][1] + gia_corrected_conf + \
                np.abs(leakage_fraction*cx1) + np.abs(ofit['error'][1]) + \
                np.abs(afit['error'][1]) + np.abs(tfit['error'][1])
            bx2[i],cx2 = 2.0*bfit['beta'][2],2.0*cfit['beta'][2]
            ex2[i] = 2.0*bfit['error'][2] + np.abs(leakage_fraction*cx2) + \
                np.abs(2.0*ofit['error'][2]) + np.abs(2.0*afit['error'][2]) + \
                np.abs(2.0*tfit['error'][2])
            args = (i,bx1[i],ex1[i],bx2[i],ex2[i])
            # add to plot with colors and label
            # \u00B1 is the unicode symbol for plus-minus
            # args = (iter_label[i],bx1[i],ex1[i])
            # plot_label = u'{0}: {1:0.1f}\u00B1{2:0.1f} Gt/yr'.format(*args)
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
            ax[h].plot(tnan, mnan, color=plot_colors[i], label=plot_label, zorder=2)
            # fill between monthly errors
            ax[h].fill_between(tnan, mnan-grace_error, y2=mnan+grace_error,
                color=plot_colors[i], alpha=0.5, zorder=1)
            # converting gigatonnes to milligrams then to mm sea level
            mm_sealevel = -1e18*(gia_corrected_mean-gia_corrected_mean[0])/ocean_area
            cumulative_gia = tstar*(gia_corrected_error[-1]-gia_corrected_error[0])
            mm_sealevel_error = 1e18*np.sqrt(grace_error**2 + cumulative_gia**2)/ocean_area
            args = (reg,iter_label[i],mm_sealevel[-1],mm_sealevel_error)
            print('{0} {1}: {2:f}+/-{3:f} mm'.format(*args),file=fid)
            # add to totals (error in quadrature)
            mm_total[i] += mm_sealevel[-1]
            mm_error[i] += mm_sealevel_error**2

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax[h].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax[h].axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax[h].axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax[h].axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[h],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax[h].add_artist(at)
        data_ticks = np.arange(np.floor(tdec[0]),np.ceil(tdec[-1])+1,2)
        ax[h].xaxis.set_ticks(data_ticks)
        data_ticks = np.arange(ylimits[1]-200,ylimits[0]-ylimits[2],-ylimits[2])
        ax[h].yaxis.set_ticks(data_ticks[::-1])
        ax[h].set_xticks(np.arange(np.floor(tdec[0]),np.ceil(tdec[-1]),2))
        # axlim = ax[h].set_xlim([np.floor(tdec[0]), np.ceil(tdec[-1])])
        # axlim = ax[h].set_xlim([np.floor(tdec[0]), np.ceil(2.0*tdec[-1])/2.0])
        axlim = ax[h].set_xlim(2002, 2021.5)
        axlim = ax[h].set_ylim(ylimits[0:2])
        ax[h].get_xaxis().set_tick_params(which='both', direction='in')
        ax[h].get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker
        ax[h].xaxis.get_major_formatter().set_useOffset(False)
        # y labels on both sides of plot
        ax1[h] = ax[h].twinx()
        # add plot to hidden mm sea level axis
        ax1[h].plot(tdec, mm_sealevel, visible=False)
        ax1[h].set_ylim(-1e18*np.array(ylimits[0:2])/ocean_area)
        #ax1[h].yaxis.set_ticks(np.arange(12,-4,-1))
        ax1[h].tick_params(axis='y',colors='black',which='both',direction='in')

    # formatted ticks on N axis
    ax1['N'].yaxis.get_major_formatter().set_useOffset(False)
    ax1['N'].set_ylabel('Equivalent Sea Level Contribution [mm]',labelpad=10,color='black')
    for tl in ax1['N'].get_yticklabels():
        tl.set_color('black')
    # hidden ticks
    ax['N'].yaxis.set_ticklabels([])
    ax1['S'].yaxis.set_ticklabels([])

    # add legend
    lgd = ax['S'].legend(loc=3,frameon=False)
    # lgd = ax['S'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
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

    # add x and y labels
    ax['N'].set_xlabel('Time [Yr]')
    ax['S'].set_xlabel('Time [Yr]')
    ax['S'].set_ylabel('Mass [Gt]')
    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.1, right=0.93, bottom=0.10, top=0.98, wspace=0.06)
    figurefile = filepath.joinpath('fig5ab_{0}_{1}.pdf'.format(PROC,DREL))
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

    # print sea level totals
    for i,F in enumerate(FLAG):
        args = ('Total',iter_label[i],mm_total[i],np.sqrt(mm_error[i]))
        print('{0} {1}: {2:f}+/-{3:f} mm'.format(*args),file=fid)


    # flags for creating plots
    FLAG = ['HEX','HEX','HEX','HEX']
    SLF = ['','_SLF1','_SLF2','_SLF3']
    plot_colors = ['black','darkorchid','darkorange','mediumseagreen']
    plot_titles = ['No SL Correction','Iteration 1','Iteration 2','Iteration 3']

    # create figure axis
    ax = {}
    fig, (ax['S'],ax['N']) = plt.subplots(num=1,ncols=2,sharey=True,figsize=(8,3.5))
    # y limits for plot
    ylimits = [-90,190,20]
    for h,reg in REGION.items():
        # leakage fraction for regions
        leakage_fraction = regional_leakage[reg]
        # read ocean bottom pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,OBP_START,OBP_END)
        OBP_file = ff.format('ECCO-GAD_OBP_Residuals',reg,'','',ocean_str,LMAX,gw_str,ds_str)
        OBP_input = np.loadtxt(mascon_dir.joinpath(subdir,OBP_file))[:nmon,:]
        # read atmospheric pressure leakage file
        subdir = sd.format('AOD1B',DREL,'',LMAX,ATM_START,ATM_END)
        # ATM_file = ff.format('ATM-GAA_Residuals',reg,'_3D','',ocean_str,LMAX,gw_str)
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
                gia_Ylms = gravtk.read_GIA_model(base_dir.joinpath(*GIA_FILE), GIA=g)
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
                ax[h].plot(tnan, rnan, color=plot_colors[i], label=plot_label)
            # save mass for iteration
            mass_m1 = gia_corrected_mean - gia_corrected_mean[0]
            error_m1 = np.copy(grace_error)

        # vertical line denoting the accelerometer shutoff
        acc = gravtk.time.convert_calendar_decimal(2016,9,
            day=3,hour=12,minute=12)
        ax[h].axvline(acc,color='0.5',ls='dashed',lw=0.5,dashes=(12,6))
        # vertical lines for end of the GRACE mission and start of GRACE-FO
        jj, = np.flatnonzero(mon == 186)
        kk, = np.flatnonzero(mon == 198)
        # ax[h].axvline(tdec[jj],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        # ax[h].axvline(tdec[kk],color='0.5',ls='dashed',lw=0.5,dashes=(8,4))
        vs = ax[h].axvspan(tdec[jj],tdec[kk],color='0.5',ls='dashed',alpha=0.15)
        vs._dashes = (6,3)
        # add horizontal line at 0
        ax[h].axhline(0.0, color='black', ls='dashed', dashes=(11,5), lw=0.5)
        # add labels
        textprops = dict(size=14,weight='bold')
        at = matplotlib.offsetbox.AnchoredText(fig_text[h],
            prop=textprops, pad=0, frameon=False, loc=2)
        ax[h].add_artist(at)
        data_ticks = np.arange(np.floor(tdec[0]),np.ceil(tdec[-1])+1,2)
        ax[h].xaxis.set_ticks(data_ticks)
        data_ticks = np.arange(ylimits[1]-10,ylimits[0]-ylimits[2],-ylimits[2])
        ax[h].yaxis.set_ticks(data_ticks[::-1])
        ax[h].set_xticks(np.arange(np.floor(tdec[0]),np.ceil(tdec[-1]),2))
        # axlim = ax[h].set_xlim([np.floor(tdec[0]), np.ceil(tdec[-1])])
        # axlim = ax[h].set_xlim([np.floor(tdec[0]), np.ceil(2.0*tdec[-1])/2.0])
        axlim = ax[h].set_xlim(2002, 2021.5)
        axlim = ax[h].set_ylim(ylimits[0:2])
        ax[h].get_xaxis().set_tick_params(which='both', direction='in')
        ax[h].get_yaxis().set_tick_params(which='both', direction='in')
        # set axis ticker to integers
        ax[h].xaxis.get_major_formatter().set_useOffset(False)

    # add legend
    lgd = ax['S'].legend(loc=3,frameon=False)
    # lgd = ax['S'].legend(loc=3,frameon=False,handletextpad=-0.2,handlelength=0)
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

    # add x and y labels
    ax['N'].set_xlabel('Time [Yr]')
    ax['S'].set_xlabel('Time [Yr]')
    ax['S'].set_ylabel('Mass Difference [Gt]')
    # adjust plot to figure dimensions
    fig.subplots_adjust(left=0.08,right=0.98,bottom=0.11,top=0.98,wspace=0.06)
    figurefile = pathlib.Path('fig6ab_{0}_{1}.pdf'.format(PROC,DREL))
    plt.savefig(figurefile, format='pdf')
    plt.cla()
    plt.clf()
    plt.close()

    # close the output table
    fid.close()

# PURPOSE: run a monte carlo type error analysis with each error component
def monte_carlo_regress(t_in, d_in, satellite_error, gia_error, leakage_model,
    ocean_error, atm_error, twc_error, RUNS=10000, ORDER=0, CYCLES=[],
    TERMS=[], STDEV=0, CONF=0, RMS=True, RELATIVE=Ellipsis):

    # remove singleton dimensions
    t_in = np.squeeze(t_in)
    d_in = np.squeeze(d_in)
    nmax = len(t_in)

    # calculate epoch for calculating relative times
    if isinstance(RELATIVE, (list, np.ndarray)):
        t_rel = t_in[RELATIVE].mean()
    elif isinstance(RELATIVE, (float, int, np.float64, np.int_)):
        t_rel = np.copy(RELATIVE)
    elif (RELATIVE == Ellipsis):
        t_rel = t_in[RELATIVE].mean()

    # create design matrix based on polynomial order and harmonics
    # with any additional fit terms
    DMAT = []
    # add polynomial orders (0=constant, 1=linear, 2=quadratic)
    for o in range(ORDER+1):
        DMAT.append((t_in-t_rel)**o)
    # add cyclical terms (0.5=semi-annual, 1=annual)
    for c in CYCLES:
        DMAT.append(np.sin(2.0*np.pi*t_in/np.float64(c)))
        DMAT.append(np.cos(2.0*np.pi*t_in/np.float64(c)))
    # add additional terms to the design matrix
    for t in TERMS:
        DMAT.append(t)
    # take the transpose of the design matrix
    DMAT = np.transpose(DMAT)

    # output beta_err range (standard deviation or confidence interval)
    if (STDEV != 0):
        # Setting output error alpha with standard deviation
        alpha = 1.0 - scipy.special.erf(STDEV/np.sqrt(2.0))
    elif (CONF != 0):
        # Setting output error alpha with confidence interval
        alpha = 1.0 - CONF
    else:
        # Setting output error alpha with default 95% confidence interval
        alpha = 1.0 - (0.95)

    # number of regression terms
    n_terms = np.shape(DMAT)[1]
    # degrees of freedom of regression
    nu = nmax - n_terms
    # Student T-Distribution with D.O.F. nu
    # t.ppf parallels tinv in matlab
    tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)

    # initiating output variables from MC run
    beta_mat = np.zeros((RUNS,n_terms))# regressed variable
    beta_hat = np.zeros((n_terms))# regressed variable with max probability
    beta_conf = np.zeros((n_terms,2))# confidence interval
    beta_pdf = np.zeros((RUNS,n_terms))# regressed variable
    st_err = np.zeros((RUNS,n_terms))# standard error
    beta_err = np.zeros((RUNS,n_terms))# error to specified std or confidence

    # Calculating Least-Squares Coefficients
    # Least-Squares fitting
    # Covariance Matrix
    # Multiplying the design matrix by itself
    Hinv = np.linalg.inv(np.dot(np.transpose(DMAT),DMAT))
    # Taking the diagonal components of the cov matrix
    hdiag = np.diag(Hinv)

    for i in range(0, RUNS):
        # allocate for data_simul
        data_simul = np.zeros((nmax))
        # create two random values between 0 and 1 for GIA and leakage
        random_value = np.random.rand(2)
        # adding the data_err to each data measurement
        # the random variable will make error between +/- the data error
        # assuming data error is uniformly distributed as a worst case
        # error is between +/- data error
        gia_error_rand = (1.0 - 2.0*random_value[0])*gia_error
        leakage_error_rand = (1.0 - 2.0*random_value[1])*leakage_model
        satellite_error_rand = (1.0 - 2.0*np.random.rand(nmax))*satellite_error
        ocean_error_rand = (1.0 - 2.0*np.random.rand(nmax))*ocean_error
        atmosphere_error_rand = (1.0 - 2.0*np.random.rand(nmax))*atm_error
        total_water_error_rand = (1.0 - 2.0*np.random.rand(nmax))*twc_error
        # adding random assessments of each measurement error to data points
        data_simul[:] = d_in + gia_error_rand + leakage_error_rand + \
            satellite_error_rand + ocean_error_rand + atmosphere_error_rand + \
            total_water_error_rand

        # Standard Least-Squares fitting (the [0] denotes coefficients output)
        beta_mat[i,:] = np.linalg.lstsq(DMAT,data_simul,rcond=-1)[0]

        # MSE = (1/nu)*sum((Y-X*B)**2)
        # Mean square error (real data - model)^2/DOF
        MSE = np.dot(np.transpose(data_simul - np.dot(DMAT,beta_mat[i,:])),
            (data_simul - np.dot(DMAT,beta_mat[i,:])))/nu

        # Standard Error (1 STD)
        st_err[i,:] = np.sqrt(MSE*hdiag)
        # beta_err is the error for each coefficient
        # beta_err = t(nu,1-alpha/2)*standard error
        beta_err[i,:] = tstar*st_err[i,:]

    # Calculating the mean of the regressed terms
    # the mean trend shouldn't be much different than the standard regression
    # if the number of runs is high enough (few thousand)
    beta_mean = np.mean(beta_mat, axis=0)
    # Variance of the regression
    # Will add to the mean of the errors from the regressions
    beta_var = np.zeros((n_terms))
    for n in range(n_terms):
        beta_var[n] = np.dot(np.transpose(beta_mat[:,n] - beta_mean[n]), \
            beta_mat[:,n] - beta_mean[n])/RUNS
        # probability distribution of the regression
        p_beta = np.exp(-0.5*(beta_mat[:,n]-beta_mean[n])**2)/np.sqrt(2.0*np.pi)
        # normalize the pdf of beta
        beta_pdf[:,n] = p_beta/np.sum(p_beta)
        # beta_hat is the beta with highest probability
        # calculated by finding the max of the probability distribution
        indices = np.argmax(beta_pdf[:,n])
        beta_hat[n] = beta_mat[:,n][indices]
        # sorting the probabilities in descending order
        im = np.argsort(beta_pdf[:,n])[::-1]
        # taking the cumulative sum of the sorted probabilities
        # will find the error (to specified confidence)
        cum_beta_pdf = np.cumsum(beta_pdf[im,n])
        # minimum beta
        beta_min = np.interp(1.0-alpha,cum_beta_pdf,beta_mat[im,n])
        # confidence interval
        # beta_hat - beta_min = error
        beta_conf[n,0] = beta_hat[n]-np.abs(beta_min-beta_hat[n])
        beta_conf[n,1] = beta_hat[n]+np.abs(beta_min-beta_hat[n])

    # Propagating RMS errors (RMS is default, worst case can be specified)
    # Also will output the variance of the regressed coefficient (as beta_var)
    # Student T-Distribution with D.O.F. (RUNS-1): removing 1 for mean
    # t.ppf parallels tinv in matlab
    tstar_mean = scipy.stats.t.ppf(1.0-(alpha/2.0),RUNS-1.0)
    if RMS:
        # mean_std is the standard error
        mean_std = np.sqrt(np.sum(st_err**2.0,axis=0)/RUNS + beta_var)
        # mean_err = t(nu,1-alpha/2)*standard error
        mean_err=np.sqrt(np.sum(beta_err**2.0,axis=0)/RUNS+tstar_mean*beta_var)
    else:
        # mean_std is the standard error
        mean_std = np.mean(st_err,axis=0) + np.sqrt(beta_var)
        # error at specified standard deviation or confidence interval
        # mean_err = t(nu,1-alpha/2)*standard error
        mean_err = np.mean(beta_err,axis=0) + tstar_mean*np.sqrt(beta_var)

    return {'matrix':beta_mat, 'beta':beta_mean, 'error':mean_err, \
        'std_err':mean_std, 'variance':beta_var, 'cov_mat':Hinv, \
        'pdf':beta_pdf, 'hat':beta_hat, 'confidence':beta_conf}


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
    parser.add_argument('--scenarios','-s',
        default=False, action='store_true',
        help='Run with all sea level scenarios')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run program for parameters
    plot_mascon_SLF_iterations(args.directory,args.center,args.release,
        args.start,args.end,args.missing,SCENARIOS=args.scenarios)

# run main program
if __name__ == '__main__':
    main()
