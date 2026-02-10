#!/usr/bin/env python
u"""
create_public_timeseries.py
Written by Tyler Sutterley (10/2023)
Creates public time series files

UPDATE HISTORY:
    Updated 10/2023: generalize mission variable to be GRACE/GRACE-FO
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: updated with public repository functions
    Written 08/2019
"""
from __future__ import print_function

import sys
import os
import io
import time
import inspect
import logging
import pathlib
import zipfile
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# program module to run with specified parameters
def create_public_timeseries(base_dir, PROC, DREL, DSET, LMAX, RAD,
    START_MON=None,
    END_MON=None,
    MISSING=None,
    MMAX=None,
    DESTRIPE=False,
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # input directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    mascon_dir = base_dir.joinpath('GRACE','mascons')
    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # Gaussian smoothing string for radius RAD (if RAD == 0: none)
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # used filtered (destriped) coefficients for stripe effects
    ds_str = '_FL' if DESTRIPE else ''
    # input GIA stokes coefficients to get titles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    gia_str = GIA_Ylms_rate.title if GIA else ''
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = 'OCN_' if REDISTRIBUTE_MASCONS else ''
    # version flags
    VERSION = ['v0','']
    DATA_VERSION = 1

    # list of all months
    months = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    nmon = len(months)

    # start and end GRACE months for correction data
    OBP_START,OBP_END = (4,254)
    ATM_START,ATM_END = (4,251)
    GLDAS_START,GLDAS_END = (4,254)

    # input directory
    FLAG={1:'_SLF1',2:'_SLF2',3:'_SLF3',4:'_SLF4',5:'_SLF5',6:'_SLF6',7:'_SLF7'}
    # subdirectory and input file formats
    sd = 'HEX_{0}_{1}{2}_SPH_CAP_MSCNS{3}_L{4:d}_{5:03d}-{6:03d}'
    ff = '{0}_{1}{2}_SPH_CAP_{3}{4}L{5:d}{6}{7}.txt'

    # data summary
    MISSION = 'GRACE/GRACE-FO'
    SUMMARY = []
    SUMMARY.append(('Regional ice mass balance time series derived from {0} '
        'mission measurements').format(MISSION))
    SUMMARY.append(('Glacial Isostatic Adjustment (GIA) estimates from '
        f'{GIA_Ylms_rate.citation} [{GIA_Ylms_rate.title}] have been removed'))
    SUMMARY.append(('Terrestrial water storage (TWS) anomalies from {0} '
        'have been removed (Rodell et al., 2004).').format('GLDAS NOAHv2.1'))
    # data project
    PROJECT = []
    PROJECT.append('NASA Gravity Recovery And Climate Experiment (GRACE)')
    PROJECT.append('GRACE Follow-On (GRACE-FO)') if (DREL == 'RL06') else None
    # data keywords
    KEYWORDS = []
    KEYWORDS.append('GRACE')
    KEYWORDS.append('GRACE-FO') if (DREL == 'RL06') else None
    # KEYWORDS.append('Level-2')
    KEYWORDS.append('Spherical Harmonic Model')
    KEYWORDS.append('Gravitational Field')
    KEYWORDS.append('Geopotential')
    KEYWORDS.append('Time Variable Gravity')
    KEYWORDS.append('Mass Transport')
    KEYWORDS.append('Satellite Geodesy')
    # data keyword vocabulary
    VOCABULARY = 'NASA Global Change Master Directory (GCMD) Science Keywords'
    # work acknowledgements
    ACKNOWLEDGEMENT = []
    ACKNOWLEDGEMENT.append(('Work was supported by an appointment to the NASA '
        'Postdoctoral Program at NASA Goddard Space Flight Center, '
        'administered by Universities Space Research Association under '
        'contract with NASA'))
    ACKNOWLEDGEMENT.append(('GRACE is a joint mission of NASA (USA) and DLR '
        '(Germany)'))
    if (DREL == 'RL06'):
        ACKNOWLEDGEMENT.append('GRACE-FO is a joint mission of NASA (USA) and '
            'GFZ (Germany)')
    # data version
    PRODUCT_VERSION = f'Release-{DREL[2:]}.{DATA_VERSION}'
    # product reference
    REFERENCE = []
    REFERENCE.append(('I. Velicogna, Y. Mohajerani, G. A, F. Landerer, '
        'J. Mouginot, B. Noel, E. Rignot and T. Sutterley, '
        '"Continuity of ice sheet mass loss in Greenland and Antarctica '
        'from the GRACE and GRACE Follow-On missions", '
        'Geophysical Research Letters, 47, (2020). '
        'https://doi.org/10.1029/2020GL087291'))
    REFERENCE.append(('T. C. Sutterley, I. Velicogna, and C.-W. Hsu, '
        '"Self-Consistent Ice Mass Balance and Regional Sea Level from '
        'Time-Variable Gravity", Earth and Space Science, 7(3), (2020). '
        'https://doi.org/10.1029/2019EA000860'))
    # append GIA reference
    if GIA_Ylms_rate.reference is not None:
        REFERENCE.append(GIA_Ylms_rate.reference)
    REFERENCE.append('M. Rodell, P. R. Houser, U. Jambor, J. Gottschalck, K. '
        'Mitchell, C.-J. Meng, K. Arsenault, B. Cosgrove, J. Radakovich, M. '
        'Bosilovich, J. K. Entin, J. P. and Walker, D. Lohmann, and D. Toll, '
        '"The Global Land Data Assimilation System." Bulletin of the American '
        'Meteorological Society, 85(3), 381-394, (2004). '
        'https://doi.org/10.1175/BAMS-85-3-381')
    # product creators and institutions
    CREATORS = 'Tyler C. Sutterley, Isabella Velicogna, and Chia-Wei Hsu'
    EMAILS = 'tsutterl@uw.edu, isabella@uci.edu, and chiaweih@email.arizona.edu'
    URL = 'https://www.ess.uci.edu/~velicogna/index.html'
    INSTITUTION = []
    INSTITUTION.append('University of Washington')
    INSTITUTION.append('University of California, Irvine')

    # open output zipfile containing regional files
    ZFILE = '{0}_SPH_CAP_{1}_L{2:d}{3}.zip'.format(gia_str,'RAD1.5',LMAX,gw_str)
    OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(ZFILE)
    zp = zipfile.ZipFile(OUTPUT_FILE, mode='w')

    # run for specific regions
    regions = []
    AIS_regions = ['AIS','WAIS','EAIS','APIS']
    GIS_regions = ['GIS','NW','NN','NE','SW','SE']
    # GIC_regions = ['CBI','CDE','ICL','SVB','FJL','SZEM','NZEM','ALK']
    HEM = ['S']*len(AIS_regions) + ['N']*len(GIS_regions) #+ ['N']*len(GIC_regions)
    remove = ['AIS']*len(AIS_regions) + ['ARC']*len(GIS_regions) #+ ['ARC']*len(GIC_regions)
    regions.extend(AIS_regions)
    regions.extend(GIS_regions)
    # regions.extend(GIC_regions)
    # for each region
    for h,reg,rem in zip(HEM,regions,remove):
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

        # input mascon file for GIA correction
        subdir = sd.format(PROC,DREL,VERSION[DATA_VERSION],
            FLAG[ITERATION],LMAX,START_MON,END_MON)
        input_file = ff.format(gia_str,reg,'',atm_str,ocean_str,LMAX,gw_str,'')
        dinput = np.loadtxt(mascon_dir.joinpath(subdir,input_file))
        mon = dinput[:nmon,0].astype(np.int64)
        tdec = dinput[:nmon,1]
        mass = dinput[:nmon,2]
        satellite_error = dinput[:nmon,3]**2
        # calculate total combined error
        grace_error = np.sqrt(np.sum(satellite_error + SLF_RMS**2 +
            OBP_RMS**2 + ATM_RMS**2 + TWC_RMS**2)/nmon)

        # open output file as in-memory object
        fid = io.StringIO()
        # write header to file
        # print header
        fid.write('{0}:\n'.format('header'))
        # data dimensions
        fid.write('  {0}:\n'.format('dimensions'))
        fid.write('    {0:22}: {1:d}\n'.format('time',nmon))
        fid.write('\n')
        fid.write('  {0}:\n'.format('global_attributes'))
        fid.write('    {0:22}: {1}\n'.format('summary','. '.join(SUMMARY)))
        fid.write('    {0:22}: {1}\n'.format('project',', '.join(PROJECT)))
        fid.write('    {0:22}: {1}\n'.format('keywords',', '.join(KEYWORDS)))
        fid.write('    {0:22}: {1}\n'.format('keywords_vocabulary',VOCABULARY))
        fid.write('    {0:22}: {1}\n'.format('acknowledgement',
            '.  '.join(ACKNOWLEDGEMENT)))
        fid.write('    {0:22}: {1}\n'.format('product_version',PRODUCT_VERSION))
        fid.write('    {0:22}:\n'.format('references'))
        for ref in REFERENCE:
            fid.write('      - {0}\n'.format(ref))
        fid.write('    {0:22}: {1}\n'.format('creator_name', CREATORS))
        fid.write('    {0:22}: {1}\n'.format('creator_email', EMAILS))
        fid.write('    {0:22}: {1}\n'.format('creator_url', URL))
        fid.write('    {0:22}: {1}\n'.format('creator_type', 'group'))
        fid.write('    {0:22}: {1}\n'.format('creator_institution',', '.join(INSTITUTION)))
        # date range and date created
        calendar_year = np.floor((mon-1)/12 + 2002.0)
        calendar_month = ((mon-1) % 12) + 1
        start_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[0],calendar_month[0])
        fid.write('    {0:22}: {1}\n'.format('time_coverage_start', start_time))
        end_time = '{0:4.0f}-{1:02.0f}'.format(calendar_year[-1],calendar_month[-1])
        fid.write('    {0:22}: {1}\n'.format('time_coverage_end', end_time))
        today = time.strftime('%Y-%m-%d',time.localtime())
        fid.write('    {0:22}: {1}\n'.format('date_created', today))
        fid.write('\n')
        # non-standard attributes
        fid.write('  {0}:\n'.format('non-standard_attributes'))
        # data format
        fid.write('    {0:22}: {1}\n'.format('formatting_string','(i4,3f11.4)'))
        fid.write('\n')
        # variables
        fid.write('  {0}:\n'.format('variables'))
        # GRACE month
        fid.write('    {0:22}:\n'.format('month'))
        long_name = 'GRACE month of each measurement epoch'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        description = 'months starting 2002-01-01'
        fid.write('      {0:20}: {1}\n'.format('description', description))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
        fid.write('      {0:20}: {1}\n'.format('units', 'month'))
        fid.write('      {0:20}: {1}\n'.format('comment', '1st column'))
        # time
        fid.write('    {0:22}:\n'.format('mid-epoch_time'))
        long_name = 'mid-date of each measurement epoch'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
        fid.write('      {0:20}: {1}\n'.format('units', 'decimal-years'))
        fid.write('      {0:20}: {1}\n'.format('comment', '2nd column'))
        # data
        fid.write('    {0:22}:\n'.format('mass'))
        fid.write('      {0:20}: {1}\n'.format('long_name', 'mass anomaly'))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
        fid.write('      {0:20}: {1}\n'.format('units', 'gigatonnes'))
        fid.write('      {0:20}: {1}\n'.format('comment', '3rd column'))
        # error
        fid.write('    {0:22}:\n'.format('mass_error'))
        fid.write('      {0:20}: {1}\n'.format('long_name', 'uncertainty'))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'single precision'))
        fid.write('      {0:20}: {1}\n'.format('units', 'gigatonnes'))
        fid.write('      {0:20}: {1}\n'.format('comment', '4th column'))
        # end of header
        fid.write('\n\n# End of YAML header\n')
        # add data
        for m in range(nmon):
            args = (mon[m],tdec[m],mass[m],grace_error)
            fid.write('{0:4d}{1:11.4f}{2:11.4f}{3:11.4f}\n'.format(*args))
        # rewind in-memory file object
        fid.seek(0)
        # write in-memory object to zip
        zp.writestr(ff.format(gia_str,reg,'','','',LMAX,gw_str,''), fid.read())

    # change output file permissions mode to MODE
    OUTPUT_FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates public data for a mascon time series
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path,
        default=filepath,
        help='Output directory for public data files')
    # GRACE/GRACE-FO data processing center
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO Level-2 data product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str, default='GSM',
        help='GRACE/GRACE-FO Level-2 data product')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-E',
        type=int, default=232,
        help='Ending GRACE/GRACE-FO month')
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,156,162,166,167,
        172,177,178,182,187,188,189,190,191,192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+', default=MISSING,
        help='Missing GRACE/GRACE-FO months')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
    # use atmospheric jump corrections from Fagiolini et al. (2015)
    parser.add_argument('--atm-correction',
        default=False, action='store_true',
        help='Apply atmospheric jump correction coefficients')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # sea level fingerprint parameters
    parser.add_argument('--iteration','-I',
        type=int, default=1,
        help='Sea level fingerprint iteration')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run program
        create_public_timeseries(
            args.directory,
            args.center,
            args.release,
            args.product,
            args.lmax,
            args.radius,
            START_MON=args.start,
            END_MON=args.end,
            MISSING=args.missing,
            MMAX=args.mmax,
            DESTRIPE=args.destripe,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
