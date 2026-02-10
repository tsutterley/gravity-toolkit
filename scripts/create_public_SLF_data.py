#!/usr/bin/env python
u"""
create_public_SLF_data.py
Written by Tyler Sutterley (10/2023)
Creates public sea level fingerprint files

UPDATE HISTORY:
    Updated 10/2023: generalize mission variable to be GRACE/GRACE-FO
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: updated with public repository functions
        updated inputs to spatial from_ascii function
    Written 08/2019
"""
from __future__ import print_function

import sys
import os
import inspect
import logging
import pathlib
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
def create_public_data(base_dir, PROC, DREL, DSET, LMAX,
    START_MON=None,
    END_MON=None,
    MISSING=None,
    GIA=None,
    GIA_FILE=None,
    REDISTRIBUTE_MASCONS=False,
    ITERATION=None,
    EXPANSION=None,
    LANDMASK=None,
    DATAFORM=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # input GIA stokes coefficients to get titles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    gia_str = f'_{GIA_Ylms_rate.title.upper()}' if GIA else ''
    # distributing mascon mass uniformly over ocean
    # mascon distribution over the ocean
    ocean_str = '_OCN' if REDISTRIBUTE_MASCONS else ''
    # version flags
    VERSION = ['v0','']
    DATA_VERSION = 1

    # Land-Sea Mask with Antarctica from Rignot (2017) and Greenland from GEUS
    # 0=Ocean, 1=Land, 2=Lake, 3=Small Island, 4=Ice Shelf
    # Open the land-sea NetCDF4 file for reading
    landsea = gravtk.spatial().from_netCDF4(LANDMASK, date=False,
        varname='LSMASK')
    dlon,dlat = landsea.spacing
    nlat, nlon = landsea.shape

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=0,
        REFERENCE='CF', FORMAT='class')

    # input directory
    FLAG={1:'',2:'_SLF1',3:'_SLF2',4:'_SLF3',5:'_SLF3',6:'_SLF3',7:'_SLF3'}
    # subdirectory and input file formats
    sd = 'HEX_{0}_{1}{2}_SPH_CAP_MSCNS{3}_L{4:d}_{5:03d}-{6:03d}'
    # input and output format
    file_format = 'SLF_ITERATION_{0}{1}{2}{3}_{4}L{5:d}{6}_{7:03d}.{8}'
    # input mascon file for GIA correction
    subdir = sd.format(PROC,DREL,VERSION[DATA_VERSION],
        FLAG[ITERATION],LMAX,START_MON,END_MON)
    # input directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    mascon_dir = base_dir.joinpath('GRACE','mascons',subdir)
    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # list of all months
    months = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    nmon = len(months)
    # output filename suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # read input sea level fingerprint files
    spatial_list = []
    for t in range(0, nmon):
        # sea level file for month (spatial fields)
        SLF = file_format.format(ITERATION, dset_str, gia_str, ocean_str, '',
            EXPANSION, '', months[t], suffix[DATAFORM])
        # read sea level file
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            dinput = gravtk.spatial().from_ascii(mascon_dir.joinpath(SLF),
                spacing=[dlon,dlat], nlat=nlat, nlon=nlon)
        elif (DATAFORM == 'netCDF4'):
            # netcdf (.nc)
            dinput = gravtk.spatial().from_netCDF4(mascon_dir.joinpath(SLF))
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            dinput = gravtk.spatial().from_HDF5(mascon_dir.joinpath(SLF))
        # append to spatial list
        spatial_list.append(dinput)

    # convert to single spatial object
    output = gravtk.spatial().from_list(spatial_list, sort=True, clear=True)

    # attributes for output file
    kwargs = {}

    # variable attributes
    kwargs['units'] = 'centimeters'
    kwargs['longname'] = 'Equivalent_Water_Thickness'
    # file-level attributes
    kwargs['attributes'] = dict(ROOT={})
    kwargs['attributes']['ROOT']['generating_institute'] = PROC
    kwargs['attributes']['ROOT']['product_release'] = DREL
    kwargs['attributes']['ROOT']['product_name'] = DSET
    kwargs['attributes']['ROOT']['title'] = 'Sea_Level_Fingerprint'
    kwargs['attributes']['ROOT']['max_degree'] = LMAX
    kwargs['attributes']['ROOT']['earth_model'] = LOVE.model
    kwargs['attributes']['ROOT']['earth_love_numbers'] = LOVE.citation
    kwargs['attributes']['ROOT']['reference_frame'] = LOVE.reference
    kwargs['attributes']['ROOT']['earth_body_tide'] = 'Wahr (1981)'
    kwargs['attributes']['ROOT']['earth_fluid_love'] = 'Han and Wahr (1989)'
    kwargs['attributes']['ROOT']['polar_motion_feedback'] = 'Kendall et al. (2005)'
    kwargs['attributes']['ROOT']['reference'] = \
        f'Output from {pathlib.Path(sys.argv[0]).name}'
    # data summary
    MISSION = 'GRACE/GRACE-FO'
    SUMMARY = []
    SUMMARY.append(('Sea level variations derived from {0} '
        'mission measurements').format(MISSION))
    SUMMARY.append(('Glacial Isostatic Adjustment (GIA) estimates from '
        f'{GIA_Ylms_rate.citation} [{GIA_Ylms_rate.title}] have been removed'))
    SUMMARY.append(('Terrestrial water storage (TWS) anomalies from {0} '
        'have been removed (Rodell et al., 2004).').format('GLDAS NOAHv2.1'))
    kwargs['attributes']['ROOT']['summary'] = '. '.join(SUMMARY)
    # data project
    PROJECT = []
    PROJECT.append('NASA Gravity Recovery And Climate Experiment (GRACE)')
    PROJECT.append('GRACE Follow-On (GRACE-FO)') if (DREL == 'RL06') else None
    kwargs['attributes']['ROOT']['project'] = ', '.join(PROJECT)
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
    kwargs['attributes']['ROOT']['keywords'] = ', '.join(KEYWORDS)
    # data keyword vocabulary
    VOCABULARY = 'NASA Global Change Master Directory (GCMD) Science Keywords'
    kwargs['attributes']['ROOT']['keywords_vocabulary'] = VOCABULARY
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
    kwargs['attributes']['ROOT']['acknowledgement'] = '. '.join(ACKNOWLEDGEMENT)
    # data version
    PRODUCT_VERSION = f'Release-{DREL[2:]}.{DATA_VERSION}'
    kwargs['attributes']['ROOT']['product_version'] = PRODUCT_VERSION
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
    kwargs['attributes']['ROOT']['references'] = '\n'.join(REFERENCE)

    # product creators and institutions
    CREATORS = 'Tyler C. Sutterley, Isabella Velicogna, and Chia-Wei Hsu'
    kwargs['attributes']['ROOT']['creator_name'] = CREATORS
    EMAILS = 'tsutterl@uw.edu, isabella@uci.edu, and chiaweih@email.arizona.edu'
    kwargs['attributes']['ROOT']['creator_email'] = EMAILS
    URL = 'https://www.ess.uci.edu/~velicogna/index.html'
    kwargs['attributes']['ROOT']['creator_url'] = URL
    kwargs['attributes']['ROOT']['creator_type'] = 'group'
    INSTITUTION = []
    INSTITUTION.append('University of Washington')
    INSTITUTION.append('University of California, Irvine')
    kwargs['attributes']['ROOT']['creator_institution'] = ', '.join(INSTITUTION)

    # output to file
    FILE = 'SLF{0}{1}{2}_L{3:d}_{4:03d}-{5:03d}.{6}'.format(dset_str, gia_str,
        ocean_str, EXPANSION, months[0], months[-1], suffix[DATAFORM])
    OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
    # save as output DATAFORM
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        # only print ocean points
        output.fill_value = 0
        output.update_mask()
        output.to_ascii(OUTPUT_FILE, date=True)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4 (.nc)
        output.to_netCDF4(OUTPUT_FILE, date=True, **kwargs)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        output.to_HDF5(OUTPUT_FILE, date=True, **kwargs)
    # set the permissions mode of the output file
    OUTPUT_FILE.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates public data for a sea level fingerprints
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
    parser.add_argument('--expansion','-e',
        type=int, default=240,
        help='Spherical harmonic expansion for sea level fingerprints')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
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
        create_public_data(
            args.directory,
            args.center,
            args.release,
            args.product,
            args.lmax,
            START_MON=args.start,
            END_MON=args.end,
            MISSING=args.missing,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            ITERATION=args.iteration,
            EXPANSION=args.expansion,
            LANDMASK=args.mask,
            DATAFORM=args.format,
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
