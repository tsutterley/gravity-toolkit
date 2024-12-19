#!/usr/bin/env python
u"""
dealiasing_global_uplift.py
Written by Tyler Sutterley (05/2023)

Reads GRACE/GRACE-FO level-1b dealiasing data files for global atmospheric
and oceanic loading and estimates anomalies in elastic crustal uplift
    atm: atmospheric loading from ECMWF
    ocn: oceanic loading from OMCT/MPIOM
    glo: global atmospheric and oceanic loading
    obp: ocean bottom pressure from OMCT/MPIOM

CALLING SEQUENCE:
    python dealiasing_global_uplift.py --release RL06 --product glo

COMMAND LINE OPTIONS:
    -D X, --directory X: Working Data Directory
    -O X, --output-directory X: output directory for spatial files
    -r X, --release X: GRACE/GRACE-FO Data Release (RL05 or RL06)
    -p X, --product X: GRACE/GRACE-FO dealiasing product (atm, ocn, glo, oba)
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -S X, --spacing X: spatial resolution of output data (dlon,dlat)
    -I X, --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    -B X, --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    -F X, --format X: input and output data format
        ascii
        netCDF4
        HDF5
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        http://www.h5py.org/
    future: Compatibility layer between Python 2 and Python 3
        http://python-future.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    spatial.py: spatial data class for reading, writing and processing data
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: attributes from units class for output netCDF4/HDF5 files
    Written 03/2023
"""
from __future__ import print_function, division

import sys
import os
import re
import gzip
import time
import logging
import pathlib
import tarfile
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: estimates global elastic uplift due to changes in atmospheric
# and oceanic loading
def dealiasing_global_uplift(base_dir,
    DREL=None,
    DSET=None,
    YEAR=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    DATAFORM=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # input directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    grace_dir = base_dir.joinpath('AOD1B', DREL)
    # output directory setup
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # list object of output files for file logs (full path)
    output_files = []

    # set number of hours in a file for a release
    # set the atmospheric and ocean model for a given release
    # set the maximum degree and order of a release
    if DREL in ('RL01','RL02','RL03','RL04','RL05'):
        # for 00, 06, 12 and 18
        nt = 4
        ATMOSPHERE = 'ECMWF'
        OCEAN_MODEL = 'OMCT'
        LMAX = 100
    elif DREL in ('RL06',):
        # for 00, 03, 06, 09, 12, 15, 18 and 21
        nt = 8
        ATMOSPHERE = 'ECMWF'
        OCEAN_MODEL = 'MPIOM'
        LMAX = 180
    else:
        raise ValueError('Invalid data release')
    # Calculating the number of cos and sin harmonics up to LMAX
    n_harm = (LMAX**2 + 3*LMAX)//2 + 1

    # AOD1B data products
    product = {}
    product['atm'] = f'Atmospheric loading from {ATMOSPHERE}'
    product['ocn'] = f'Oceanic loading from {OCEAN_MODEL}'
    product['glo'] = 'Global atmospheric and oceanic loading'
    product['oba'] = f'Ocean bottom pressure from {OCEAN_MODEL}'

    # output spatial units
    UNITS = 'mmCU'
    # attributes for output files
    attributes = dict(ROOT={})
    attributes['ROOT']['title'] = product[DSET]
    attributes['ROOT']['project_version'] = DREL
    attributes['ROOT']['product_name'] = DSET
    attributes['ROOT']['product_type'] = 'gravity_field'
    attributes['ROOT']['reference'] = \
        f'Output from {pathlib.Path(sys.argv[0]).name}'
    # output suffix for data formats
    suffix = dict(ascii='txt',netCDF4='nc',HDF5='H5')

    # compile regular expressions operators for file dates
    # will extract the year and month from the tar file (.tar.gz)
    regex_year = r'\d{4}' if YEAR is None else r'|'.join(f'{y:d}' for y in YEAR)
    tx = re.compile(rf'AOD1B_({regex_year})-(\d+)_\d+\.(tar\.gz|tgz)$', re.X)
    # and the calendar day from the ascii file (.asc or gzipped .asc.gz)
    fx = re.compile(r'AOD1B_\d+-\d+-(\d+)_X_\d+.asc(.gz)?$', re.X)
    # compile regular expressions operator for the clm/slm headers
    # for the specific AOD1b product
    hx = re.compile(rf'^DATA.*SET.*{DSET}', re.X)
    # compile regular expression operator to find numerical instances
    # will extract the data from the file
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?)|NaN)(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    # finding all of the tar files in the AOD1b directory
    input_tar_files = [tf for tf in grace_dir.iterdir() if tx.match(tf.name)]

    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output spatial data
    grid = gravtk.spatial()
    # Output Degree Interval
    if (INTERVAL == 1):
        # (0:360,90:-90)
        n_lon = np.int64((360.0/dlon)+1.0)
        n_lat = np.int64((180.0/dlat)+1.0)
        grid.lon = dlon*np.arange(0,n_lon)
        grid.lat = 90.0 - dlat*np.arange(0,n_lat)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        grid.lon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        n_lon = len(grid.lon)
        n_lat = len(grid.lat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0, maxlon+dlon/2.0, dlon)
        grid.lat = np.arange(maxlat-dlat/2.0, minlat-dlat/2.0, -dlat)
        n_lon = len(grid.lon)
        n_lat = len(grid.lat)

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth parameters
    attributes['ROOT']['earth_model'] = LOVE.model
    attributes['ROOT']['earth_love_numbers'] = LOVE.citation
    attributes['ROOT']['reference_frame'] = LOVE.reference
    attributes['ROOT']['max_degree'] = LMAX
    # add attributes for earth parameters
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE)
    attributes['ROOT']['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['ROOT']['earth_density'] = f'{factors.rho_e:0.3f} g/cm^3'
    attributes['ROOT']['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'
    # degree dependent factors for converting to output units
    dfactor = factors.get(UNITS)

    # attributes for output variables
    # Defining attributes for longitude and latitude
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    # Defining attributes for dataset
    attributes['z'] = {}
    units_name, units_longname = gravtk.units.get_attributes(UNITS)
    attributes['z']['long_name'] = units_longname
    attributes['z']['units'] = units_name
    # Defining attributes for date
    attributes['time'] = {}
    attributes['time']['long_name'] = 'time'
    attributes['time']['calendar'] = 'standard'
    attributes['time']['standard_name'] = 'time'

    # Computing plms for converting to spatial domain
    theta = (90.0 - grid.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # for each tar file
    for input_file in sorted(input_tar_files):
        # extract the year and month from the file
        YY,MM,SFX = tx.findall(input_file.name).pop()
        # number of days per month
        dpm = gravtk.time.calendar_days(int(YY))
        # output monthly spatial file
        FILE = f'AOD1B_{DREL}_{DSET}_{UNITS}_{YY}_{MM}.{suffix[DATAFORM]}'
        output_file = OUTPUT_DIRECTORY.joinpath(FILE)
        # if output file exists: check if input tar file is newer
        TEST = False
        # check if output file exists
        if output_file.exists():
            # check last modification time of input and output files
            input_mtime = input_file.stat().st_mtime
            output_mtime = output_file.stat().st_mtime
            # if input tar file is newer: overwrite the output file
            if (input_mtime > output_mtime):
                TEST = True
        else:
            TEST = True

        # As there are so many files.. this will only read the new files
        # or will rewrite if CLOBBER is set (if wanting something changed)
        if TEST:
            # track tar file
            logging.debug(str(input_file))
            # open the AOD1B monthly tar file
            tar = tarfile.open(name=str(input_file), mode='r:gz')
            # number of time points
            n_time = int(nt*dpm[int(MM)-1])
            # flattened harmonics object
            YLMS = gravtk.harmonics(lmax=LMAX, mmax=LMAX,
                flattened=True)
            YLMS.l = np.zeros((n_harm), dtype=int)
            YLMS.m = np.zeros((n_harm), dtype=int)
            YLMS.clm = np.zeros((n_harm, n_time))
            YLMS.slm = np.zeros((n_harm, n_time))
            # calendar dates
            years = np.zeros((n_time), dtype=int)
            months = np.zeros((n_time), dtype=int)
            days = np.zeros((n_time), dtype=int)
            hours = np.zeros((n_time), dtype=int)
            # Iterate over every member within the tar file
            for member in tar.getmembers():
                # track tar file members
                logging.debug(member.name)
                # get calendar day from file
                DD,SFX = fx.findall(member.name).pop()
                # open data file for day
                if (SFX == '.gz'):
                    fid = gzip.GzipFile(fileobj=tar.extractfile(member))
                else:
                    fid = tar.extractfile(member)
                # create counter for hour in dataset
                c = 0
                # while loop ends when dataset is read
                while (c < nt):
                    # read line
                    file_contents = fid.readline().decode('ISO-8859-1')
                    # find file header for data product
                    if bool(hx.search(file_contents)):
                        # track file header lines
                        logging.debug(file_contents)
                        # extract hour from header
                        HH, = re.findall(r'(\d+):\d+:\d+',file_contents)
                        # convert dates to int and save to arrays
                        i = (int(DD)-1)*nt + c
                        years[i] = np.int64(YY)
                        months[i] = np.int64(MM)
                        days[i] = np.int64(DD)
                        hours[i] = np.int64(HH)
                        # read each line of spherical harmonics
                        for k in range(0,n_harm):
                            file_contents = fid.readline().decode('ISO-8859-1')
                            # find numerical instances in the data line
                            line_contents = rx.findall(file_contents)
                            # spherical harmonic degree and order
                            YLMS.l[k] = np.int64(line_contents[0])
                            YLMS.m[k] = np.int64(line_contents[1])
                            # extract spherical harmonics
                            YLMS.clm[k,i] = np.float64(line_contents[2])
                            YLMS.slm[k,i] = np.float64(line_contents[3])
                        # add 1 to hour counter
                        c += 1
                # close the input file for day
                fid.close()
            # calculate times for flattened harmonics
            YLMS.time = gravtk.time.convert_calendar_decimal(
                years, months, day=days, hour=hours)
            YLMS.month = gravtk.time.calendar_to_grace(YLMS.time)
            # convert to expanded form in output units
            Ylms = YLMS.expand(date=True).convolve(dfactor)
            # convert harmonics to spatial domain
            grid.data = np.zeros((n_lat,n_lon,n_time))
            grid.mask = np.zeros((n_lat,n_lon,n_time), dtype=bool)
            # calculate delta times for output spatial grids
            grid.time = np.array(hours + 24*(days-1), dtype=int)
            # for each date in the harmonics object
            for i,iYlm in enumerate(Ylms):
                # convert to spatial domain
                grid.data[:,:,i] = gravtk.harmonic_summation(
                    iYlm.clm, iYlm.slm, grid.lon, grid.lat,
                    LMAX=LMAX, PLM=PLM).T
            # update attributes for time
            attributes['time']['units'] = \
                f'hours since {YY}-{MM}-01T00:00:00'
            # output spatial data to file
            grid.to_file(output_file, format=DATAFORM,
                attributes=attributes)
            # set the permissions mode of the output file
            output_file.chmod(mode=MODE)
            # append output file to list
            output_files.append(output_file)
            # close the tar file
            tar.close()

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the AOD1b spatial analysis
def output_log_file(input_arguments, output_files):
    # format: aod1b_spatial_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'aod1b_spatial_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the AOD1b spatial analysis
def output_error_log_file(input_arguments):
    # format: aod1b_spatial_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'aod1b_spatial_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads GRACE/GRACE-FO level-1b dealiasing data files
            for global atmospheric and oceanic loading and estimates anomalies
            in elastic crustal uplift
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Output directory for spatial files')
    # GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='',
        help='GRACE/GRACE-FO Data Release')
    # GRACE/GRACE-FO level-1b dealiasing product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str.lower, default='glo',
        choices=['atm','ocn','glo','oba'],
        help='GRACE/GRACE-FO Level-1b data product')
    # years to run
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,2024),
        help='Years of data to run')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds','-B',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # Output log file for each job in forms
    # aod1b_spatial_run_2002-04-01_PID-00000.log
    # aod1b_spatial_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the output files (octal)
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
        # run AOD1b uplift program with parameters
        output_files = dealiasing_global_uplift(args.directory,
            DREL=args.release,
            DSET=args.product,
            YEAR=args.year,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            DATAFORM=args.format,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
