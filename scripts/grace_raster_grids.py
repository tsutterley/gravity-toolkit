#!/usr/bin/env python
u"""
grace_raster_grids.py
Written by Tyler Sutterley (06/2024)

Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports
    projected spatial fields

Will correct with the specified GIA model group, destripe/smooth/process,
    and export the data in specified units

Spatial output units: cm w.e., mm geoid height, mm elastic uplift,
    microgal gravity perturbation or surface pressure (mbar)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for raster files
    -P X, --file-prefix X: prefix string for input and output files
    -c X, --center X: GRACE/GRACE-FO processing center
    -r X, --release X: GRACE/GRACE-FO data release
    -p X, --product X: GRACE/GRACE-FO Level-2 data product
    -S X, --start X: starting GRACE/GRACE-FO month
    -E X, --end X: ending GRACE/GRACE-FO month
    -N X, --missing X: Missing GRACE/GRACE-FO months
    --lmin X: minimum spherical harmonic degree
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use decorrelation filter (destriping filter)
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
    -F X, --format X: input/output data format
        netCDF4
        HDF5
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        AW13-IJ05: Geruo A IJ05-R2 GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    --atm-correction: Apply atmospheric jump correction coefficients
    --pole-tide: Correct for pole tide drift
    --geocenter X: Update Degree 1 coefficients with SLR or derived values
        Tellus: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
        SLR: satellite laser ranging coefficients from CSR
        UCI: Sutterley and Velicogna coefficients, Remote Sensing (2019)
        Swenson: GRACE-derived coefficients from Sean Swenson
        GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    --interpolate-geocenter: Least-squares model missing Degree 1 coefficients
    --slr-c20 X: Replace C20 coefficients with SLR values
        CSR: use values from CSR (TN-07,TN-09,TN-11)
        GFZ: use values from GFZ
        GSFC: use values from GSFC (TN-14)
    --slr-21 X: Replace C21 and S21 coefficients with SLR values
        CSR: use values from CSR
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC
    --slr-22 X: Replace C22 and S22 coefficients with SLR values
        CSR: use values from CSR
    --slr-c30 X: Replace C30 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GFZ: use values from GFZ GravIS
        GSFC: use values from GSFC (TN-14)
    --slr-c40 X: Replace C40 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    --slr-c50 X: Replace C50 coefficients with SLR values
        CSR: use values from CSR (5x5 with 6,1)
        GSFC: use values from GSFC
    -U X, --units X: output units
        1: cm of water thickness
        2: mm of geoid height
        3: mm of elastic crustal deformation [Davis 2004]
        4: microGal gravitational perturbation
        5: mbar equivalent surface pressure
    --spacing X: output grid spacing
    --bounds X: output grid extents [xmin,xmax,ymin,ymax]
    --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
        ascii
        netCDF4
        HDF5
        gfc
    --mask X: Land-sea mask for redistributing land water flux
    --remove-file X: Monthly files to be removed from the GRACE/GRACE-FO data
    --remove-format X: Input data format for files to be removed
        ascii
        netCDF4
        HDF5
        index-ascii
        index-netCDF4
        index-HDF5
    --redistribute-removed: redistribute removed mass fields over the ocean
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified spherical
            harmonic degree and order and for a specified date range
        Includes degree 1 with with Swenson values (if specified)
        Replaces low-degree harmonics with SLR values (if specified)
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: converts a land-sea mask to a series of spherical harmonics
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    geocenter.py: converts between spherical harmonics and geocenter variations
    clenshaw_summation.py: calculate spatial field from spherical harmonics
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 06/2024: use wrapper to importlib for optional dependencies
    Updated 03/2024: increase mask buffer to twice the smoothing radius
    Written 08/2023
"""
from __future__ import print_function

import sys
import os
import time
import logging
import pathlib
import numpy as np
import argparse
import traceback
import collections
import gravity_toolkit as gravtk

# attempt imports
geoidtk = gravtk.utilities.import_dependency('geoid_toolkit')
pyproj = gravtk.utilities.import_dependency('pyproj')

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: try to get the projection information
def get_projection(PROJECTION):
    # EPSG projection code
    try:
        crs = pyproj.CRS.from_epsg(int(PROJECTION))
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # no projection can be made
    raise pyproj.exceptions.CRSError

# PURPOSE: import GRACE/GRACE-FO files for a given months range
# Converts the GRACE/GRACE-FO harmonics applying the specified procedures
def grace_raster_grids(base_dir, PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    LMIN=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DESTRIPE=False,
    UNITS=None,
    BOUNDS=[],
    SPACING=[],
    PROJECTION='4326',
    GIA=None,
    GIA_FILE=None,
    ATM=False,
    POLE_TIDE=False,
    DEG1=None,
    DEG1_FILE=None,
    MODEL_DEG1=False,
    SLR_C20=None,
    SLR_21=None,
    SLR_22=None,
    SLR_C30=None,
    SLR_C40=None,
    SLR_C50=None,
    DATAFORM=None,
    MEAN_FILE=None,
    MEANFORM=None,
    REMOVE_FILES=None,
    REMOVE_FORMAT=None,
    REDISTRIBUTE_REMOVED=False,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    MODE=0o775):

    # recursively create output directory if not currently existing
    OUTPUT_DIRECTORY = pathlib.Path(OUTPUT_DIRECTORY).expanduser().absolute()
    if not OUTPUT_DIRECTORY.exists():
        OUTPUT_DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output attributes for raster files
    attributes = dict(ROOT=collections.OrderedDict())
    attributes['ROOT']['generating_institute'] = PROC
    attributes['ROOT']['product_release'] = DREL
    attributes['ROOT']['product_name'] = DSET
    attributes['ROOT']['product_type'] = 'gravity_field'
    attributes['ROOT']['title'] = 'GRACE/GRACE-FO Spatial Data'
    attributes['ROOT']['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # list object of output files for file logs (full path)
    output_files = []

    # file information
    suffix = dict(netCDF4='nc', HDF5='H5')[DATAFORM]

    # read arrays of kl, hl, and ll Love Numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth model and love numbers
    attributes['ROOT']['earth_model'] = LOVE.model
    attributes['ROOT']['earth_love_numbers'] = LOVE.citation
    attributes['ROOT']['reference_frame'] = LOVE.reference

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        gw_str = f'_r{RAD:0.0f}km'
        attributes['ROOT']['smoothing_radius'] = f'{RAD:0.0f} km'
    else:
        # else = 1
        gw_str = ''

    # flag for spherical harmonic order
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # add attributes for LMAX and MMAX
    attributes['ROOT']['max_degree'] = LMAX
    attributes['ROOT']['max_order'] = MMAX

    # reading GRACE months for input date range
    # replacing low-degree harmonics with SLR values if specified
    # include degree 1 (geocenter) harmonics if specified
    # correcting for Pole-Tide and Atmospheric Jumps if specified
    Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX, SLR_21=SLR_21,
        SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C40=SLR_C40, SLR_C50=SLR_C50,
        DEG1_FILE=DEG1_FILE, MODEL_DEG1=MODEL_DEG1, ATM=ATM,
        POLE_TIDE=POLE_TIDE)
    # convert to harmonics object and remove mean if specified
    GRACE_Ylms = gravtk.harmonics().from_dict(Ylms)
    nt = len(GRACE_Ylms.time)
    # add attributes for input GRACE/GRACE-FO spherical harmonics
    for att_name, att_val in Ylms['attributes'].items():
        attributes['ROOT'][att_name] = att_val

    # use a mean file for the static field to remove
    if MEAN_FILE:
        # read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        MEAN_FILE = pathlib.Path(MEAN_FILE).expanduser().absolute()
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=MEANFORM, date=False)
        # remove the input mean
        GRACE_Ylms.subtract(mean_Ylms)
        attributes['ROOT']['lineage'].append(MEAN_FILE.name)
    else:
        GRACE_Ylms.mean(apply=True)

    # filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        # destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
        attributes['ROOT']['filtering'] = 'Destriped'
    else:
        # using standard GRACE/GRACE-FO harmonics
        ds_str = ''

    # input GIA spherical harmonic datafiles
    GIA_FILE = pathlib.Path(GIA_FILE).expanduser().absolute() if GIA else None
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
    # output GIA string for filename
    if GIA:
        gia_str = f'_{GIA_Ylms_rate.title}'
        attributes['ROOT']['GIA'] = (str(GIA_Ylms_rate.citation), GIA_FILE.name)
    else:
        gia_str = ''
    # monthly GIA calculated by gia_rate*time elapsed
    # finding change in GIA each month
    GIA_Ylms = GIA_Ylms_rate.drift(GRACE_Ylms.time, epoch=2003.3)
    GIA_Ylms.month[:] = np.copy(GRACE_Ylms.month)

    # default file prefix
    if not FILE_PREFIX:
        fargs = (PROC,DREL,DSET,Ylms['title'],gia_str)
        FILE_PREFIX = '{0}_{1}_{2}{3}{4}_'.format(*fargs)

    # read Land-Sea Mask and convert to spherical harmonics
    land_Ylms = gravtk.land_stokes(LANDMASK, LMAX,
        MMAX=MMAX, LOVE=LOVE)
    # Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_REMOVED:
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX,
            MMAX=MMAX, LOVE=LOVE)
        ocean_str = '_OCN'
    else:
        ocean_str = ''

    # input spherical harmonic datafiles to be removed from the GRACE data
    # Remove sets of Ylms from the GRACE data before returning
    remove_Ylms = GRACE_Ylms.zeros_like()
    remove_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    remove_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    if REMOVE_FILES:
        # extend list if a single format was entered for all files
        if len(REMOVE_FORMAT) < len(REMOVE_FILES):
            REMOVE_FORMAT = REMOVE_FORMAT*len(REMOVE_FILES)
        # for each file to be removed
        for REMOVE_FILE,REMOVEFORM in zip(REMOVE_FILES,REMOVE_FORMAT):
            if REMOVEFORM in ('ascii','netCDF4','HDF5'):
                # ascii (.txt)
                # netCDF4 (.nc)
                # HDF5 (.H5)
                Ylms = gravtk.harmonics().from_file(REMOVE_FILE,
                    format=REMOVEFORM)
                attributes['ROOT']['lineage'].append(Ylms.name)
            elif REMOVEFORM in ('index-ascii','index-netCDF4','index-HDF5'):
                # read from index file
                _,removeform = REMOVEFORM.split('-')
                # index containing files in data format
                Ylms = gravtk.harmonics().from_index(REMOVE_FILE,
                    format=removeform)
                attributes['ROOT']['lineage'].extend([f.name for f in Ylms.filename])
            # reduce to GRACE/GRACE-FO months and truncate to degree and order
            Ylms = Ylms.subset(GRACE_Ylms.month).truncate(lmax=LMAX,mmax=MMAX)
            # distribute removed Ylms uniformly over the ocean
            if REDISTRIBUTE_REMOVED:
                # calculate ratio between total removed mass and
                # a uniformly distributed cm of water over the ocean
                ratio = Ylms.clm[0,0,:]/ocean_Ylms.clm[0,0]
                # for each spherical harmonic
                for m in range(0,MMAX+1):# MMAX+1 to include MMAX
                    for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                        # remove the ratio*ocean Ylms from Ylms
                        # note: x -= y is equivalent to x = x - y
                        Ylms.clm[l,m,:] -= ratio*ocean_Ylms.clm[l,m]
                        Ylms.slm[l,m,:] -= ratio*ocean_Ylms.slm[l,m]
            # filter removed coefficients
            if DESTRIPE:
                Ylms = Ylms.destripe()
            # add data for month t and INDEX_FILE to the total
            # remove_clm and remove_slm matrices
            # redistributing the mass over the ocean if specified
            remove_Ylms.add(Ylms)

    # converting x,y from projection to latitude/longitude
    crs1 = get_projection(PROJECTION)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    # dictionary of coordinate reference system variables
    crs_to_dict = crs1.to_dict()
    # Climate and Forecast (CF) Metadata Conventions
    if (crs1.to_epsg() == 4326):
        y_cf,x_cf = crs1.cs_to_cf()
    else:
        x_cf,y_cf = crs1.cs_to_cf()

    # output spatial units
    # Setting units factor for output
    # dfactor computes the degree dependent coefficients
    factors = gravtk.units(lmax=LMAX)
    # 1: cmwe, centimeters water equivalent
    # 2: mmGH, millimeters geoid height
    # 3: mmCU, millimeters elastic crustal deformation
    # 4: micGal, microGal gravity perturbations
    # 5: mbar, millibars equivalent surface pressure
    units = gravtk.units.bycode(UNITS)
    # output spatial units and descriptive units longname
    units_name, units_longname = gravtk.units.get_attributes(units)
    # add attributes for earth parameters
    attributes['ROOT']['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['ROOT']['earth_density'] = f'{factors.rho_e:0.3f} g/cm^3'
    attributes['ROOT']['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'

    # projection attributes
    attributes['crs'] = {}
    # add projection attributes
    attributes['crs']['standard_name'] = \
        crs1.to_cf()['grid_mapping_name'].title()
    attributes['crs']['spatial_epsg'] = crs1.to_epsg()
    attributes['crs']['spatial_ref'] = crs1.to_wkt()
    attributes['crs']['proj4_params'] = crs1.to_proj4()
    for att_name,att_val in crs1.to_cf().items():
        attributes['crs'][att_name] = att_val
    if ('lat_0' in crs_to_dict.keys() and (crs1.to_epsg() != 4326)):
        attributes['crs']['latitude_of_projection_origin'] = \
            crs_to_dict['lat_0']
    # x and y
    attributes['x'],attributes['y'] = ({},{})
    for att_name in ['long_name','standard_name','units']:
        attributes['x'][att_name] = x_cf[att_name]
        attributes['y'][att_name] = y_cf[att_name]
    # time
    attributes['time'] = {}
    attributes['time']['units'] = 'years'
    attributes['time']['long_name'] = 'Date_in_Decimal_Years'
    # output gridded data
    fill_value = -9999.0
    attributes['z'] = {}
    attributes['z']['units'] = units_name
    attributes['z']['long_name'] = units_longname
    attributes['z']['degree_of_truncation'] = LMAX
    attributes['z']['_FillValue'] = fill_value
    # set grid mapping attribute
    attributes['z']['grid_mapping'] = 'crs'

    # output data variables
    output = {}
    # projection variable
    output['crs'] = np.array((),dtype=np.byte)
    # spacing and bounds of output grid
    dx,dy = np.broadcast_to(np.atleast_1d(SPACING),(2,))
    xmin,xmax,ymin,ymax = np.copy(BOUNDS)
    # create x and y from spacing and bounds
    output['x'] = np.arange(xmin + dx/2.0, xmax + dx, dx)
    output['y'] = np.arange(ymin + dx/2.0, ymax + dy, dy)
    ny,nx = (len(output['y']),len(output['x']))
    gridx, gridy = np.meshgrid(output['x'],output['y'])
    gridlon, gridlat = transformer.transform(gridx, gridy)

    # semimajor axis of ellipsoid [m]
    a_axis = crs1.ellipsoid.semi_major_metre
    # ellipsoidal flattening
    flat = 1.0/crs1.ellipsoid.inverse_flattening
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = geoidtk.spatial.geocentric_latitude(
        gridlon, gridlat, a_axis=a_axis, flat=flat)

    # calculate spatial mask with an extended radius
    THRESHOLD = 0.025
    mask = gravtk.clenshaw_summation(land_Ylms.clm, land_Ylms.slm,
        gridlon.flatten(), latitude_geocentric.flatten(), RAD=2*RAD,
        UNITS=1, LMAX=LMAX, LOVE=LOVE)
    ii,jj = np.nonzero(mask.reshape(ny,nx) > THRESHOLD)

    # output gridded raster data
    output['z'] = np.ma.zeros((ny,nx,nt), fill_value=fill_value)
    output['z'].mask = np.ones((ny,nx,nt), dtype=bool)
    output['time'] = np.zeros((nt))

    # converting harmonics to truncated, smoothed coefficients in units
    # combining harmonics to calculate output raster grids
    for i,grace_month in enumerate(GRACE_Ylms.month):
        logging.debug(grace_month)
        # GRACE/GRACE-FO harmonics for time t
        Ylms = GRACE_Ylms.index(i)
        # Remove GIA rate for time
        Ylms.subtract(GIA_Ylms.index(i))
        # Remove monthly files to be removed
        Ylms.subtract(remove_Ylms.index(i))
        # truncate to degree and order LMAX and MMAX
        # truncate minimum degree to LMIN
        Ylms.truncate(LMAX, lmin=LMIN, mmax=MMAX)
        # convert spherical harmonics to output raster grid
        output['z'].data[ii,jj,i] = gravtk.clenshaw_summation(
            Ylms.clm, Ylms.slm, gridlon[ii,jj], latitude_geocentric[ii,jj],
            RAD=RAD, UNITS=units, LMAX=LMAX, LOVE=LOVE)
        output['z'].mask[ii,jj,i] = False
        # copy time variables for month
        output['time'][i] = np.copy(Ylms.time)
    # convert masked values to fill value
    output['z'].data[output['z'].mask] = output['z'].fill_value

    # output raster files to netCDF4 or HDF5
    FILE = (f'{FILE_PREFIX}{units}_L{LMAX:d}{order_str}{gw_str}{ds_str}_'
        f'{START:03d}-{END:03d}.{suffix}')
    output_file = OUTPUT_DIRECTORY.joinpath(FILE)
    # use spatial functions from geoid toolkit to write rasters
    if (DATAFORM == 'netCDF4'):
        geoidtk.spatial.to_netCDF4(output, attributes, output_file,
            data_type='grid')
    elif (DATAFORM == 'HDF5'):
        geoidtk.spatial.to_HDF5(output, attributes, output_file,
            data_type='grid')
    # set the permissions mode of the output files
    output_file.chmod(mode=MODE)
    # add file to list
    output_files.append(output_file)

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the GRACE analysis
def output_log_file(input_arguments, output_files):
    # format: GRACE_processing_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_run_{0}_PID-{1:d}.log'.format(*args)
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

# PURPOSE: print a error file log for the GRACE analysis
def output_error_log_file(input_arguments):
    # format: GRACE_processing_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GRACE_processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Calculates monthly spatial raster grids from
            GRACE/GRACE-FO spherical harmonic coefficients
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
        help='Output directory for raster files')
    parser.add_argument('--file-prefix','-P',
        type=str,
        help='Prefix string for input and output files')
    # Data processing center or satellite mission
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
    # minimum spherical harmonic degree
    parser.add_argument('--lmin',
        type=int, default=1,
        help='Minimum spherical harmonic degree')
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
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
    # output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,3,4,5],
        help='Output units')
    # output grid spacing
    parser.add_argument('--spacing',
        type=float, default=1.0, nargs='+',
        help='Output grid spacing')
    # bounds of output grid
    parser.add_argument('--bounds', type=float,
        nargs=4, default=[-180.0,180.0,-90.0,90.0],
        metavar=('xmin','xmax','ymin','ymax'),
        help='Output grid extents')
    # spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
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
    # correct for pole tide drift follow Wahr et al. (2015)
    parser.add_argument('--pole-tide',
        default=False, action='store_true',
        help='Correct for pole tide drift')
    # Update Degree 1 coefficients with SLR or derived values
    # Tellus: GRACE/GRACE-FO TN-13 from PO.DAAC
    #     https://grace.jpl.nasa.gov/data/get-data/geocenter/
    # SLR: satellite laser ranging from CSR
    #     ftp://ftp.csr.utexas.edu/pub/slr/geocenter/
    # UCI: Sutterley and Velicogna, Remote Sensing (2019)
    #     https://www.mdpi.com/2072-4292/11/18/2108
    # Swenson: GRACE-derived coefficients from Sean Swenson
    #     https://doi.org/10.1029/2007JB005338
    # GFZ: GRACE/GRACE-FO coefficients from GFZ GravIS
    #     http://gravis.gfz-potsdam.de/corrections
    parser.add_argument('--geocenter',
        metavar='DEG1', type=str,
        choices=['Tellus','SLR','SLF','UCI','Swenson','GFZ'],
        help='Update Degree 1 coefficients with SLR or derived values')
    parser.add_argument('--geocenter-file',
        type=pathlib.Path,
        help='Specific geocenter file if not default')
    parser.add_argument('--interpolate-geocenter',
        default=False, action='store_true',
        help='Least-squares model missing Degree 1 coefficients')
    # replace low degree harmonics with values from Satellite Laser Ranging
    parser.add_argument('--slr-c20',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C20 coefficients with SLR values')
    parser.add_argument('--slr-21',
        type=str, default=None, choices=['CSR','GFZ','GSFC'],
        help='Replace C21 and S21 coefficients with SLR values')
    parser.add_argument('--slr-22',
        type=str, default=None, choices=['CSR','GSFC'],
        help='Replace C22 and S22 coefficients with SLR values')
    parser.add_argument('--slr-c30',
        type=str, default=None, choices=['CSR','GFZ','GSFC','LARES'],
        help='Replace C30 coefficients with SLR values')
    parser.add_argument('--slr-c40',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C40 coefficients with SLR values')
    parser.add_argument('--slr-c50',
        type=str, default=None, choices=['CSR','GSFC','LARES'],
        help='Replace C50 coefficients with SLR values')
    # input data format (netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['netCDF4','HDF5'],
        help='Input/output data format')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=pathlib.Path,
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    # input data format for mean file (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    # monthly files to be removed from the GRACE/GRACE-FO data
    parser.add_argument('--remove-file',
        type=pathlib.Path, nargs='+',
        help='Monthly files to be removed from the GRACE/GRACE-FO data')
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--remove-format',
        type=str, nargs='+', choices=choices,
        help='Input data format for files to be removed')
    parser.add_argument('--redistribute-removed',
        default=False, action='store_true',
        help='Redistribute removed mass fields over the ocean')
    # land-sea mask for redistributing fluxes
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=pathlib.Path, default=lsmask,
        help='Land-sea mask for redistributing land water flux')
    # Output log file for each job in forms
    # GRACE_processing_run_2002-04-01_PID-00000.log
    # GRACE_processing_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
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
        # run grace_raster_grids algorithm with parameters
        output_files = grace_raster_grids(
            args.directory,
            args.center,
            args.release,
            args.product,
            args.lmax,
            args.radius,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            LMIN=args.lmin,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DESTRIPE=args.destripe,
            UNITS=args.units,
            SPACING=args.spacing,
            BOUNDS=args.bounds,
            PROJECTION=args.projection,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            ATM=args.atm_correction,
            POLE_TIDE=args.pole_tide,
            DEG1=args.geocenter,
            DEG1_FILE=args.geocenter_file,
            MODEL_DEG1=args.interpolate_geocenter,
            SLR_C20=args.slr_c20,
            SLR_21=args.slr_21,
            SLR_22=args.slr_22,
            SLR_C30=args.slr_c30,
            SLR_C40=args.slr_c40,
            SLR_C50=args.slr_c50,
            DATAFORM=args.format,
            MEAN_FILE=args.mean_file,
            MEANFORM=args.mean_format,
            REMOVE_FILES=args.remove_file,
            REMOVE_FORMAT=args.remove_format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
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
