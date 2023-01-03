#!/usr/bin/env python
u"""
scale_grace_maps.py
Written by Tyler Sutterley (01/2023)

Reads in GRACE/GRACE-FO spherical harmonic coefficients and exports
    monthly scaled spatial fields, estimated scaling errors,
    and estimated scaled delta errors

Will correct with the specified GIA model group, destripe/smooth/process,
    and export the data in centimeters water equivalent (cm w.e.)

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for spatial files
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
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: input/output data format
        ascii
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
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
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
    --scale-file X: scaling factor file
    --error-file X: scaling factor error file
    --power-file X: scaling factor power file
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
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    grace_input_months.py: Reads GRACE/GRACE-FO files for a specified date range
        Includes degree 1 values (if specified)
        Replaces low-degree harmonics with SLR values (if specified)
    read_GIA_model.py: reads spherical harmonics for glacial isostatic adjustment
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    ocean_stokes.py: reads a land-sea mask and converts to spherical harmonics
    gen_stokes.py: converts a spatial field into spherical harmonic coefficients
    geocenter.py: converts between spherical harmonics and geocenter variations
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    time_series.smooth.py: smoothes a time-series using a Loess-type algorithm
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    C-W Hsu and I Velicogna, "Detection of Sea Level Fingerprints derived from
        GRACE gravity data", Geophysical Research Letters, 44(17), (2017).
        https://doi.org/10.1002/2017GL074070

    F W Landerer and S C Swenson, "Accuracy of scaled GRACE terrestrial water
        storage estimates", Water Resources Research, 48(4), W04531, (2012).
        https://doi.org/10.1029/2011WR011453

    J Wahr, S C Swenson, and I Velicogna, "Accuracy of GRACE mass estimates",
        Geophysical Research Letters, 33(6), L06401, (2006).
        https://doi.org/10.1029/2005GL025305

UPDATE HISTORY:
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
        option to specify a specific geocenter correction file
        fix default file prefix to include center and release information
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 10/2021: using python logging for handling verbose output
        add more choices for setting input format of the removed files
    Updated 07/2021: switch from parameter files to argparse arguments
        simplified file imports using wrappers in harmonics
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
    Updated 02/2021: changed remove index to files with specified formats
    Updated 02/2021: for public release
"""
from __future__ import print_function, division

import sys
import os
import copy
import time
import logging
import argparse
import numpy as np
import traceback
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: import GRACE/GRACE-FO files for a given months range
# Calculates monthly scaled spatial maps from GRACE/GRACE-FO
# spherical harmonic coefficients
def scale_grace_maps(base_dir, PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    LMIN=None,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DESTRIPE=False,
    DDEG=None,
    INTERVAL=None,
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
    SCALE_FILE=None,
    ERROR_FILE=None,
    POWER_FILE=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    FILE_PREFIX=None,
    VERBOSE=0,
    MODE=0o775):

    # recursively create output Directory if not currently existing
    if not os.access(OUTPUT_DIRECTORY, os.F_OK):
        os.makedirs(OUTPUT_DIRECTORY, mode=MODE, exist_ok=True)

    # list object of output files for file logs (full path)
    output_files = []

    # file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # output file format
    file_format = '{0}{1}{2}_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'

    # read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # output spatial units
    unit_str = 'cmwe'
    unit_name = 'Equivalent_Water_Thickness'
    # invalid value
    fill_value = -9999.0

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # Read Ocean function and convert to Ylms for redistribution
    if REDISTRIBUTE_REMOVED:
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX, MMAX=MMAX,
            LOVE=(hl,kl,ll))

    # Grid spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Grid dimensions
    if (INTERVAL == 1):# (0:360, 90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):# degree spacing/2
        nlon = np.int64((360.0/dlon))
        nlat = np.int64((180.0/dlat))

    # read data for input scale files (ascii, netCDF4, HDF5)
    if (DATAFORM == 'ascii'):
        kfactor = gravtk.spatial(spacing=[dlon,dlat], nlat=nlat,
            nlon=nlon).from_ascii(SCALE_FILE, date=False)
        k_error = gravtk.spatial(spacing=[dlon,dlat], nlat=nlat,
            nlon=nlon).from_ascii(ERROR_FILE, date=False)
        k_power = gravtk.spatial(spacing=[dlon,dlat], nlat=nlat,
            nlon=nlon).from_ascii(POWER_FILE, date=False)
    elif (DATAFORM == 'netCDF4'):
        kfactor = gravtk.spatial().from_netCDF4(SCALE_FILE, date=False)
        k_error = gravtk.spatial().from_netCDF4(ERROR_FILE, date=False)
        k_power = gravtk.spatial().from_netCDF4(POWER_FILE, date=False)
    elif (DATAFORM == 'HDF5'):
        kfactor = gravtk.spatial().from_HDF5(SCALE_FILE, date=False)
        k_error = gravtk.spatial().from_HDF5(ERROR_FILE, date=False)
        k_power = gravtk.spatial().from_HDF5(POWER_FILE, date=False)
    # input data shape
    nlat,nlon = kfactor.shape

    # input GRACE/GRACE-FO spherical harmonic datafiles for date range
    # replacing low-degree harmonics with SLR values if specified
    # include degree 1 (geocenter) harmonics if specified
    # correcting for Pole-Tide and Atmospheric Jumps if specified
    Ylms = gravtk.grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        START, END, MISSING, SLR_C20, DEG1, MMAX=MMAX, SLR_21=SLR_21,
        SLR_22=SLR_22, SLR_C30=SLR_C30, SLR_C40=SLR_C40, SLR_C50=SLR_C50,
        DEG1_FILE=DEG1_FILE, MODEL_DEG1=MODEL_DEG1, ATM=ATM,
        POLE_TIDE=POLE_TIDE)
    # create harmonics object from GRACE/GRACE-FO data
    GRACE_Ylms = gravtk.harmonics().from_dict(Ylms)
    GRACE_Ylms.directory = Ylms['directory']
    # use a mean file for the static field to remove
    if MEAN_FILE:
        # read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=MEANFORM, date=False)
        # remove the input mean
        GRACE_Ylms.subtract(mean_Ylms)
    else:
        GRACE_Ylms.mean(apply=True)
    # date information of GRACE/GRACE-FO coefficients
    nfiles = len(GRACE_Ylms.time)

    # filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        # destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
    else:
        # using standard GRACE/GRACE-FO harmonics
        ds_str = ''

    # input GIA spherical harmonic datafiles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
    gia_str = f'_{GIA_Ylms_rate.title}' if GIA else ''
    # monthly GIA calculated by gia_rate*time elapsed
    # finding change in GIA each month
    GIA_Ylms = GIA_Ylms_rate.drift(GRACE_Ylms.time, epoch=2003.3)
    GIA_Ylms.month[:] = np.copy(GRACE_Ylms.month)

    # default file prefix
    if not FILE_PREFIX:
        fargs = (PROC,DREL,DSET,Ylms['title'],gia_str)
        FILE_PREFIX = '{0}_{1}_{2}{3}{4}_'.format(*fargs)

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
            elif REMOVEFORM in ('index-ascii','index-netCDF4','index-HDF5'):
                # read from index file
                _,removeform = REMOVEFORM.split('-')
                # index containing files in data format
                Ylms = gravtk.harmonics().from_index(REMOVE_FILE,
                    format=removeform)
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

    # calculating GRACE/GRACE-FO error (Wahr et al. 2006)
    # output GRACE error file (for both LMAX==MMAX and LMAX != MMAX cases)
    fargs = (PROC,DREL,DSET,LMAX,order_str,ds_str,atm_str,GRACE_Ylms.month[0],
        GRACE_Ylms.month[-1], suffix[DATAFORM])
    delta_format = '{0}_{1}_{2}_DELTA_CLM_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    DELTA_FILE = os.path.join(GRACE_Ylms.directory,delta_format.format(*fargs))
    # check full path of the GRACE directory for delta file
    # if file was previously calculated: will read file
    # else: will calculate the GRACE/GRACE-FO error
    if not os.access(DELTA_FILE, os.F_OK):
        # add output delta file to list object
        output_files.append(DELTA_FILE)

        # Delta coefficients of GRACE time series (Error components)
        delta_Ylms = gravtk.harmonics(lmax=LMAX,mmax=MMAX)
        delta_Ylms.clm = np.zeros((LMAX+1,MMAX+1))
        delta_Ylms.slm = np.zeros((LMAX+1,MMAX+1))
        # Smoothing Half-Width (CNES is a 10-day solution)
        # All other solutions are monthly solutions (HFWTH for annual = 6)
        if ((PROC == 'CNES') and (DREL in ('RL01','RL02'))):
            HFWTH = 19
        else:
            HFWTH = 6
        # Equal to the noise of the smoothed time-series
        # for each spherical harmonic order
        for m in range(0,MMAX+1):# MMAX+1 to include MMAX
            # for each spherical harmonic degree
            for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                # Delta coefficients of GRACE time series
                for cs,csharm in enumerate(['clm','slm']):
                    # calculate GRACE Error (Noise of smoothed time-series)
                    # With Annual and Semi-Annual Terms
                    val1 = getattr(GRACE_Ylms, csharm)
                    smth = gravtk.time_series.smooth(GRACE_Ylms.time,
                        val1[l,m,:], HFWTH=HFWTH)
                    # number of smoothed points
                    nsmth = len(smth['data'])
                    tsmth = np.mean(smth['time'])
                    # GRACE/GRACE-FO delta Ylms
                    # variance of data-(smoothed+annual+semi)
                    val2 = getattr(delta_Ylms, csharm)
                    val2[l,m] = np.sqrt(np.sum(smth['noise']**2)/nsmth)

        # attributes for output files
        attributes = {}
        attributes['title'] = 'GRACE/GRACE-FO Spherical Harmonic Errors'
        attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
        # save GRACE/GRACE-FO delta harmonics to file
        delta_Ylms.time = np.copy(tsmth)
        delta_Ylms.month = np.int64(nsmth)
        delta_Ylms.to_file(DELTA_FILE, format=DATAFORM, **attributes)
        # set the permissions mode of the output harmonics file
        os.chmod(DELTA_FILE, MODE)
        # append delta harmonics file to output files list
        output_files.append(DELTA_FILE)
    else:
        # read GRACE/GRACE-FO delta harmonics from file
        delta_Ylms = gravtk.harmonics().from_file(DELTA_FILE, format=DATAFORM)
        # copy time and number of smoothed fields
        tsmth = np.squeeze(delta_Ylms.time)
        nsmth = np.int64(delta_Ylms.month)

    # Output spatial data object
    grid = gravtk.spatial()
    grid.lon = np.copy(kfactor.lon)
    grid.lat = np.copy(kfactor.lat)
    grid.time = np.zeros((nfiles))
    grid.month = np.zeros((nfiles),dtype=np.int64)
    grid.data = np.zeros((nlat,nlon,nfiles))
    grid.mask = np.zeros((nlat,nlon,nfiles),dtype=bool)

    # Computing plms for converting to spatial domain
    phi = grid.lon[np.newaxis,:]*np.pi/180.0
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))
    # square of legendre polynomials truncated to order MMAX
    mm = np.arange(0,MMAX+1)
    PLM2 = PLM[:,mm,:]**2

    # dfactor is the degree dependent coefficients
    # for converting to centimeters water equivalent (cmwe)
    dfactor = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll).cmwe

    # converting harmonics to truncated, smoothed coefficients in units
    # combining harmonics to calculate output spatial fields
    for i,gm in enumerate(GRACE_Ylms.month):
        # GRACE/GRACE-FO harmonics for time t
        Ylms = GRACE_Ylms.index(i)
        # Remove GIA rate for time
        Ylms.subtract(GIA_Ylms.index(i))
        # Remove monthly files to be removed
        Ylms.subtract(remove_Ylms.index(i))
        # smooth harmonics and convert to output units
        Ylms.convolve(dfactor*wt)
        # convert spherical harmonics to output spatial grid
        grid.data[:,:,i] = gravtk.harmonic_summation(Ylms.clm, Ylms.slm,
            grid.lon, grid.lat, LMAX=LMAX, MMAX=MMAX, PLM=PLM).T
        # copy time variables for month
        grid.time[i] = np.copy(Ylms.time)
        grid.month[i] = np.copy(Ylms.month)
    # update spacing and dimensions
    grid.update_spacing()
    grid.update_extents()
    grid.update_dimensions()

    # scale output data with kfactor
    grid = grid.scale(kfactor.data)
    grid.replace_invalid(fill_value, mask=kfactor.mask)

    # output monthly files to ascii, netCDF4 or HDF5
    fargs = (FILE_PREFIX, '', unit_str, LMAX, order_str, gw_str,
        ds_str, grid.month[0], grid.month[-1], suffix[DATAFORM])
    FILE = os.path.join(OUTPUT_DIRECTORY,file_format.format(*fargs))
    # attributes for output files
    attributes = {}
    attributes['units'] = copy.copy(unit_str)
    attributes['longname'] = copy.copy(unit_name)
    attributes['title'] = 'GRACE/GRACE-FO Spatial Data'
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        grid.to_ascii(FILE, date=True, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4
        grid.to_netCDF4(FILE, date=True, verbose=VERBOSE, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5
        grid.to_HDF5(FILE, date=True, verbose=VERBOSE, **attributes)
    # set the permissions mode of the output files
    os.chmod(FILE, MODE)
    # add file to list
    output_files.append(FILE)

    # calculate power of scaled GRACE/GRACE-FO data
    scaled_power = grid.sum(power=2.0).power(0.5)
    # calculate residual leakage errors
    # scaled by ratio of GRACE and synthetic power
    ratio = scaled_power.scale(k_power.power(-1).data)
    error = k_error.scale(ratio.data)

    # output monthly error files to ascii, netCDF4 or HDF5
    fargs = (FILE_PREFIX, 'ERROR_', unit_str, LMAX, order_str, gw_str,
        ds_str, grid.month[0], grid.month[-1], suffix[DATAFORM])
    FILE = os.path.join(OUTPUT_DIRECTORY,file_format.format(*fargs))
    # attributes for output files
    attributes['title'] = 'GRACE/GRACE-FO Scaling Error'
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        error.to_ascii(FILE, date=False, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4
        error.to_netCDF4(FILE, date=False, verbose=VERBOSE, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5
        error.to_HDF5(FILE, date=False, verbose=VERBOSE, **attributes)
    # set the permissions mode of the output files
    os.chmod(FILE, MODE)
    # add file to list
    output_files.append(FILE)

    # Output spatial data object
    delta = gravtk.spatial()
    delta.lon = np.copy(kfactor.lon)
    delta.lat = np.copy(kfactor.lat)
    delta.time = np.copy(tsmth)
    delta.month = np.copy(nsmth)
    delta.data = np.zeros((nlat,nlon))
    delta.mask = np.zeros((nlat,nlon),dtype=bool)
    # calculate scaled spatial error
    # Calculating cos(m*phi)^2 and sin(m*phi)^2
    m = delta_Ylms.m[:,np.newaxis]
    ccos = np.cos(np.dot(m,phi))**2
    ssin = np.sin(np.dot(m,phi))**2

    # truncate delta harmonics to spherical harmonic range
    Ylms = delta_Ylms.truncate(LMAX,lmin=LMIN,mmax=MMAX)
    # convolve delta harmonics with degree dependent factors
    # smooth harmonics and convert to output units
    Ylms = Ylms.convolve(dfactor*wt).power(2.0).scale(1.0/nsmth)
    # Calculate fourier coefficients
    d_cos = np.zeros((MMAX+1,nlat))# [m,th]
    d_sin = np.zeros((MMAX+1,nlat))# [m,th]
    # Calculating delta spatial values
    for k in range(0,nlat):
        # summation over all spherical harmonic degrees
        d_cos[:,k] = np.sum(PLM2[:,:,k]*Ylms.clm, axis=0)
        d_sin[:,k] = np.sum(PLM2[:,:,k]*Ylms.slm, axis=0)
    # Multiplying by c/s(phi#m) to get spatial error map
    delta.data[:] = np.sqrt(np.dot(ccos.T,d_cos) + np.dot(ssin.T,d_sin)).T
    # update spacing and dimensions
    delta.update_spacing()
    delta.update_extents()
    delta.update_dimensions()

    # scale output harmonic errors with kfactor
    delta = delta.scale(kfactor.data)
    delta.replace_invalid(fill_value, mask=kfactor.mask)

    # output monthly files to ascii, netCDF4 or HDF5
    fargs = (FILE_PREFIX, 'DELTA_', unit_str, LMAX, order_str, gw_str,
        ds_str, grid.month[0], grid.month[-1], suffix[DATAFORM])
    FILE = os.path.join(OUTPUT_DIRECTORY,file_format.format(*fargs))
    # attributes for output files
    attributes['title'] = 'GRACE/GRACE-FO Spatial Error'
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        delta.to_ascii(FILE, date=True, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4
        delta.to_netCDF4(FILE, date=True, verbose=VERBOSE, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5
        delta.to_HDF5(FILE, date=True, verbose=VERBOSE, **attributes)
    # set the permissions mode of the output files
    os.chmod(FILE, MODE)
    # add file to list
    output_files.append(FILE)

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the GRACE/GRACE-FO analysis
def output_log_file(input_arguments, output_files):
    # format: scale_GRACE_maps_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'scale_GRACE_maps_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
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

# PURPOSE: print a error file log for the GRACE/GRACE-FO analysis
def output_error_log_file(input_arguments):
    # format: scale_GRACE_maps_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'scale_GRACE_maps_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = os.path.expanduser(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(os.path.join(DIRECTORY,LOGFILE))
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
        description="""Calculates scaled spatial maps from
            GRACE/GRACE-FO spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    parser.add_argument('--output-directory','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Output directory for spatial files')
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
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
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
    # output grid parameters
    parser.add_argument('--spacing',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval',
        type=int, default=2, choices=[1,2],
        help=('Output grid interval (1: global, 2: centered global)'))
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
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
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
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
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
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    # monthly files to be removed from the GRACE/GRACE-FO data
    parser.add_argument('--remove-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
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
    # scaling factor files
    parser.add_argument('--scale-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        required=True, help='Scaling factor file')
    parser.add_argument('--error-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        required=True, help='Scaling factor error file')
    parser.add_argument('--power-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        required=True, help='Scaling factor power file')
    # land-sea mask for redistributing fluxes
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing land water flux')
    # Output log file for each job in forms
    # scale_GRACE_maps_run_2002-04-01_PID-00000.log
    # scale_GRACE_maps_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run scale_grace_maps algorithm with parameters
        output_files = scale_grace_maps(
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
            DDEG=args.spacing,
            INTERVAL=args.interval,
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
            SCALE_FILE=args.scale_file,
            ERROR_FILE=args.error_file,
            POWER_FILE=args.power_file,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
            FILE_PREFIX=args.file_prefix,
            VERBOSE=args.verbose,
            MODE=args.mode)
    except Exception as e:
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
