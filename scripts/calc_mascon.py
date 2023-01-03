#!/usr/bin/env python
u"""
calc_mascon.py
Written by Tyler Sutterley (01/2023)

Calculates a time-series of regional mass anomalies through a least-squares
    mascon procedure from GRACE/GRACE-FO time-variable gravity data

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    -O X, --output-directory X: output directory for mascon files
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
    -F X, --format X: input data format for auxiliary files
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
    --mean-file X: GRACE/GRACE-FO mean file to remove from the harmonic data
    --mean-format X: Input data format for GRACE/GRACE-FO mean file
        ascii
        netCDF4
        HDF5
        gfc
    --mask X: Land-sea mask for redistributing mascon mass and land water flux
    --mascon-file X: index file of mascons spherical harmonics
    --mascon-format X: input format for mascon files
        ascii
        netCDF4
        HDF5
    --redistribute-mascons: redistribute mascon mass over the ocean
    --fit-method X: method for fitting sensitivity kernel to harmonics
        1: mass coefficients
        2: geoid coefficients
    --remove-file X: Monthly files to be removed from the GRACE/GRACE-FO data
    --remove-format X: Input data format for files to be removed
        ascii
        netCDF4
        HDF5
        index-ascii
        index-netCDF4
        index-HDF5
    --redistribute-removed: redistribute removed mass fields over the ocean
    --reconstruct-file X: reconstructed mascon time series file to be removed
    --remove-reconstruct: remove reconstructed mascon time series fields
    --log: Output log of files created for each job
    -V, --verbose: Verbose output of processing run
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
    time_series.smooth.py: smoothes a time-series using a Loess-type algorithm
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    time.py: utilities for calculating time operations
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

REFERENCES:
    I Velicogna, T C Sutterley and M R van den Broeke. "Regional acceleration
        in ice mass loss from Greenland and Antarctica using GRACE
        time-variable gravity data". Geophysical Research Letters,
        41(22):8130-8137, 2014. https://doi.org/10.1002/2014GL061052

    T Jacob, J Wahr, W Pfeffer, and S C Swenson "Recent contributions of
        glaciers and ice caps to sea level rise". Nature, 482, 514-518 (2012).
        https://doi.org/10.1038/nature10847

    V M Tiwari, J Wahr, S and Swenson, "Dwindling groundwater resources in
        northern India, from satellite gravity observations",
        Geophysical Research Letters, 36(18), L18401, (2009).
        https://doi.org/10.1029/2009GL039401

    J Wahr, S C Swenson, and I Velicogna, "Accuracy of GRACE mass estimates",
        Geophysical Research Letters, 33(6), L06401, (2006).
        https://doi.org/10.1029/2005GL025305

UPDATE HISTORY:
    Updated 01/2023: refactored time series analysis functions
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 09/2022: add option to replace degree 4 zonal harmonics with SLR
    Updated 04/2022: use wrapper function for reading load Love numbers
        include utf-8 encoding in reads to be windows compliant
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
        option to specify a specific geocenter correction file
    Updated 11/2021: add GSFC low-degree harmonics
    Updated 10/2021: using python logging for handling verbose output
        fix choices for setting input format of the removed files
    Updated 08/2021: reorganize GRACE/GRACE-FO file import
        add option for setting input format of the mascon files
    Updated 07/2021: simplified file imports using wrappers in harmonics
        added path to default land-sea mask for mass redistribution
        remove choices for argparse processing centers
    Updated 06/2021: switch from parameter files to argparse arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 04/2021: include parameters for replacing C21/S21 and C22/S22
        add parser object for removing commented or empty lines
    Updated 02/2021: changed remove index to files with specified formats
    Updated 01/2021: harmonics object output from gen_stokes.py/ocean_stokes.py
    Updated 12/2020: added more love number options and from gfc for mean files
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: use utilities to define path to load love numbers file
    Updated 04/2020: using harmonics class for spherical harmonic operations
        updated load love numbers read function
    Updated 03/2020: switched to destripe_harmonics for filtering harmonics
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: can replace C30 with coefficients from SLR
    Updated 05/2019: use exact date for calculating GIA drift
    Updated 10/2018: verify integers for python3 compatibility
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2018: include a flag if using atmospheric ECMWF "jump" files
        added extrapolation of load love numbers if LMAX is greater than 696
    Updated 09/2017: use a different land-sea mask for calculating ocean_Ylms
        use rcond=-1 in numpy least-squares algorithm
    Updated 02/2017: input spherical harmonic mapping indices as integers
        clean up ocean redistribution read function
    Updated 01/2017: minor clean up of file name definitions
    Updated 06/2016: output destripe strings for sensitivity kernels and error
    Updated 05/2016: using __future__ print function
    Updated 02/2016: direct calculation of number of harmonics n_harm
        added POLE_TIDE and ATM parameters inputs if using command-line entry
        use getopt parameters to set number of PROCESSES to run in parallel,
            whether or not to output a log file, added new help module
    Updated 12/2015: fixed DELTA ascii file output
    Updated 11/2015: create unique log filenames
    Updated 08/2015: changed sys.exit to a raise exception instance
        Added pole tide and GAE/GAF/GAG correction parameters
    Updated 07/2015: added output of the sensitivity kernel Ylms in addition
        to the spatial fields (rather than just the spatial fields)
    Updated 06/2015: will output logs with parameters and output_files
    Updated 06/2015: changed RECONSTRUCT to take the summation of several
        mascons listed by an index rather than a precomputed field
    Updated 05/2015: added parameter MMAX for LMAX != MMAX
    Updated 03/2015: changed remove Ylms from individual (GLDAS,RACMO, etc)
        to indices similar to the GRACE processing programs
        updated comments and header
    Updated 02/2015: added portion to redistribute GLDAS and mascon
        mass uniformly over the ocean
    Updated 01/2015: added error handling for multiprocessing threads
    Updated 11/2014: added INTERVAL parameter for (degree spacing)/2
        input/output file type (ascii, netCDF4, HDF5)
    Updated 10/2014: distributed computing with the multiprocessing module
    Updated 05/2014: added import functions
    Updated 02/2014: updated comments and added os.path.joins for connecting
        directories and files (generalizing code)
        some general updates to the program code
    Updated 09/2013: saving GRACE DELTA file (won't calculate each time)
        added option to remove RACMO data
    Updated 08/2013: general updates to inputting data
    Updated 03/2012: edited to use new gen_stokes time-series option
    Updated 02/2012: Added sensitivity kernels
    Written 02/2012
"""
from __future__ import print_function, division

import sys
import os
import re
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

# PURPOSE: calculate a regional time-series through a least
# squares mascon process
def calc_mascon(base_dir, PROC, DREL, DSET, LMAX, RAD,
    START=None,
    END=None,
    MISSING=None,
    LMIN=None,
    MMAX=None,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    REFERENCE=None,
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
    MASCON_FILE=None,
    MASCON_FORMAT=None,
    REDISTRIBUTE_MASCONS=False,
    FIT_METHOD=0,
    REMOVE_FILES=None,
    REMOVE_FORMAT=None,
    REDISTRIBUTE_REMOVED=False,
    RECONSTRUCT=False,
    RECONSTRUCT_FILE=None,
    LANDMASK=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # recursively create output Directory if not currently existing
    if (not os.access(OUTPUT_DIRECTORY, os.F_OK)):
        os.makedirs(OUTPUT_DIRECTORY, mode=MODE, exist_ok=True)

    # list object of output files for file logs (full path)
    output_files = []

    # file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # file parser for reading index files
    # removes commented lines (can comment out files in the index)
    # removes empty lines (if there are extra empty lines)
    parser = re.compile(r'^(?!\#|\%|$)', re.VERBOSE)

    # read arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # Earth Parameters
    factors = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll)
    # Average Density of the Earth [g/cm^3]
    rho_e = factors.rho_e
    # Average Radius of the Earth [cm]
    rad_e = factors.rad_e

    # for datasets not GSM: will add a label for the dataset
    dset_str = '' if (DSET == 'GSM') else f'_{DSET}'
    # atmospheric ECMWF "jump" flag (if ATM)
    atm_str = '_wATM' if ATM else ''
    # output string for both LMAX==MMAX and LMAX != MMAX cases
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # Read Ocean function and convert to Ylms for redistribution
    if (REDISTRIBUTE_MASCONS | REDISTRIBUTE_REMOVED):
        # read Land-Sea Mask and convert to spherical harmonics
        ocean_Ylms = gravtk.ocean_stokes(LANDMASK, LMAX, MMAX=MMAX,
            LOVE=(hl,kl,ll))
        ocean_str = '_OCN'
    else:
        # not distributing uniformly over ocean
        ocean_str = ''

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
    # use a mean file for the static field to remove
    if MEAN_FILE:
        # read data form for input mean file (ascii, netCDF4, HDF5, gfc)
        mean_Ylms = gravtk.harmonics().from_file(MEAN_FILE,
            format=MEANFORM, date=False)
        # remove the input mean
        GRACE_Ylms.subtract(mean_Ylms)
    else:
        GRACE_Ylms.mean(apply=True)
    # filter GRACE/GRACE-FO coefficients
    if DESTRIPE:
        # destriping GRACE/GRACE-FO coefficients
        ds_str = '_FL'
        GRACE_Ylms = GRACE_Ylms.destripe()
    else:
        # using standard GRACE/GRACE-FO harmonics
        ds_str = ''
    # full path to directory for specific GRACE/GRACE-FO product
    GRACE_Ylms.directory = Ylms['directory']
    # date information of GRACE/GRACE-FO coefficients
    n_files = len(GRACE_Ylms.time)

    # input GIA spherical harmonic datafiles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA, mmax=MMAX)
    gia_str = f'_{GIA_Ylms_rate.title}' if GIA else ''
    # monthly GIA calculated by gia_rate*time elapsed
    # finding change in GIA each month
    GIA_Ylms = GIA_Ylms_rate.drift(GRACE_Ylms.time, epoch=2003.3)
    GIA_Ylms.month[:] = np.copy(GRACE_Ylms.month)

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

    # input reconstructed spherical harmonic datafiles
    construct_Ylms = GRACE_Ylms.zeros_like()
    construct_Ylms.time[:] = np.copy(GRACE_Ylms.time)
    construct_Ylms.month[:] = np.copy(GRACE_Ylms.month)
    if RECONSTRUCT:
        # input index for reconstructed spherical harmonic datafiles
        with open(RECONSTRUCT_FILE, mode='r', encoding='utf8') as f:
            file_list = [l for l in f.read().splitlines() if parser.match(l)]
        # for each valid file in the index (iterate over mascons)
        for reconstruct_file in file_list:
            # read reconstructed spherical harmonics
            Ylms = gravtk.harmonics().from_file(reconstruct_file,
                format=DATAFORM)
            # truncate clm and slm matrices to LMAX/MMAX
            # add harmonics object to total
            construct_Ylms.add(Ylms.truncate(lmax=LMAX, mmax=MMAX))
        # filter reconstructed coefficients
        if DESTRIPE:
            construct_Ylms = construct_Ylms.destripe()
        # set flag for removing reconstructed coefficients
        construct_str = '_LEAKAGE'
    else:
        # set flag for not removing the reconstructed coefficients
        construct_str = ''

    # input mascon spherical harmonic datafiles
    with open(MASCON_FILE, mode='r', encoding='utf8') as f:
        mascon_files = [l for l in f.read().splitlines() if parser.match(l)]
    # number of mascons
    n_mas = len(mascon_files)
    # spatial area of the mascon
    total_area = np.zeros((n_mas))
    # name of each mascon
    mascon_name = []
    # for each valid file in the index (iterate over mascons)
    mascon_list = []
    for k,fi in enumerate(mascon_files):
        # read mascon spherical harmonics
        Ylms = gravtk.harmonics().from_file(os.path.expanduser(fi),
            format=MASCON_FORMAT, date=False)
        # Calculating the total mass of each mascon (1 cmwe uniform)
        total_area[k] = 4.0*np.pi*(rad_e**3)*rho_e*Ylms.clm[0,0]/3.0
        # distribute mascon mass uniformly over the ocean
        if REDISTRIBUTE_MASCONS:
            # calculate ratio between total mascon mass and
            # a uniformly distributed cm of water over the ocean
            ratio = Ylms.clm[0,0]/ocean_Ylms.clm[0,0]
            # for each spherical harmonic
            for m in range(0,MMAX+1):# MMAX+1 to include MMAX
                for l in range(m,LMAX+1):# LMAX+1 to include LMAX
                    # remove ratio*ocean Ylms from mascon Ylms
                    # note: x -= y is equivalent to x = x - y
                    Ylms.clm[l,m] -= ratio*ocean_Ylms.clm[l,m]
                    Ylms.slm[l,m] -= ratio*ocean_Ylms.slm[l,m]
        # truncate mascon spherical harmonics to d/o LMAX/MMAX and add to list
        mascon_list.append(Ylms.truncate(lmax=LMAX, mmax=MMAX))
        # mascon base is the file without directory or suffix
        mascon_base = os.path.basename(mascon_files[k])
        mascon_base = os.path.splitext(mascon_base)[0]
        # if lower case, will capitalize
        mascon_base = mascon_base.upper()
        # if mascon name contains degree and order info, remove
        mascon_name.append(mascon_base.replace(f'_L{LMAX:d}', ''))
    # create single harmonics object from list
    mascon_Ylms = gravtk.harmonics().from_list(mascon_list, date=False)
    # clear mascon list variable
    del mascon_list

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
        delta_Ylms = gravtk.harmonics().from_file(DELTA_FILE,
            format=DATAFORM)
        # truncate GRACE/GRACE-FO delta clm and slm to d/o LMAX/MMAX
        delta_Ylms = delta_Ylms.truncate(lmax=LMAX, mmax=MMAX)
        tsmth = np.squeeze(delta_Ylms.time)
        nsmth = np.int64(delta_Ylms.month)

    # Calculating the number of cos and sin harmonics between LMIN and LMAX
    # taking into account MMAX (if MMAX == LMAX then LMAX-MMAX=0)
    n_harm=np.int64(LMAX**2 - LMIN**2 + 2*LMAX + 1 - (LMAX-MMAX)**2 - (LMAX-MMAX))

    # Initialing harmonics for least squares fitting
    # mascon kernel
    M_lm = np.zeros((n_harm,n_mas))
    # mascon kernel converted to output unit
    MA_lm = np.zeros((n_harm,n_mas))
    # corrected clm and slm
    Y_lm = np.zeros((n_harm,n_files))
    # sensitivity kernel
    A_lm = np.zeros((n_harm,n_mas))
    # Satellite error harmonics
    delta_lm = np.zeros((n_harm))
    # Initializing output Mascon time-series
    mascon = np.zeros((n_mas,n_files))
    # Mascon satellite error component
    M_delta = np.zeros((n_mas))
    # Initializing conversion factors
    # factor for converting to coefficients of mass
    fact = np.zeros((n_harm))
    # smoothing factor
    wt_lm = np.zeros((n_harm))

    # ii is a counter variable for building the mascon column array
    ii = 0
    # Creating column array of clm/slm coefficients
    # Order is [C00...C6060,S11...S6060]
    # Calculating factor to convert geoid spherical harmonic coefficients
    # to coefficients of mass (Wahr, 1998)
    coeff = rho_e*rad_e/3.0
    # Switching between Cosine and Sine Stokes
    for cs,csharm in enumerate(['clm','slm']):
        # copy cosine and sin harmonics
        mascon_harm = getattr(mascon_Ylms, csharm)
        grace_harm = getattr(GRACE_Ylms, csharm)
        GIA_harm = getattr(GIA_Ylms, csharm)
        remove_harm = getattr(remove_Ylms, csharm)
        construct_harm = getattr(construct_Ylms, csharm)
        delta_harm = getattr(delta_Ylms, csharm)
        # for each spherical harmonic degree
        # +1 to include LMAX
        for l in range(LMIN,LMAX+1):
            # for each spherical harmonic order
            # Sine Stokes for (m=0) = 0
            mm = np.min([MMAX,l])
            # +1 to include l or MMAX (whichever is smaller)
            for m in range(cs,mm+1):
                # Mascon Spherical Harmonics
                M_lm[ii,:] = np.copy(mascon_harm[l,m,:])
                # GRACE Spherical Harmonics
                # Correcting GRACE Harmonics for GIA and Removed Terms
                Y_lm[ii,:] = grace_harm[l,m,:] - GIA_harm[l,m,:] - \
                    remove_harm[l,m,:] - construct_harm[l,m,:]
                # GRACE delta spherical harmonics
                delta_lm[ii] = np.copy(delta_harm[l,m])
                # degree dependent factor to convert to mass
                fact[ii] = (2.0*l + 1.0)/(1.0 + kl[l])
                # degree dependent smoothing
                wt_lm[ii] = np.copy(wt[l])
                # add 1 to counter
                ii += 1

    # Converting mascon coefficients to fit method
    if (FIT_METHOD == 1):
        # Fitting Sensitivity Kernel as mass coefficients
        # converting M_lm to mass coefficients of the kernel
        for i in range(n_harm):
            MA_lm[i,:] = M_lm[i,:]*wt_lm[i]*fact[i]
        fit_factor = wt_lm*fact
    else:
        # Fitting Sensitivity Kernel as geoid coefficients
        for i in range(n_harm):
            MA_lm[:,:] = M_lm[i,:]*wt_lm[i]
        fit_factor = wt_lm*np.ones((n_harm))

    # Fitting the sensitivity kernel from the input kernel
    for i in range(n_harm):
        # setting kern_i equal to 1 for d/o
        kern_i = np.zeros((n_harm))
        # converting to mass coefficients if specified
        kern_i[i] = 1.0*fit_factor[i]
        # spherical harmonics solution for the
        # mascon sensitivity kernels
        # Least Squares Solutions: Inv(X'.X).(X'.Y)
        kern_lm = np.linalg.lstsq(MA_lm,kern_i,rcond=-1)[0]
        for k in range(n_mas):
            A_lm[i,k] = kern_lm[k]*total_area[k]

    # for each mascon
    for k in range(n_mas):
        # Multiply the Satellite error (noise of a smoothed time-series
        # with annual and semi-annual components) by the sensitivity kernel
        # Converting to Gigatonnes
        M_delta[k] = np.sqrt(np.sum((delta_lm*A_lm[:,k])**2))/1e15

        # output filename format (for both LMAX==MMAX and LMAX != MMAX cases):
        # mascon name, GRACE dataset, GIA model, LMAX, (MMAX,)
        # Gaussian smoothing, filter flag, remove reconstructed fields flag
        # output GRACE error file
        fargs = (mascon_name[k], dset_str, gia_str.upper(), atm_str, ocean_str,
            LMAX, order_str, gw_str, ds_str, construct_str)
        file_format = '{0}{1}{2}{3}{4}_L{5:d}{6}{7}{8}{9}.txt'
        output_file = os.path.join(OUTPUT_DIRECTORY,file_format.format(*fargs))

        # Output mascon datafiles
        # Will output each mascon time series
        # month, date, mascon mass [Gt], satellite error [Gt], mascon area [km^2]
        # open output mascon time-series file
        fid = open(output_file, mode='w', encoding='utf8')
        # for each date
        formatting_string = '{0:03d} {1:12.4f} {2:16.10f} {3:16.10f} {4:16.5f}'
        for t,mon in enumerate(GRACE_Ylms.month):
            # Summing over all spherical harmonics for mascon k, and time t
            # multiplies by the degree dependent factor to convert
            # the harmonics into mass coefficients
            # Converting mascon mass time-series from g to gigatonnes
            mascon[k,t] = np.sum(A_lm[:,k]*Y_lm[:,t])/1e15
            # output to file
            args=(mon,GRACE_Ylms.time[t],mascon[k,t],M_delta[k],total_area[k]/1e10)
            print(formatting_string.format(*args), file=fid)
        # close the output file
        fid.close()
        # change the permissions mode
        os.chmod(output_file, MODE)
        # add output files to list object
        output_files.append(output_file)

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the GRACE mascon analysis
def output_log_file(input_arguments, output_files):
    # format: calc_mascon_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_mascon_run_{0}_PID-{1:d}.log'.format(*args)
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

# PURPOSE: print a error file log for the GRACE mascon analysis
def output_error_log_file(input_arguments):
    # format: calc_mascon_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'calc_mascon_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Calculates a time-series of regional mass anomalies
            through a least-squares mascon procedure from GRACE/GRACE-FO
            time-variable gravity data
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
        help='Output directory for mascon files')
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
        help='Input data format for auxiliary files')
    # mean file to remove
    parser.add_argument('--mean-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GRACE/GRACE-FO mean file to remove from the harmonic data')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--mean-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5','gfc'],
        help='Input data format for GRACE/GRACE-FO mean file')
    # mascon index file and parameters
    parser.add_argument('--mascon-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Index file of mascons spherical harmonics')
    parser.add_argument('--mascon-format',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input data format for mascon files')
    parser.add_argument('--redistribute-mascons',
        default=False, action='store_true',
        help='Redistribute mascon mass over the ocean')
    # 1: mass coefficients
    # 2: geoid coefficients
    parser.add_argument('--fit-method',
        type=int, default=1, choices=(1,2),
        help='Method for fitting sensitivity kernel to harmonics')
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
    # mascon reconstruct parameters
    parser.add_argument('--remove-reconstruct',
        default=False, action='store_true',
        help='Remove reconstructed mascon time series fields')
    parser.add_argument('--reconstruct-file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Reconstructed mascon time series file to be removed')
    # land-sea mask for redistributing mascon mass and land water flux
    lsmask = gravtk.utilities.get_data_path(['data','landsea_hd.nc'])
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), default=lsmask,
        help='Land-sea mask for redistributing mascon mass and land water flux')
    # Output log file for each job in forms
    # calc_mascon_run_2002-04-01_PID-00000.log
    # calc_mascon_failed_run_2002-04-01_PID-00000.log
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
        # run calc_mascon algorithm with parameters
        output_files = calc_mascon(
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
            DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
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
            MASCON_FILE=args.mascon_file,
            MASCON_FORMAT=args.mascon_format,
            REDISTRIBUTE_MASCONS=args.redistribute_mascons,
            FIT_METHOD=args.fit_method,
            REMOVE_FILES=args.remove_file,
            REMOVE_FORMAT=args.remove_format,
            REDISTRIBUTE_REMOVED=args.redistribute_removed,
            RECONSTRUCT=args.remove_reconstruct,
            RECONSTRUCT_FILE=args.reconstruct_file,
            LANDMASK=args.mask,
            OUTPUT_DIRECTORY=args.output_directory,
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
