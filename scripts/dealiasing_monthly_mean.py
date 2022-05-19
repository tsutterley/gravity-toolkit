#!/usr/bin/env python
u"""
dealiasing_monthly_mean.py
Written by Tyler Sutterley (04/2022)

Reads GRACE/GRACE-FO AOD1B datafiles for a specific product and outputs monthly
    the mean for a specific GRACE/GRACE-FO processing center and data release
    GAA: atmospheric loading from ECMWF
    GAB: oceanic loading from OMCT/MPIOM
    GAC: global atmospheric and oceanic loading
    GAD: ocean bottom pressure from OMCT/MPIOM

Creates a file for each month (such as making GAA and GAB files for CSR)

CALLING SEQUENCE:
    python dealiasing_monthly_mean.py --center CSR --release RL06 --product GAA

COMMAND LINE OPTIONS:
    -D X, --directory X: Working Data Directory
    -c X, --center X: GRACE/GRACE-FO Processing Center
    -r X, --release X: GRACE/GRACE-FO Data Release (RL05 or RL06)
    -p X, --product X: GRACE/GRACE-FO dealiasing product (GAA, GAB, GAC, GAD)
    -l X, --lmax X: Maximum spherical harmonic degree and order for output
    -F X, --format X: Output data format
        ascii
        netCDF4
        HDF5
        SHM
    -C, --clobber: Overwrite existing data
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format
        http://www.h5py.org/

PROGRAM DEPENDENCIES:
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 04/2022: use argparse descriptions within documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use default argument files to define options
        added option to output in spherical harmonic model (SHM) format
        remove choices for argparse processing centers
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: using utilities from time module
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 04/2020: using harmonics class for operations and outputting to file
        reduce output date file to only months with AOD data
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2019: using python3 compatible regular expression patterns
    Updated 10/2018: using future division for python3 Compatibility
    Updated 08/2018: using full release string (RL05 instead of 5)
    Updated 03/2018: copy date file from input GSM directory to output directory
    Written 03/2018
"""
from __future__ import print_function, division

import sys
import os
import re
import gzip
import time
import logging
import tarfile
import argparse
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.utilities as utilities
from gravity_toolkit.harmonics import harmonics

#-- PURPOSE: calculate the Julian day from the year and the day of the year
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def calc_julian_day(YEAR, DAY_OF_YEAR):
    JD = 367.0*YEAR - np.floor(7.0*(YEAR + np.floor(10.0/12.0))/4.0) - \
        np.floor(3.0*(np.floor((YEAR + 8.0/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0/9.0) + np.float64(DAY_OF_YEAR) + 1721028.5
    return JD

#-- PURPOSE: reads the AOD1B data and outputs a monthly mean
def dealiasing_monthly_mean(base_dir, PROC=None, DREL=None, DSET=None,
    LMAX=None, DATAFORM=None, CLOBBER=False, MODE=0o775):

    #-- output data suffix
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    #-- aod1b data products
    aod1b_products = dict(GAA='atm',GAB='ocn',GAC='glo',GAD='oba')
    #-- compile regular expressions operator for the clm/slm headers
    #-- for the specific AOD1b product
    hx = re.compile(r'^DATA.*SET.*{0}'.format(aod1b_products[DSET]),re.VERBOSE)
    #-- compile regular expression operator to find numerical instances
    #-- will extract the data from the file
    regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
    rx = re.compile(regex_pattern, re.VERBOSE)

    #-- set number of hours in a file
    #-- set the ocean model for a given release
    if DREL in ('RL01','RL02','RL03','RL04','RL05'):
        #-- for 00, 06, 12 and 18
        n_time = 4
        ATMOSPHERE = 'ECMWF'
        OCEAN_MODEL = 'OMCT'
        default_center = 'EIGEN'
        default_lmax = 100
    elif DREL in ('RL06',):
        #-- for 00, 03, 06, 09, 12, 15, 18 and 21
        n_time = 8
        ATMOSPHERE = 'ECMWF'
        OCEAN_MODEL = 'MPIOM'
        default_center = 'GFZOP'
        default_lmax = 180
    else:
        raise ValueError('Invalid data release')
    #-- Maximum spherical harmonic degree (LMAX)
    LMAX = default_lmax if not LMAX else LMAX
    #-- Calculating the number of cos and sin harmonics up to d/o of file
    n_harm = (default_lmax**2 + 3*default_lmax)//2 + 1

    #-- AOD1B data products
    product = {}
    product['atm'] = 'Atmospheric loading from {0}'.format(ATMOSPHERE)
    product['ocn'] = 'Oceanic loading from {0}'.format(OCEAN_MODEL)
    product['glo'] = 'Global atmospheric and oceanic loading'
    product['oba'] = 'Ocean bottom pressure from {0}'.format(OCEAN_MODEL)

    #-- GRACE AOD1B directory for data release
    aod1b_dir = os.path.join(base_dir,'AOD1B',DREL)
    #-- GRACE data directory for data release and processing center
    grace_dir = os.path.join(base_dir,PROC,DREL)
    #-- recursively create output directory if not currently existing
    if not os.access(os.path.join(grace_dir,DSET),os.F_OK):
        os.makedirs(os.path.join(grace_dir,DSET), MODE)
    #-- file formatting string if outputting to SHM format
    shm = '{0}-2_{1:4.0f}{2:03.0f}-{3:4.0f}{4:03.0f}_{5}_{6}_{7}_{8}00.gz'
    #-- center name if outputting to SHM format
    if (PROC == 'CSR'):
        CENTER = 'UTCSR'
    elif (PROC == 'GFZ'):
        CENTER = default_center
    elif (PROC == 'JPL'):
        CENTER = 'JPLEM'
    else:
        CENTER = default_center

    #-- read input DATE file from GSM data product
    grace_datefile = '{0}_{1}_DATES.txt'.format(PROC, DREL)
    date_input = np.loadtxt(os.path.join(grace_dir,'GSM',grace_datefile),
        skiprows=1)
    grace_month = date_input[:,1].astype(np.int64)
    start_yr = date_input[:,2]
    start_day = date_input[:,3].astype(np.int64)
    end_yr = date_input[:,4]
    end_day = date_input[:,5].astype(np.int64)
    #-- output date file reduced to months with complete AOD
    f_out = open(os.path.join(grace_dir,DSET,grace_datefile), 'w')
    #-- date file header information
    args = ('Mid-date','Month','Start_Day','End_Day','Total_Days')
    print('{0} {1:>10} {2:>11} {3:>10} {4:>13}'.format(*args),file=f_out)

    #-- for each GRACE/GRACE-FO month
    for t,gm in enumerate(grace_month):
        #-- check if GRACE/GRACE-FO month crosses years
        if (start_yr[t] != end_yr[t]):
            #-- check if start_yr is a Leap Year or Standard Year
            dpy = gravity_toolkit.time.calendar_days(start_yr[t]).sum()
            #-- list of Julian Days to read from both start and end year
            julian_days_to_read = []
            #-- add days to read from start and end years
            julian_days_to_read.extend([calc_julian_day(start_yr[t],D)
                for D in range(start_day[t],dpy+1)])
            julian_days_to_read.extend([calc_julian_day(end_yr[t],D)
                for D in range(1,end_day[t]+1)])
        else:
            #-- Julian Days to read going from start_day to end_day
            julian_days_to_read = [calc_julian_day(start_yr[t],D)
                for D in range(start_day[t],end_day[t]+1)]

        #-- output filename for GRACE/GRACE-FO month
        if (DATAFORM == 'SHM'):
            MISSION = 'GRAC' if (gm <= 186) else 'GRFO'
            FILE = shm.format(DSET.upper(),start_yr[t],start_day[t],
                end_yr[t],end_day[t],MISSION,CENTER,'BC01',DREL[2:])
        else:
            args = (PROC,DREL,DSET.upper(),LMAX,gm,suffix[DATAFORM])
            FILE = '{0}_{1}_{2}_CLM_L{3:d}_{4:03d}.{5}'.format(*args)

        #-- calendar dates to read
        JD = np.array(julian_days_to_read)
        Y,M,D,h,m,s = gravity_toolkit.time.convert_julian(JD,
            ASTYPE='i', FORMAT='tuple')
        #-- find unique year and month pairs to read
        rx1='|'.join(['{0:d}-{1:02d}'.format(*p) for p in set(zip(Y,M))])
        rx2='|'.join(['{0:0d}-{1:02d}-{2:02d}'.format(*p) for p in set(zip(Y,M,D))])
        #-- compile regular expressions operators for finding tar files
        tx = re.compile(r'AOD1B_({0})_\d+.(tar.gz|tgz)$'.format(rx1),re.VERBOSE)
        #-- finding all of the tar files in the AOD1b directory
        input_tar_files = [tf for tf in os.listdir(aod1b_dir) if tx.match(tf)]
        #-- compile regular expressions operators for file dates
        #-- will extract year and month and calendar day from the ascii file
        fx = re.compile(r'AOD1B_({0})_X_\d+.asc(.gz)?$'.format(rx2),re.VERBOSE)

        #-- check the last modified times of the tar file members
        input_mtime = np.zeros_like(julian_days_to_read,dtype=np.int64)
        input_file_check = np.zeros_like(julian_days_to_read,dtype=bool)
        c = 0
        #-- for each tar file
        for fi in sorted(input_tar_files):
            #-- open the AOD1B monthly tar file
            tar = tarfile.open(name=os.path.join(aod1b_dir,fi), mode='r:gz')
            #-- for each ascii file within the tar file that matches fx
            monthly_members = [m for m in tar.getmembers() if fx.match(m.name)]
            for member in monthly_members:
                #-- check last modification time of input tar file members
                input_mtime[c] = member.mtime
                input_file_check[c] = True
                c += 1

        #-- check if all files exist
        COMPLETE = input_file_check.all()
        #-- if output file exists: check if input tar file is newer
        TEST = False
        OVERWRITE = 'clobber'
        if os.access(os.path.join(grace_dir,DSET,FILE), os.F_OK):
            #-- check last modification time of input and output files
            output_mtime = os.stat(os.path.join(grace_dir,DSET,FILE)).st_mtime
            #-- if input tar file is newer: overwrite the output file
            if (input_mtime > output_mtime).any():
                TEST = True
                OVERWRITE = 'overwrite'
        else:
            TEST = True
            OVERWRITE = 'new'

        #-- print GRACE/GRACE-FO dates if there is a complete month of AOD
        if COMPLETE:
            #-- print GRACE/GRACE-FO dates to file
            print(('{0:13.8f} {1:03d} {2:8.0f} {3:03d} {4:8.0f} {5:03d} '
                '{6:8.0f}').format(date_input[t,0],gm,start_yr[t],start_day[t],
                end_yr[t],end_day[t],date_input[t,6]),file=f_out)

        #-- if there are new files, files to be rewritten or clobbered
        if COMPLETE and (TEST or CLOBBER):
            #-- if verbose: output information about the output file
            logging.info('{0} ({1})'.format(FILE,OVERWRITE))
            #-- allocate for the mean output harmonics
            Ylms = harmonics(lmax=LMAX, mmax=LMAX)
            nt = len(julian_days_to_read)*n_time
            Ylms.clm = np.zeros((LMAX+1,LMAX+1,nt))
            Ylms.slm = np.zeros((LMAX+1,LMAX+1,nt))
            Ylms.time = np.zeros((nt))
            count = 0
            #-- for each tar file
            for fi in sorted(input_tar_files):
                #-- open the AOD1B monthly tar file
                tar = tarfile.open(name=os.path.join(aod1b_dir,fi), mode='r:gz')
                #-- for each ascii file within the tar file that matches fx
                monthly_members=[m for m in tar.getmembers() if fx.match(m.name)]
                for member in monthly_members:
                    #-- extract member name
                    YMD,SFX = fx.findall(member.name).pop()
                    #-- open datafile for day
                    if (SFX == '.gz'):
                        fid = gzip.GzipFile(fileobj=tar.extractfile(member))
                    else:
                        fid = tar.extractfile(member)
                    #-- create counters for hour in dataset
                    hours = np.zeros((n_time))
                    c = 0
                    #-- while loop ends when dataset is read
                    while (c < n_time):
                        #-- read line
                        file_contents=fid.readline().decode('ISO-8859-1')
                        #-- find file header for data product
                        if bool(hx.search(file_contents)):
                            #-- extract hour from header and convert to float
                            HH, = re.findall(r'(\d+):\d+:\d+',file_contents)
                            hours[c] = np.int64(HH)
                            #-- read each line of spherical harmonics
                            for k in range(0,n_harm):
                                file_contents=fid.readline().decode('ISO-8859-1')
                                #-- find numerical instances in the data line
                                line_contents = rx.findall(file_contents)
                                #-- spherical harmonic degree and order
                                l1 = np.int64(line_contents[0])
                                m1 = np.int64(line_contents[1])
                                #-- spherical harmonic data saved to output Ylms
                                if (l1 <= LMAX) & (m1 <= LMAX):
                                    Ylms.clm[l1,m1,c]+=np.float64(line_contents[2])
                                    Ylms.slm[l1,m1,c]+=np.float64(line_contents[3])
                            #-- add 1 to hour counter
                            c += 1
                    #-- close the input file for day
                    fid.close()
                    #-- year fraction of the particular date and times
                    YEAR = np.repeat(Y[count//n_time], n_time).astype('f')
                    MONTH = np.repeat(M[count//n_time], n_time).astype('f')
                    DAY = np.repeat(D[count//n_time], n_time).astype('f')
                    Ylms.time[count:count+n_time] = \
                        gravity_toolkit.time.convert_calendar_decimal(YEAR,
                        MONTH, day=DAY, hour=hours)
                    #-- add to day counter
                    count += n_time

            #-- calculate mean harmonics for GRACE/GRACE-FO month
            #-- convert from harmonics object to dealiasing object
            mean_Ylms = dealiasing().from_harmonics(Ylms.mean())
            mean_Ylms.time = np.mean(Ylms.time)
            mean_Ylms.month = np.int64(gm)
            #-- product information
            mean_Ylms.center = PROC
            mean_Ylms.release = DREL
            mean_Ylms.product = DSET
            #-- start and end time for month
            start_time = gravity_toolkit.time.convert_julian(np.min(JD))
            mean_Ylms.start_time = ['{0:4.0f}'.format(start_time['year']),
                '{0:02.0f}'.format(start_time['month']),
                '{0:02.0f}'.format(start_time['day'])]
            end_time = gravity_toolkit.time.convert_julian(np.max(JD))
            mean_Ylms.end_time = ['{0:4.0f}'.format(end_time['year']),
                '{0:02.0f}'.format(end_time['month']),
                '{0:02.0f}'.format(end_time['day'])]
            #-- output mean Ylms to file
            if (DATAFORM == 'ascii'):
                #-- ascii (.txt)
                mean_Ylms.to_ascii(os.path.join(grace_dir,DSET,FILE))
            elif (DATAFORM == 'netCDF4'):
                #-- netcdf (.nc)
                mean_Ylms.to_netCDF4(os.path.join(grace_dir,DSET,FILE))
            elif (DATAFORM == 'HDF5'):
                #-- HDF5 (.H5)
                mean_Ylms.to_HDF5(os.path.join(grace_dir,DSET,FILE))
            elif (DATAFORM == 'SHM'):
                mean_Ylms.to_SHM(os.path.join(grace_dir,DSET,FILE))
            #-- set the permissions mode of the output file
            os.chmod(os.path.join(grace_dir,DSET,FILE), MODE)

        elif not COMPLETE:
            logging.info('File {0} not output (incomplete)'.format(FILE))

    #-- if outputting as spherical harmonic model files
    if (DATAFORM == 'SHM'):
        #-- Create an index file for each GRACE product
        grace_files = [fi for fi in os.listdir(os.path.join(grace_dir,DSET)) if
            re.match(r'{0}-2(.*?)\.gz'.format(DSET),fi)]
        #-- outputting GRACE filenames to index
        with open(os.path.join(grace_dir,DSET,'index.txt'),'w') as fid:
            for fi in sorted(grace_files):
                print(fi, file=fid)
        #-- change permissions of index file
        os.chmod(os.path.join(grace_dir,DSET,'index.txt'), MODE)

    #-- print completion flag
    logging.info('Complete: {0}/{1}/{2}'.format(PROC,DREL,DSET))
    #-- close the output date file
    f_out.close()

#-- PURPOSE: additional routines for the harmonics module
class dealiasing(harmonics):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.center=None
        self.release='RLxx'
        self.product=None
        self.start_time=[None]*3
        self.end_time=[None]*3

    def from_harmonics(self, temp):
        """
        Convert a harmonics object to a new dealiasing object
        """
        self = dealiasing(lmax=temp.lmax, mmax=temp.mmax)
        #-- try to assign variables to self
        for key in ['clm','slm','time','month','shape','ndim','filename',
            'center','release','product','start_time','end_time']:
            try:
                val = getattr(temp, key)
                setattr(self, key, np.copy(val))
            except AttributeError:
                pass
        #-- assign ndim and shape attributes
        self.update_dimensions()
        return self

    def to_SHM(self, filename, **kwargs):
        """
        Write a harmonics object to SHM file
        Inputs: full path of output SHM file
        Options:
            harmonics objects contain date information
            keyword arguments for SHM output
        """
        self.filename = os.path.expanduser(filename)
        #-- set default verbosity
        kwargs.setdefault('verbose',False)
        logging.info(self.filename)
        #-- open the output file
        fid = gzip.open(self.filename, 'wt')
        #-- print the header informat
        self.print_header(fid)
        self.print_harmonic(fid)
        self.print_global(fid)
        self.print_variables(fid,'double precision')
        #-- output file format
        file_format = ('{0:6} {1:4d} {2:4d} {3:+18.12E} {4:+18.12E} '
            '{5:10.4E} {6:10.4E} {7} {8} {9}')
        #-- start and end time in line format
        start_date = '{0}{1}{2}.0000'.format(*self.start_time)
        end_date = '{0}{1}{2}.0000'.format(*self.end_time)
        #-- write to file for each spherical harmonic degree and order
        for m in range(0, self.mmax+1):
            for l in range(m, self.lmax+1):
                args = ('GRCOF2', l, m, self.clm[l,m], self.slm[l,m],
                    0, 0, start_date, end_date, 'nnnn')
                print(file_format.format(*args), file=fid)
        #-- close the output file
        fid.close()

    #-- PURPOSE: print YAML header to top of file
    def print_header(self, fid):
        #-- print header
        fid.write('{0}:\n'.format('header'))
        #-- data dimensions
        fid.write('  {0}:\n'.format('dimensions'))
        fid.write('    {0:22}: {1:d}\n'.format('degree',self.lmax))
        fid.write('    {0:22}: {1:d}\n'.format('order',self.lmax))
        fid.write('\n')

    #-- PURPOSE: print spherical harmonic attributes to YAML header
    def print_harmonic(self, fid):
        #-- non-standard attributes
        fid.write('  {0}:\n'.format('non-standard_attributes'))
        #-- product id
        product_id = '{0}-2'.format(self.product)
        fid.write('    {0:22}: {1}\n'.format('product_id',product_id))
        #-- format id
        fid.write('    {0:22}:\n'.format('format_id'))
        short_name = 'SHM'
        fid.write('      {0:20}: {1}\n'.format('short_name',short_name))
        long_name = 'Earth Gravity Spherical Harmonic Model Format'
        fid.write('      {0:20}: {1}\n'.format('long_name',long_name))
        #-- harmonic normalization
        normalization = 'fully normalized'
        fid.write('    {0:22}: {1}\n'.format('normalization',
            normalization))
        #-- earth parameters
        #-- gravitational constant
        fid.write('    {0:22}:\n'.format('earth_gravity_param'))
        long_name = 'gravitational constant times mass of Earth'
        fid.write('      {0:20}: {1}\n'.format('long_name',long_name))
        units = 'm3/s2'
        fid.write('      {0:20}: {1}\n'.format('units',units))
        value = '3.9860044180E+14'
        fid.write('      {0:20}: {1}\n'.format('value',value))
        #-- equatorial radius
        fid.write('    {0:22}:\n'.format('mean_equator_radius'))
        long_name = 'mean equator radius'
        fid.write('      {0:20}: {1}\n'.format('long_name',long_name))
        units = 'meters'
        fid.write('      {0:20}: {1}\n'.format('units',units))
        value = '6.3781366000E+06'
        fid.write('      {0:20}: {1}\n'.format('value',value))
        fid.write('\n')

    #-- PURPOSE: print global attributes to YAML header
    def print_global(self,fid):
        fid.write('  {0}:\n'.format('global_attributes'))
        #-- product title
        if (self.month <= 186):
            MISSION = 'GRACE'
            PROJECT = 'NASA Gravity Recovery And Climate Experiment (GRACE)'
            ACKNOWLEDGEMENT = ('GRACE is a joint mission of NASA (USA) and '
                'DLR (Germany).')
        else:
            MISSION = 'GRACE-FO'
            PROJECT = ('NASA Gravity Recovery And Climate Experiment '
                'Follow-On (GRACE-FO)')
            ACKNOWLEDGEMENT = ('GRACE-FO is a joint mission of the US National '
                'Aeronautics and Space Administration and the German Research '
                'Center for Geosciences.')
        args = (MISSION,self.product,self.center,self.release)
        title = '{0} Geopotential {1} Coefficients {2} {3}'.format(*args)
        fid.write('    {0:22}: {1}\n'.format('title',title))
        #-- product summaries
        summaries = {}
        summaries['GAA'] = ("Spherical harmonic coefficients that represent "
            "anomalous contributions of the non-tidal atmosphere to the Earth's "
            "mean gravity field during the specified timespan. This includes the "
            "contribution of atmospheric surface pressure over the continents, "
            "the static contribution of atmospheric pressure to ocean bottom "
            "pressure elsewhere, and the contribution of upper-air density "
            "anomalies above both the continents and the oceans.")
        summaries['GAB'] = ("Spherical harmonic coefficients that represent "
            "anomalous contributions of the non-tidal dynamic ocean to ocean "
            "bottom pressure during the specified timespan.")
        summaries['GAC'] = ("Spherical harmonic coefficients that represent "
            "the sum of the ATM (or GAA) and OCN (or GAB) coefficients during "
            "the specified timespan. These coefficients represent anomalous "
            "contributions of the non-tidal dynamic ocean to ocean bottom "
            "pressure, the non-tidal atmospheric surface pressure over the "
            "continents, the static contribution of atmospheric pressure to "
            "ocean bottom pressure, and the upper-air density anomalies above "
            "both the continents and the oceans.")
        summaries['GAD'] = ("Spherical harmonic coefficients that are zero "
            "over the continents, and provide the anomalous simulated ocean "
            "bottom pressure that includes non-tidal air and water "
            "contributions elsewhere during the specified timespan. These "
            "coefficients differ from GLO (or GAC) coefficients over the "
            "ocean domain by disregarding upper air density anomalies.")
        summary = summaries[self.product]
        fid.write('    {0:22}: {1}\n'.format('summary',''.join(summary)))
        fid.write('    {0:22}: {1}\n'.format('project',PROJECT))
        keywords = []
        keywords.append('GRACE')
        keywords.append('GRACE-FO') if (self.month > 186) else None
        keywords.append('Level-2')
        keywords.append('SHM')
        keywords.append('Spherical Harmonic Model')
        keywords.append('Gravitational Field')
        keywords.append('GSM')
        keywords.append('Geopotential')
        keywords.append('Gravity Field')
        keywords.append('Mass')
        keywords.append('Mass Transport')
        keywords.append('Total Water Storage')
        keywords.append('Time Variable Gravity')
        keywords.append('Mass Balance')
        keywords.append('Gravity Anomaly')
        keywords.append('Satellite Geodesy')
        keywords.append('Ocean')
        keywords.append('Ocean Bottom Pressure')
        keywords.append('AOD')
        keywords.append('Atmosphere')
        keywords.append('Non-tidal Atmosphere')
        keywords.append('Dealiasing Product')
        fid.write('    {0:22}: {1}\n'.format('keywords',', '.join(keywords)))
        vocabulary = 'NASA Global Change Master Directory (GCMD) Science Keywords'
        fid.write('    {0:22}: {1}\n'.format('keywords_vocabulary',vocabulary))
        if (self.center == 'CSR'):
            institution = 'UT-AUSTIN/CSR'
        elif (self.center == 'GFZ'):
            institution = 'GFZ German Research Centre for Geosciences'
        elif (self.center == 'JPL'):
            institution = 'NASA/JPL'
        else:
            #-- default to GFZ
            institution = 'GFZ German Research Centre for Geosciences'
        fid.write('    {0:22}: {1}\n'.format('institution',institution))
        src = 'All data from AOD1B {0}'.format(self.release)
        fid.write('    {0:22}: {1}\n'.format('source',src))
        fid.write('    {0:22}: {1:d}\n'.format('processing_level',2))
        fid.write('    {0:22}: {1}\n'.format('acknowledgement',ACKNOWLEDGEMENT))
        PRODUCT_VERSION = 'Release-{0}'.format(self.release[2:])
        fid.write('    {0:22}: {1}\n'.format('product_version',PRODUCT_VERSION))
        fid.write('    {0:22}:\n'.format('references'))
        #-- date range and date created
        start_date = '{0}-{1}-{2}'.format(*self.start_time)
        fid.write('    {0:22}: {1}\n'.format('time_coverage_start',start_date))
        end_date = '{0}-{1}-{2}'.format(*self.end_time)
        fid.write('    {0:22}: {1}\n'.format('time_coverage_end',end_date))
        today = time.strftime('%Y-%m-%d',time.localtime())
        fid.write('    {0:22}: {1}\n'.format('date_created', today))
        fid.write('\n')

    #-- PURPOSE: print variable descriptions to YAML header
    def print_variables(self,fid,data_precision):
        #-- variables
        fid.write('  {0}:\n'.format('variables'))
        #-- record_key
        fid.write('    {0:22}:\n'.format('record_key'))
        long_name = 'Earth Gravity Spherical Harmonic Model Format Type 2'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'string'))
        fid.write('      {0:20}: {1}\n'.format('comment', '1st column'))
        #-- degree_index
        fid.write('    {0:22}:\n'.format('degree_index'))
        long_name = 'spherical harmonic degree l'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
        fid.write('      {0:20}: {1}\n'.format('comment', '2nd column'))
        #-- order_index
        fid.write('    {0:22}:\n'.format('order_index'))
        long_name = 'spherical harmonic order m'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'integer'))
        fid.write('      {0:20}: {1}\n'.format('comment', '3rd column'))
        #-- clm
        fid.write('    {0:22}:\n'.format('clm'))
        long_name = 'Clm coefficient; cosine coefficient for degree l and order m'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '4th column'))
        #-- slm
        fid.write('    {0:22}:\n'.format('slm'))
        long_name = 'Slm coefficient; sine coefficient for degree l and order m'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '5th column'))
        #-- clm_std_dev
        fid.write('    {0:22}:\n'.format('clm_std_dev'))
        long_name = 'standard deviation of Clm'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '6th column'))
        #-- slm_std_dev
        fid.write('    {0:22}:\n'.format('slm_std_dev'))
        long_name = 'standard deviation of Slm'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '7th column'))
        #-- epoch_begin_time
        fid.write('    {0:22}:\n'.format('epoch_begin_time'))
        long_name = 'epoch begin of Clm, Slm coefficients'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('time_format', 'yyyymmdd.hhmm'))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '8th column'))
        #-- epoch_stop_time
        fid.write('    {0:22}:\n'.format('epoch_stop_time'))
        long_name = 'epoch stop of Clm, Slm coefficients'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('time_format', 'yyyymmdd.hhmm'))
        fid.write('      {0:20}: {1}\n'.format('data_type', data_precision))
        fid.write('      {0:20}: {1}\n'.format('comment', '9th column'))
        #-- solution_flags
        fid.write('    {0:22}:\n'.format('solution_flags'))
        long_name = 'Coefficient adjustment and a priori flags'
        fid.write('      {0:20}: {1}\n'.format('long_name', long_name))
        fid.write('      {0:20}: {1}\n'.format('coverage_content_type',
            'auxiliaryInformation'))
        fid.write('      {0:20}: {1}\n'.format('data_type', 'byte'))
        fid.write('      {0:20}:\n'.format('flag_meanings'))
        #-- solution flag meanings
        m = []
        m.append('Clm adjusted, y for yes and n for no')
        m.append('Slm adjusted, y for yes and n for no')
        m.append('stochastic a priori info for Clm, y for yes and n for no')
        m.append('stochastic a priori info for Slm, y for yes and n for no')
        for i,meaning in enumerate(m):
            fid.write('        - char {0:d} = {1}\n'.format(i, meaning))
        fid.write('      {0:20}: {1}\n'.format('comment', '10th column'))
        #-- end of header
        fid.write('\n\n# End of YAML header\n')

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads GRACE/GRACE-FO AOD1B datafiles for a
            specific product and outputs monthly mean for a specific
            GRACE/GRACE-FO processing center and data release
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Data processing center or satellite mission
    parser.add_argument('--center','-c',
        metavar='PROC', type=str, required=True,
        help='GRACE/GRACE-FO Processing Center')
    #-- GRACE/GRACE-FO data release
    parser.add_argument('--release','-r',
        metavar='DREL', type=str, default='RL06',
        help='GRACE/GRACE-FO Data Release')
    #-- GRACE/GRACE-FO dealiasing product
    parser.add_argument('--product','-p',
        metavar='DSET', type=str.upper, nargs='+',
        choices=['GAA','GAB','GAC','GAD'],
        help='GRACE/GRACE-FO dealiasing product')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=180,
        help='Maximum spherical harmonic degree')
    #-- input and output data format (ascii, netCDF4, HDF5, SHM)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4',
        choices=['ascii','netCDF4','HDF5','SHM'],
        help='Output data format')
    #-- clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    for DSET in args.product:
        #-- run monthly mean AOD1b program with parameters
        dealiasing_monthly_mean(args.directory,
            PROC=args.center,
            DREL=args.release,
            DSET=DSET,
            LMAX=args.lmax,
            DATAFORM=args.format,
            CLOBBER=args.clobber,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
