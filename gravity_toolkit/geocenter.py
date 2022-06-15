#!/usr/bin/env python
u"""
geocenter.py
Written by Tyler Sutterley (06/2022)
Data class for reading and processing geocenter data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    PyYAML: YAML parser and emitter for Python
        https://github.com/yaml/pyyaml

UPDATE HISTORY:
    Updated 06/2022: drop external reader dependency for UCI format
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 03/2022: add try/except for read_GRACE_geocenter
    Updated 12/2021: added netCDF4 reader for UCI iteration files
        add cartesian and surface mass density conversions for errors
        logging case_insensitive_filename output for debugging
    Updated 11/2021: converted to class with all data readers and converters
    Updated 07/2020: added function docstrings
    Updated 06/2019: added option RADIUS to manually set the Earth's radius
    Updated 04/2017: changed option from INV to INVERSE and made True/False
    Updated 04/2015: calculate radius of the Earth directly in program
    Updated 02/2014: minor update to if statement
    Updated 03/2013: converted to python
"""
import os
import re
import io
import copy
import gzip
import time
import uuid
import yaml
import logging
import netCDF4
import numpy as np
import gravity_toolkit.time

class geocenter(object):
    """
    Data class for reading and processing geocenter data

    Attributes
    ----------
    C10: float
        cosine spherical harmonics of degree 1 and order 0
    C11: float
        cosine spherical harmonics of degree 1 and order 1
    S11: float
        sine spherical harmonics of degree 1 and order 1
    X: float
        X-component of Cartesian geocenter coordinates
    Y: float
        Y-component of Cartesian geocenter coordinates
    Z: float
        Z-component of Cartesian geocenter coordinates
    time: float
        time variable of the spherical harmonics
    month: int
        GRACE/GRACE-FO months variable of the spherical harmonics
    radius: float, default 6371000.790009159
        Average Radius of the Earth [mm]
    """
    np.seterr(invalid='ignore')
    def __init__(self, **kwargs):
        #-- WGS84 ellipsoid parameters
        a_axis = 6378137.0#-- [m] semimajor axis of the ellipsoid
        flat = 1.0/298.257223563#-- flattening of the ellipsoid
        #-- Mean Earth's Radius in mm having the same volume as WGS84 ellipsoid
        kwargs.setdefault('radius', 1000.0*a_axis*(1.0 - flat)**(1.0/3.0))
        #-- cartesian coordinates
        kwargs.setdefault('X',None)
        kwargs.setdefault('Y',None)
        kwargs.setdefault('Z',None)
        #-- set default class attributes
        self.C10=None
        self.C11=None
        self.S11=None
        self.X=copy.copy(kwargs['X'])
        self.Y=copy.copy(kwargs['Y'])
        self.Z=copy.copy(kwargs['Z'])
        self.time=None
        self.month=None
        self.filename=None
        #-- Average Radius of the Earth [mm]
        self.radius=copy.copy(kwargs['radius'])

    def case_insensitive_filename(self,filename):
        """
        Searches a directory for a filename without case dependence

        Parameters
        ----------
        filename: str
            input filename
        """
        #-- check if filename is open file object
        if isinstance(filename, io.IOBase):
            self.filename = copy.copy(filename)
        else:
            #-- tilde-expand input filename
            self.filename = os.path.expanduser(filename)
            #-- check if file presently exists with input case
            if not os.access(self.filename,os.F_OK):
                #-- search for filename without case dependence
                basename = os.path.basename(filename)
                directory = os.path.dirname(os.path.expanduser(filename))
                f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
                #-- check that geocenter file exists
                if not f:
                    errmsg = '{0} not found in file system'.format(filename)
                    raise FileNotFoundError(errmsg)
                self.filename = os.path.join(directory,f.pop())
        #-- print filename
        logging.debug(self.filename)
        return self

    #-- PURPOSE: read AOD1b geocenter for month and calculate the mean harmonics
    #-- need to run aod1b_geocenter.py to write these monthly geocenter files
    def from_AOD1B(self, release, calendar_year, calendar_month):
        """
        Reads monthly non-tidal ocean and atmospheric variation geocenter files

        Parameters
        ----------
        release: str
            GRACE/GRACE-FO/Swarm data release for dealiasing product
        calendar_year: int
            calendar year of data
        calendar_month: int
            calendar month of data
        """

        #-- full path to AOD geocenter for month (using glo coefficients)
        args = (release,'glo',calendar_year,calendar_month)
        AOD1B_file = 'AOD1B_{0}_{1}_{2:4.0f}_{3:02.0f}.txt'.format(*args)
        #-- check that file exists
        if not os.access(os.path.join(self.directory,AOD1B_file), os.F_OK):
            errmsg = 'AOD1B File {0} not in File System'.format(AOD1B_file)
            raise FileNotFoundError(errmsg)
        #-- read AOD1b geocenter skipping over commented header text
        with open(os.path.join(self.directory,AOD1B_file), mode='r', encoding='utf8') as f:
            file_contents=[i for i in f.read().splitlines() if not re.match(r'#',i)]
        #-- extract X,Y,Z from each line in the file
        n_lines = len(file_contents)
        temp = geocenter()
        temp.X = np.zeros((n_lines))
        temp.Y = np.zeros((n_lines))
        temp.Z = np.zeros((n_lines))
        for i,line in enumerate(file_contents):
            line_contents = line.split()
            #-- first column: ISO-formatted date and time
            cal_date = time.strptime(line_contents[0],r'%Y-%m-%dT%H:%M:%S')
            #-- verify that dates are within year and month
            assert (cal_date.tm_year == calendar_year)
            assert (cal_date.tm_mon == calendar_month)
            #-- second-fourth columns: X, Y and Z geocenter variations
            temp.X[i],temp.Y[i],temp.Z[i] = np.array(line_contents[1:],dtype='f')
        #-- convert X,Y,Z into spherical harmonics
        temp.from_cartesian()
        #-- return the spherical harmonic coefficients
        return temp

    def from_gravis(self, geocenter_file, **kwargs):
        """
        Reads monthly geocenter spherical harmonic data files from
        `GFZ GravIS calculated using GRACE/GRACE-FO measurements
        and Ocean Models of degree 1
        <ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/GRAVIS-2B_GFZOP_GEOCENTER_0002.dat>`_


        Parameters
        ----------
        geocenter_file: str
            degree 1 file
        header: bool, default True
            file contains header text to be skipped

        References
        ----------
        .. [Dahle2019] Dahle and Murboeck, "Post-processed GRACE/GRACE-FO
            Geopotential GSM Coefficients GFZ RL06 (Level-2B Product)."
            V. 0002. *GFZ Data Services*, (2019).
            `doi: 10.5880/GFZ.GRAVIS_06_L2B <https://doi.org/10.5880/GFZ.GRAVIS_06_L2B>`_
        """

        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- set default keyword arguments
        kwargs.setdefault('header',True)

        #-- Combined GRACE/SLR geocenter solution file produced by GFZ GravIS
        #-- Column  1: MJD of BEGINNING of solution data span
        #-- Column  2: Year and fraction of year of BEGINNING of solution data span
        #-- Column  3: Coefficient C(1,0)
        #-- Column  4: Coefficient C(1,0) - mean C(1,0) (1.0E-10)
        #-- Column  5: C(1,0) uncertainty (1.0E-10)
        #-- Column  6: Coefficient C(1,1)
        #-- Column  7: Coefficient C(1,1) - mean C(1,1) (1.0E-10)
        #-- Column  8: C(1,1) uncertainty (1.0E-10)
        #-- Column  9: Coefficient S(1,1)
        #-- Column 10: Coefficient S(1,1) - mean S(1,1) (1.0E-10)
        #-- Column 11: S(1,1) uncertainty (1.0E-10)

        with open(self.filename, mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        #-- number of lines contained in the file
        file_lines = len(file_contents)

        #-- counts the number of lines in the header
        count = 0
        #-- Reading over header text
        while kwargs['header']:
            #-- file line at count
            line = file_contents[count]
            #-- find PRODUCT: within line to set HEADER flag to False when found
            kwargs['header'] = not bool(re.match(r'PRODUCT:+',line))
            #-- add 1 to counter
            count += 1

        #-- output dictionary with spherical harmonic solutions
        dinput = {}
        #-- number of months within the file
        n_mon = file_lines - count
        #-- date and GRACE/GRACE-FO month
        dinput['time'] = np.zeros((n_mon))
        dinput['month'] = np.zeros((n_mon),dtype=int)
        #-- monthly spherical harmonic replacement solutions
        dinput['C10'] = np.zeros((n_mon))
        dinput['C11'] = np.zeros((n_mon))
        dinput['S11'] = np.zeros((n_mon))
        #-- monthly spherical harmonic formal standard deviations
        dinput['eC10'] = np.zeros((n_mon))
        dinput['eC11'] = np.zeros((n_mon))
        dinput['eS11'] = np.zeros((n_mon))
        #-- time count
        t = 0
        #-- for every other line:
        for line in file_contents[count:]:
            #-- find numerical instances in line including exponents,
            #-- decimal points and negatives
            line_contents = re.findall(r'[-+]?\d*\.\d*(?:[eE][-+]?\d+)?',line)
            count = len(line_contents)
            #-- check for empty lines
            if (count > 0):
                #-- reading decimal year for start of span
                dinput['time'][t] = np.float64(line_contents[1])
                #-- Spherical Harmonic data for line
                dinput['C10'][t] = np.float64(line_contents[2])
                dinput['C11'][t] = np.float64(line_contents[5])
                dinput['S11'][t] = np.float64(line_contents[8])
                #-- monthly spherical harmonic formal standard deviations
                dinput['eC10'][t] = np.float64(line_contents[4])*1e-10
                dinput['eC11'][t] = np.float64(line_contents[7])*1e-10
                dinput['eS11'][t] = np.float64(line_contents[10])*1e-10
                #-- GRACE/GRACE-FO month of geocenter solutions
                dinput['month'][t] = gravity_toolkit.time.calendar_to_grace(
                    dinput['time'][t], around=np.round)
                #-- add to t count
                t += 1
        #-- truncate variables if necessary
        for key,val in dinput.items():
            dinput[key] = val[:t]

        #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
        #-- Accelerometer shutoffs make the relation between month number
        #-- and date more complicated as days from other months are used
        #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
        #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
        #-- For all: May 2015 (161) is centered in Apr 2015 (160)
        #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
        dinput['month'] = gravity_toolkit.time.adjust_months(dinput['month'])

        #-- return the GFZ GravIS geocenter solutions
        return self.from_dict(dinput)


    def from_SLR(self, geocenter_file, **kwargs):
        """
        Reads monthly geocenter files from satellite laser ranging corrected
        for non-tidal ocean and atmospheric variation

        Reads monthly geocenter files from `satellite laser ranging
        provided by CSR <http://download.csr.utexas.edu/pub/slr/geocenter/>`_

            - `RL04`: GCN_RL04.txt
            - `RL05`: GCN_RL05.txt

        `New CF-CM geocenter dataset
        <http://download.csr.utexas.edu/pub/slr/geocenter/GCN_L1_L2_30d_CF-CM.txt>`_
        to reflect the `true degree-1 mass variations
        <http://download.csr.utexas.edu/pub/slr/geocenter/geocenter/README_L1_L2>`_

        `New geocenter solutions from Minkang Cheng
        <http://download.csr.utexas.edu/outgoing/cheng/gct2est.220_5s>`_

        Parameters
        ----------
        geocenter_file: str
            Satellite Laser Ranging file
        AOD: bool, default False
            remove Atmospheric and Oceanic Dealiasing products
        release: str or NoneType, default None
            GRACE/GRACE-FO/Swarm data release for AOD
        header: int, default 0
            rows of data to skip when importing data
        columns: list, default []
            column names of ascii file

                - ``'time'``: date in decimal-years
                - ``'X'``: X-component of geocenter variation
                - ``'Y'``: Y-component of geocenter variation
                - ``'Z'``: Z-component of geocenter variation
                - ``'X_sigma'``: X-component uncertainty
                - ``'Y_sigma'``: Y-component uncertainty
                - ``'Z_sigma'``: Z-component uncertainty
        """

        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- set default keyword arguments
        kwargs.setdefault('AOD',False)
        kwargs.setdefault('columns',[])
        kwargs.setdefault('header',0)
        kwargs.setdefault('release',None)
        #-- copy keyword arguments to variables
        COLUMNS = copy.copy(kwargs['columns'])
        HEADER = copy.copy(kwargs['header'])

        #-- directory setup for AOD1b data starting with input degree 1 file
        #-- this will verify that the input paths work
        base_dir = os.path.join(os.path.dirname(self.filename),os.path.pardir)
        self.directory = os.path.abspath(os.path.join(base_dir,'AOD1B',
            kwargs['release'],'geocenter'))
        #-- check that AOD1B directory exists
        if not os.access(self.directory, os.F_OK):
            errmsg = '{0} not found in file system'.format(self.directory)
            raise FileNotFoundError(errmsg)

        #-- Input geocenter file and split lines
        with open(os.path.expanduser(geocenter_file), mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        ndate = len(file_contents) - HEADER

        #-- compile regular expression operator to find numerical instances
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)

        #-- initializing output data
        #-- Degree 1 Stokes Coefficients
        self.C10 = np.zeros((ndate))
        self.C11 = np.zeros((ndate))
        self.S11 = np.zeros((ndate))
        #-- Degree 1 Stokes Coefficient Errors
        self.eC10 = np.zeros((ndate))
        self.eC11 = np.zeros((ndate))
        self.eS11 = np.zeros((ndate))
        #-- Date information
        self.time = np.zeros((ndate))
        self.month = np.zeros((ndate), dtype=np.int32)
        JD = np.zeros((ndate))

        #-- for each date
        for t,file_line in enumerate(file_contents[HEADER:]):
            #-- find numerical instances in line
            #-- replacing fortran double precision exponential
            line_contents = rx.findall(file_line.replace('D','E'))
            #-- extract date
            self.time[t] = np.float64(line_contents[COLUMNS.index('time')])

            #-- extract geocenter variations
            temp = geocenter(radius=self.radius)
            temp.X = np.float64(line_contents[COLUMNS.index('X')])
            temp.Y = np.float64(line_contents[COLUMNS.index('Y')])
            temp.Z = np.float64(line_contents[COLUMNS.index('Z')])
            temp.from_cartesian()
            #-- copy spherical harmonics to output
            self.C10[t] = np.copy(temp.C10)
            self.C11[t] = np.copy(temp.C11)
            self.S11[t] = np.copy(temp.S11)

            #-- extract geocenter uncertainties
            temp = geocenter(radius=self.radius)
            temp.X = np.float64(line_contents[COLUMNS.index('X_sigma')])
            temp.Y = np.float64(line_contents[COLUMNS.index('Y_sigma')])
            temp.Z = np.float64(line_contents[COLUMNS.index('Z_sigma')])
            temp.from_cartesian()
            #-- copy spherical harmonic uncertainties to output
            self.eC10[t] = np.copy(temp.C10)
            self.eC11[t] = np.copy(temp.C11)
            self.eS11[t] = np.copy(temp.S11)

            #-- Calculation of the Julian date from calendar date
            JD[t] = gravity_toolkit.time.calendar_to_julian(self.time[t])
            #-- convert the julian date into calendar dates
            YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(JD[t],
                FORMAT='tuple')
            #-- calculate the GRACE/GRACE-FO month (Apr02 == 004)
            #-- https://grace.jpl.nasa.gov/data/grace-months/
            self.month[t] = gravity_toolkit.time.calendar_to_grace(YY,month=MM)

            #-- if removing the Atmospheric and Oceanic dealiasing
            if kwargs['AOD']:
                #-- read the AOD1B file for the month and year
                temp = self.from_AOD1B(kwargs['release'], YY, MM)
                #-- remove the monthly mean AOD
                self.C10[t] -= np.mean(temp.C10)
                self.C11[t] -= np.mean(temp.C11)
                self.S11[t] -= np.mean(temp.S11)

        #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
        #-- Accelerometer shutoffs make the relation between month number
        #-- and date more complicated as days from other months are used
        #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
        #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
        #-- For all: May 2015 (161) is centered in Apr 2015 (160)
        #-- For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
        self.month = gravity_toolkit.time.adjust_months(self.month)
        #-- return the geocenter harmonics
        return self

    def from_UCI(self, geocenter_file, **kwargs):
        """
        Reads monthly geocenter files computed using GRACE/GRACE-FO
        measurements and ocean models [Swenson2008]_ [Sutterley2019]_

        Parameters
        ----------
        geocenter_file: str
            input datafile with geocenter coefficients

        References
        ----------
        .. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr,
            "Estimating geocenter variations from a combination
            of GRACE and ocean model output", *Journal of Geophysical
            Research*, 113(B08410), (2008).
            `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
        .. [Sutterley2019] T. C. Sutterley, and I. Velicogna, "Improved
            estimates of geocenter variability from time-variable gravity
            and ocean model outputs", *Remote Sensing*, 11(18), 2108, (2019).
            `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_
        """
        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- read geocenter file and get contents
        with open(os.path.expanduser(geocenter_file), mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        #-- number of lines contained in the file
        file_lines = len(file_contents)

        #-- counts the number of lines in the header
        HEADER = False
        count = 0
        #-- Reading over header text
        while (HEADER is False) and (count < file_lines):
            #-- file line at count
            line = file_contents[count]
            #--if End of YAML Header is found: set HEADER flag
            HEADER = bool(re.search("\# End of YAML header",line))
            #-- add 1 to counter
            count += 1

        #-- verify HEADER flag was set
        if not HEADER:
            raise IOError('Data not found in file:\n\t{0}'.format(geocenter_file))

        #-- number of months within the file
        n_mon = np.int64(file_lines - count)
        #-- output time variables
        DEG1 = {}
        DEG1['time'] = np.zeros((n_mon))
        DEG1['JD'] = np.zeros((n_mon))
        DEG1['month'] = np.zeros((n_mon), dtype=np.int64)
        #-- parse the YAML header (specifying yaml loader)
        DEG1.update(yaml.load('\n'.join(file_contents[:count]),
            Loader=yaml.BaseLoader))

        #-- compile numerical expression operator
        regex_pattern = '[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)

        #-- get names and columns of input variables
        variables = copy.copy(DEG1['header']['variables'])
        variables.pop('mid-epoch_time')
        variables.pop('month')
        columns = {}
        #-- for each output data variable
        for key in variables:
            DEG1[key] = np.zeros((n_mon))
            comment_text, = rx.findall(variables[key]['comment'])
            columns[key] = int(comment_text) - 1

        #-- for every other line:
        for t, line in enumerate(file_contents[count:]):
            #-- find numerical instances in line including integers, exponents,
            #-- decimal points and negatives
            line_contents = rx.findall(line)
            #-- extacting mid-date time and GRACE/GRACE-FO "month"
            DEG1['time'][t] = np.float64(line_contents[0])
            DEG1['month'][t] = np.int64(line_contents[-1])
            #-- calculate mid-date as Julian dates
            #-- calendar year of date
            year = np.floor(DEG1['time'][t])
            #-- check if year is a leap year
            days_per_year = np.sum(gravity_toolkit.time.calendar_days(year))
            #-- calculation of day of the year
            day_of_the_year = days_per_year*(DEG1['time'][t] % 1)
            #-- calculate Julian day
            DEG1['JD'][t] = np.float64(367.0*year - np.floor(7.0*(year)/4.0) -
                np.floor(3.0*(np.floor((year - 8.0/7.0)/100.0) + 1.0)/4.0) +
                np.floor(275.0/9.0) + day_of_the_year + 1721028.5)
            #-- extract fully-normalized degree one spherical harmonics
            for key,val in columns.items():
                DEG1[key][t] = np.float64(line_contents[val])

        #-- return the geocenter harmonics
        return self.from_dict(DEG1)

    def from_swenson(self, geocenter_file, **kwargs):
        """
        Reads `monthly geocenter coefficients
        <https://github.com/swensosc/GRACE_Tiles/blob/master/ancillary_data/gad_gsm.rl05.txt>`_
        computed by Sean Swenson using GRACE/GRACE-FO measurements
        and Ocean Models of degree 1

        Parameters
        ----------
        geocenter_file: str
            degree 1 file
        header: bool, default True
            file contains header text to be skipped

        References
        ----------
        .. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr,
            "Estimating geocenter variations from a combination
            of GRACE and ocean model output", *Journal of Geophysical
            Research*, 113(B08410), (2008).
            `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
        """
        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- set default keyword arguments
        kwargs.setdefault('header',True)

        #-- read degree 1 file and get contents
        with open(self.filename, mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        #-- number of lines contained in the file
        file_lines = len(file_contents)

        #-- counts the number of lines in the header
        count = 0
        #-- Reading over header text
        while kwargs['header'] and (count < file_lines):
            #-- file line at count
            line = file_contents[count]
            #-- find Time within line to set HEADER flag to False when found
            kwargs['header'] = not bool(re.search(r"Time",line))
            #-- add 1 to counter
            count += 1

        #-- catch to see if HEADER flag was not set to false
        if kwargs['header']:
            raise IOError('Data lines not found in file {0}'.format(geocenter_file))

        #-- number of months within the file
        n_mon = np.int64(file_lines - count)
        self.C10 = np.zeros((n_mon))
        self.C11 = np.zeros((n_mon))
        self.S11 = np.zeros((n_mon))
        self.time = np.zeros((n_mon))
        JD = np.zeros((n_mon))
        self.month = np.zeros((n_mon), dtype=np.int64)

        #-- compile numerical expression operator
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)

        #-- for every other line:
        for t, line in enumerate(file_contents[count:]):
            #-- find numerical instances in line including integers, exponents,
            #-- decimal points and negatives
            line_contents = rx.findall(line)

            #-- extacting time
            self.time[t]=np.float64(line_contents[0])
            #-- extracting spherical harmonics and convert to cmwe
            self.C10[t]=0.1*np.float64(line_contents[1])
            self.C11[t]=0.1*np.float64(line_contents[2])
            self.S11[t]=0.1*np.float64(line_contents[3])

            #-- calculate the GRACE months
            if (len(line_contents) == 5):
                #-- months are included as last column
                self.month[t] = np.int64(line_contents[4])
            else:
                #-- months to be calculated from date
                #-- Calculation of the Julian date from calendar date
                JD[t] = gravity_toolkit.time.calendar_to_julian(self.time[t])
                #-- convert the julian date into calendar dates (day, month, year)
                cal_date = gravity_toolkit.time.convert_julian(JD[t])
                #-- calculate the GRACE month (Apr02 == 004)
                #-- https://grace.jpl.nasa.gov/data/grace-months/
                #-- Notes on special months (e.g. 119, 120) below
                self.month[t] = gravity_toolkit.time.calendar_to_grace(
                    cal_date['year'], month=cal_date['month'])

        #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
        #-- Accelerometer shutoffs make the relation between month number
        #-- and date more complicated as days from other months are used
        #-- For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
        #-- For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
        #-- For all: May 2015 (161) is centered in Apr 2015 (160)
        self.month = gravity_toolkit.time.adjust_months(self.month)
        #-- converts from cm water equivalent to fully-normalized
        self.from_cmwe()
        #-- return the geocenter harmonics
        return self

    def from_tellus(self, geocenter_file, **kwargs):
        """
        Reads monthly geocenter spherical harmonic data files from GRACE Tellus
        Technical Notes (TN-13) calculated using GRACE/GRACE-FO measurements and
        Ocean Models of Degree 1

        `Datasets distributed by NASA PO.DAAC
        <https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1>`_

        Parameters
        ----------
        geocenter_file: str
            degree 1 file

                - ``CSR``: TN-13_GEOC_CSR_RL06.txt
                - ``GFZ``: TN-13_GEOC_GFZ_RL06.txt
                - ``JPL``: TN-13_GEOC_JPL_RL06.txt
        header: bool, default True
            file contains header text to be skipped
        JPL: bool, default True
            use JPL TN-13 geocenter files with self-attraction and loading

        References
        ----------
        .. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr,
            "Estimating geocenter variations from a combination
            of GRACE and ocean model output", *Journal of Geophysical
            Research*, 113(B08410), (2008).
            `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_

        .. [Sun2016a] Y. Sun, R. Riva, and P. Ditmar, "Observed changes
            in the Earth's dynamic oblateness from GRACE data and
            geophysical models", *Journal of Geodesy*, 90(1), 81-89, (2016).
            `doi: 10.1007/s00190-015-0852-y <https://doi.org/10.1007/s00190-015-0852-y>`_

        .. [Sun2016b] Y. Sun, R. Riva, and P. Ditmar, "Optimizing estimates of
            annual variations and trends in geocenter motion and J2 from
            a combination of GRACE data and geophysical models",
            *Journal of Geophysical Research: Solid Earth*, 121, (2016).
            `doi: 10.1002/2016JB013073 <https://doi.org/10.1002/2016JB013073>`_
        """
        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- set default keyword arguments
        kwargs.setdefault('header',True)
        kwargs.setdefault('JPL',True)

        #-- read degree 1 file and get contents
        with open(self.filename, mode='r', encoding='utf8') as f:
            file_contents = f.read().splitlines()
        #-- number of lines contained in the file
        file_lines = len(file_contents)

        #-- counts the number of lines in the header
        count = 0
        #-- Reading over header text
        header_flag = r"end\sof\sheader" if kwargs['JPL'] else r"'\(a6,"
        while kwargs['header']:
            #-- file line at count
            line = file_contents[count]
            #-- find header_flag within line to set HEADER flag to False when found
            kwargs['header'] = not bool(re.match(header_flag,line))
            #-- add 1 to counter
            count += 1

        #-- number of months within the file
        n_mon = (file_lines - count)//2
        #-- GRACE/GRACE-FO months
        self.month = np.zeros((n_mon),dtype=np.int64)
        #-- calendar dates in year-decimal
        self.time = np.zeros((n_mon))
        #-- spherical harmonic data
        self.C10 = np.zeros((n_mon))
        self.C11 = np.zeros((n_mon))
        self.S11 = np.zeros((n_mon))
        #-- spherical harmonic uncertainties
        self.eC10 = np.zeros((n_mon))
        self.eC11 = np.zeros((n_mon))
        self.eS11 = np.zeros((n_mon))

        #-- compile numerical expression operator
        regex_pattern = r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?'
        rx = re.compile(regex_pattern, re.VERBOSE)

        #-- time count
        t = 0
        #-- for every line of data
        for line in file_contents[count:]:
            #-- find numerical instances in line including integers, exponents,
            #-- decimal points and negatives
            line_contents = rx.findall(line)

            #-- spherical harmonic order
            m = np.int64(line_contents[2])
            #-- extract spherical harmonic data for order
            if (m == 0):
                self.C10[t] = np.float64(line_contents[3])
                self.eC10[t] = np.float64(line_contents[5])
            elif (m == 1):
                self.C11[t] = np.float64(line_contents[3])
                self.S11[t] = np.float64(line_contents[4])
                self.eC11[t] = np.float64(line_contents[5])
                self.eS11[t] = np.float64(line_contents[6])
            else:
                raise Exception('Unknown harmonic order {0:d}'.format(m))

            #-- calendar year and month
            if kwargs['JPL']:
                #-- start and end date of month
                start_date = time.strptime(line_contents[7][:8],r'%Y%m%d')
                end_date = time.strptime(line_contents[8][:8],r'%Y%m%d')
                #-- convert date to year decimal
                ts = gravity_toolkit.time.convert_calendar_decimal(start_date.tm_year,
                    start_date.tm_mon, day=start_date.tm_mday)
                te = gravity_toolkit.time.convert_calendar_decimal(end_date.tm_year,
                    end_date.tm_mon, day=end_date.tm_mday)
                #-- calculate mean time
                self.time[t] = np.mean([ts,te])
                #-- calculate year and month for estimating GRACE/GRACE-FO month
                year = np.floor(self.time[t])
                month = np.int64(12*(self.time[t] % 1) + 1)
            else:
                #-- dates of month
                cal_date = time.strptime(line_contents[0][:6],r'%Y%m')
                #-- calculate year and month for estimating GRACE/GRACE-FO month
                year = cal_date.tm_year
                month = cal_date.tm_mon
                #-- convert date to year decimal
                self.time[t], = gravity_toolkit.time.convert_calendar_decimal(
                    cal_date.tm_year, cal_date.tm_mon)
            #-- estimated GRACE/GRACE-FO month
            #-- Accelerometer shutoffs complicate the month number calculation
            self.month[t] = gravity_toolkit.time.calendar_to_grace(year,month)

            #-- will only advance in time after reading the
            #-- order 1 coefficients (t+0=t)
            t += m

        #-- The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
        #-- Accelerometer shutoffs make the relation between month number
        #-- and date more complicated as days from other months are used
        self.month = gravity_toolkit.time.adjust_months(self.month)
        #-- return the geocenter harmonics
        return self

    def from_netCDF4(self, geocenter_file, **kwargs):
        """
        Reads geocenter file and extracts dates and spherical harmonic data
        from a netCDF4 file

        Parameters
        ----------
        geocenter_file: str
            degree 1 netCDF4 file
        compression: str or NoneType, default None
            file compression type

        References
        ----------
        .. [Sutterley2019] T. C. Sutterley, and I. Velicogna, "Improved
            estimates of geocenter variability from time-variable gravity
            and ocean model outputs", *Remote Sensing*, 11(18), 2108, (2019).
            `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_
        """
        kwargs.setdefault('compression',None)
        #-- set filename
        self.case_insensitive_filename(geocenter_file)
        #-- Open the netCDF4 file for reading
        if (kwargs['compression'] == 'gzip'):
            #-- read gzipped file as in-memory (diskless) netCDF4 dataset
            with gzip.open(self.filename,'r') as f:
                fileID = netCDF4.Dataset(uuid.uuid4().hex,
                    memory=f.read())
        elif (kwargs['compression'] == 'bytes'):
            #-- read as in-memory (diskless) netCDF4 dataset
            fileID = netCDF4.Dataset(uuid.uuid4().hex,
                memory=self.filename.read())
        else:
            fileID = netCDF4.Dataset(self.filename, 'r')
        #-- Getting the data from each netCDF4 variable
        DEG1 = {}
        #-- converting netCDF4 objects into numpy arrays
        for key,val in fileID.variables.items():
            DEG1[key] = val[:].copy()
        #-- close the netCDF4 file
        fileID.close()
        #-- return the geocenter harmonics
        return self.from_dict(DEG1)

    def copy(self, **kwargs):
        """
        Copy a geocenter object to a new geocenter object

        Parameters
        ----------
            fields: list
                default keys in geocenter object
        """
        #-- set default keyword arguments
        kwargs.setdefault('fields',['time','month',
            'C10','C11','S11','eC10','eC11','eS11',
            'X','Y','Z'])
        temp = geocenter()
        #-- try to assign variables to self
        for key in kwargs['fields']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.copy(val))
            except AttributeError:
                pass
        return temp

    def from_dict(self, temp, **kwargs):
        """
        Convert a dictionary object to a geocenter object

        Parameters
        ----------
        temp: obj
            dictionary object to be converted
        fields: list
            default keys in dictionary
        """
        #-- set default keyword arguments
        kwargs.setdefault('fields',['time','month',
            'C10','C11','S11','eC10','eC11','eS11',
            'X','Y','Z','X_sigma','Y_sigma','Z_sigma'])
        #-- assign dictionary variables to self
        for key in kwargs['fields']:
            try:
                setattr(self, key, temp[key].copy())
            except (AttributeError, KeyError):
                pass
        return self

    def from_harmonics(self, temp, **kwargs):
        """
        Convert a harmonics object to a geocenter object

        Parameters
        ----------
        temp: obj
            harmonics object to be converted
        fields: list
            default keys in harmonics object
        """
        #-- reassign shape and ndim attributes
        temp.update_dimensions()
        #-- set default keyword arguments
        kwargs.setdefault('fields',['time','month','filename'])
        #-- try to assign variables to self
        for key in kwargs['fields']:
            try:
                val = getattr(temp, key)
                setattr(self, key, np.copy(val))
            except AttributeError:
                pass
        #-- get spherical harmonic objects
        if (temp.ndim == 2):
            self.C10 = np.copy(temp.clm[1,0])
            self.C11 = np.copy(temp.clm[1,1])
            self.S11 = np.copy(temp.slm[1,1])
        elif (temp.ndim == 3):
            self.C10 = np.copy(temp.clm[1,0,:])
            self.C11 = np.copy(temp.clm[1,1,:])
            self.S11 = np.copy(temp.slm[1,1,:])
        #-- return the geocenter object
        return self

    def from_matrix(self, clm, slm):
        """
        Converts spherical harmonic matrices to a geocenter object

        Parameters
        ----------
        clm: float
            cosine spherical harmonics of degree 1
        slm: float
            sine spherical harmonics of degree 1
        """
        #-- verify dimensions
        clm = np.atleast_3d(clm)
        slm = np.atleast_3d(slm)
        #-- output geocenter object
        self.C10 = np.copy(clm[1,0,:])
        self.C11 = np.copy(clm[1,1,:])
        self.S11 = np.copy(slm[1,1,:])
        return self

    def to_dict(self, **kwargs):
        """
        Convert a geocenter object to a dictionary object

        Parameters
        ----------
        fields: obj
            default attributes in geocenter object
        """
        #-- output dictionary
        temp = {}
        #-- set default keyword arguments
        kwargs.setdefault('fields',['time','month',
            'C10','C11','S11','eC10','eC11','eS11',
            'X','Y','Z','X_sigma','Y_sigma','Z_sigma'])
        #-- assign dictionary variables to self
        for key in kwargs['fields']:
            try:
                val = getattr(self, key)
            except (AttributeError, KeyError):
                pass
            else:
                temp[key] = copy.copy(val)
        #-- return the dictionary object
        return temp

    def to_matrix(self):
        """
        Converts a geocenter object to spherical harmonic matrices
        """
        #-- verify dimensions
        _,nt = np.shape(np.atleast_2d(self.C10))
        #-- output spherical harmonics
        clm = np.zeros((2,2,nt))
        slm = np.zeros((2,2,nt))
        #-- copy geocenter harmonics to matrices
        clm[1,0,:] = np.atleast_2d(self.C10)
        clm[1,1,:] = np.atleast_2d(self.C11)
        slm[1,1,:] = np.atleast_2d(self.S11)
        return dict(clm=clm, slm=slm)

    def to_cartesian(self, kl=0.0):
        """
        Converts normalized spherical harmonics to cartesian geocenter variations

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        #-- Stokes Coefficients to cartesian geocenter
        try:
            self.Z = self.C10*self.radius*np.sqrt(3.0)/(1.0 + kl)
            self.X = self.C11*self.radius*np.sqrt(3.0)/(1.0 + kl)
            self.Y = self.S11*self.radius*np.sqrt(3.0)/(1.0 + kl)
        except Exception as e:
            pass
        #-- convert errors to cartesian geocenter
        try:
            self.Z_sigma = self.eC10*self.radius*np.sqrt(3.0)/(1.0 + kl)
            self.X_sigma = self.eC11*self.radius*np.sqrt(3.0)/(1.0 + kl)
            self.Y_sigma = self.eS11*self.radius*np.sqrt(3.0)/(1.0 + kl)
        except Exception as e:
            pass
        return self

    def to_cmwe(self, kl=0.0):
        """
        Converts normalized spherical harmonics to centimeters water equivalent

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        #-- Average Density of the Earth [g/cm^3]
        rho_e = 5.517
        #-- Average Radius of the Earth [cm]
        rad_e = 6.371e8
        #-- convert to centimeters water equivalent
        self.C10 *= (rho_e*rad_e)/(1.0 + kl)
        self.C11 *= (rho_e*rad_e)/(1.0 + kl)
        self.S11 *= (rho_e*rad_e)/(1.0 + kl)
        #-- convert errors to centimeters water equivalent
        try:
            self.eC10 *= (rho_e*rad_e)/(1.0 + kl)
            self.eC11 *= (rho_e*rad_e)/(1.0 + kl)
            self.eS11 *= (rho_e*rad_e)/(1.0 + kl)
        except Exception as e:
            pass
        return self

    def to_mmwe(self, kl=0.0):
        """
        Converts normalized spherical harmonics to millimeters water equivalent

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        self.to_cmwe(kl=kl)
        #-- convert to millimeters water equivalent
        self.C10 *= 10.0
        self.C11 *= 10.0
        self.S11 *= 10.0
        #-- convert errors to millimeters water equivalent
        try:
            self.eC10 *= 10.0
            self.eC11 *= 10.0
            self.eS11 *= 10.0
        except Exception as e:
            pass
        return self

    def from_cartesian(self, kl=0.0):
        """
        Converts cartesian geocenter variations to normalized spherical harmonics

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        #-- cartesian geocenter to Stokes Coefficients
        self.C10 = (1.0 + kl)*self.Z/(self.radius*np.sqrt(3.0))
        self.C11 = (1.0 + kl)*self.X/(self.radius*np.sqrt(3.0))
        self.S11 = (1.0 + kl)*self.Y/(self.radius*np.sqrt(3.0))
        #-- convert cartesian geocenter to stokes coefficients
        try:
            self.eC10 = (1.0 + kl)*self.Z_sigma/(self.radius*np.sqrt(3.0))
            self.eC11 = (1.0 + kl)*self.X_sigma/(self.radius*np.sqrt(3.0))
            self.eS11 = (1.0 + kl)*self.Y_sigma/(self.radius*np.sqrt(3.0))
        except Exception as e:
            pass
        return self

    def from_cmwe(self, kl=0.0):
        """
        Normalizes spherical harmonics from centimeters water equivalent (cmwe)

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        #-- Average Density of the Earth [g/cm^3]
        rho_e = 5.517
        #-- Average Radius of the Earth [cm]
        rad_e = 6.371e8
        #-- convert from centimeters water equivalent
        self.C10 *= (1.0 + kl)/(rho_e*rad_e)
        self.C11 *= (1.0 + kl)/(rho_e*rad_e)
        self.S11 *= (1.0 + kl)/(rho_e*rad_e)
        #-- convert errors from centimeters water equivalent
        try:
            self.eC10 *= (1.0 + kl)/(rho_e*rad_e)
            self.eC11 *= (1.0 + kl)/(rho_e*rad_e)
            self.eS11 *= (1.0 + kl)/(rho_e*rad_e)
        except Exception as e:
            pass
        return self

    def from_mmwe(self, kl=0.0):
        """
        Normalizes spherical harmonics from millimeters water equivalent (mmwe)

        Parameters
        ----------
        kl: float
            gravitational load love number of degree 1
        """
        self.from_cmwe(kl=kl)
        #-- convert from millimeters water equivalent
        self.C10 /= 10.0
        self.C11 /= 10.0
        self.S11 /= 10.0
        #-- convert errors from centimeters water equivalent
        try:
            self.eC10 /= 10.0
            self.eC11 /= 10.0
            self.eS11 /= 10.0
        except Exception as e:
            pass
        return self

    def mean(self, apply=False, indices=Ellipsis):
        """
        Compute mean gravitational field and remove from data if specified

        Parameters
        ----------
        apply: bool, default False
            remove the mean field from the input harmonics
        indices: int, default Ellipsis
            indices of input harmonics object to compute mean
        """
        temp = geocenter()
        #-- calculate mean static field
        temp.C10 = np.mean(self.C10[indices])
        temp.C11 = np.mean(self.C11[indices])
        temp.S11 = np.mean(self.S11[indices])
        #-- calculating the time-variable gravity field by removing
        #-- the static component of the gravitational field
        if apply:
            self.C10 -= temp.C10
            self.C11 -= temp.C11
            self.S11 -= temp.S11
        #-- calculate mean of temporal variables
        for key in ['time','month']:
            try:
                val = getattr(self, key)
                setattr(temp, key, np.mean(val[indices]))
            except:
                continue
        #-- return the mean field
        return temp

    def add(self, temp):
        """
        Add two geocenter objects

        Parameters
        ----------
        temp: obj
            geocenter object to be added
        """
        self.C10 += temp.C10
        self.C11 += temp.C11
        self.S11 += temp.S11
        return self

    def subtract(self, temp):
        """
        Subtract one geocenter object from another

        Parameters
        ----------
        temp: obj
            geocenter object to be subtracted
        """
        self.C10 -= temp.C10
        self.C11 -= temp.C11
        self.S11 -= temp.S11
        return self

    def multiply(self, temp):
        """
        Multiply two geocenter objects

        Parameters
        ----------
        temp: obj
            geocenter object to be multiplied
        """
        self.C10 *= temp.C10
        self.C11 *= temp.C11
        self.S11 *= temp.S11
        return self

    def divide(self, temp):
        """
        Divide one geocenter object from another

        Parameters
        ----------
        temp: obj
            geocenter object to be divided
        """
        self.C10 /= temp.C10
        self.C11 /= temp.C11
        self.S11 /= temp.S11
        return self

    def scale(self, var):
        """
        Multiply a geocenter object by a constant

        Parameters
        ----------
        var: float
            scalar value to which the geocenter object will be multiplied
        """
        temp = geocenter()
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        #-- multiply by a single constant or a time-variable scalar
        temp.C10 = var*self.C10
        temp.C11 = var*self.C11
        temp.S11 = var*self.S11
        return temp

    def power(self, power):
        """
        Raise a geocenter object to a power

        Parameters
        ----------
        power: float
            power to which the geocenter object will be raised
        """
        temp = geocenter()
        temp.time = np.copy(self.time)
        temp.month = np.copy(self.month)
        temp.C10 = np.power(self.C10,power)
        temp.C11 = np.power(self.C11,power)
        temp.S11 = np.power(self.S11,power)
        return temp
