#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (11/2022)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added more time parsing for longer periods
    Updated 08/2022: added file parsing functions from GRACE date utilities
        added function to dynamically select newest version of granules
        output variables to unit conversion to seconds and the number of days
        per month for both leap and standard years
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2021: added function for calendar year (decimal) to Julian Day
    Updated 09/2021: add functions for converting to and from GRACE months
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: added adjust_months function to fix special months cases
    Updated 01/2021: add date parser for cases when only a date and no units
    Updated 12/2020: merged with convert_julian and convert_calendar_decimal
        added calendar_days routine to get number of days per month
    Updated 09/2020: parse date strings "time-units since yyyy-mm-dd hh:mm:ss"
    Updated 08/2020: added NASA Earthdata routines for downloading from CDDIS
    Written 07/2020
"""
import os
import re
import copy
import warnings
import datetime
import numpy as np
import dateutil.parser

# conversion factors between time units and seconds
_to_sec = {'microseconds': 1e-6, 'microsecond': 1e-6,
           'microsec': 1e-6, 'microsecs': 1e-6,
           'milliseconds': 1e-3, 'millisecond': 1e-3,
           'millisec': 1e-3, 'millisecs': 1e-3,
           'msec': 1e-3, 'msecs': 1e-3, 'ms': 1e-3,
           'seconds': 1.0, 'second': 1.0, 'sec': 1.0,
           'secs': 1.0, 's': 1.0,
           'minutes': 60.0, 'minute': 60.0,
           'min': 60.0, 'mins': 60.0,
           'hours': 3600.0, 'hour': 3600.0,
           'hr': 3600.0, 'hrs': 3600.0, 'h': 3600.0,
           'day': 86400.0, 'days': 86400.0, 'd': 86400.0}
# approximate conversions for longer periods
_to_sec['mon'] = 30.0 * 86400.0
_to_sec['month'] = 30.0 * 86400.0
_to_sec['months'] = 30.0 * 86400.0
_to_sec['common_year'] = 365.0 * 86400.0
_to_sec['common_years'] = 365.0 * 86400.0
_to_sec['year'] = 365.25 * 86400.0
_to_sec['years'] = 365.25 * 86400.0

# PURPOSE: parse a date string into epoch and units scale
def parse_date_string(date_string):
    """
    parse a date string of the form

    - time-units since ``yyyy-mm-dd hh:mm:ss``
    - ``yyyy-mm-dd hh:mm:ss`` for exact calendar dates

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss

    Returns
    -------
    epoch: list
        epoch of delta time
    conversion_factor: float
        multiplication factor to convert to seconds
    """
    # try parsing the original date string as a date
    try:
        epoch = dateutil.parser.parse(date_string)
    except ValueError:
        pass
    else:
        # return the epoch (as list)
        return (datetime_to_list(epoch),0.0)
    # split the date string into units and epoch
    units, epoch = split_date_string(date_string)
    if units not in _to_sec.keys():
        raise ValueError(f'Invalid units: {units}')
    # return the epoch (as list) and the time unit conversion factors
    return (datetime_to_list(epoch), _to_sec[units])

# PURPOSE: split a date string into units and epoch
def split_date_string(date_string):
    """
    split a date string into units and epoch

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss
    """
    try:
        units,_,epoch = date_string.split(None,2)
    except ValueError:
        raise ValueError(f'Invalid format: {date_string}')
    else:
        return (units.lower(),dateutil.parser.parse(epoch))

# PURPOSE: convert a datetime object into a list
def datetime_to_list(date):
    """
    convert a datetime object into a list

    Parameters
    ----------
    date: datetime object

    Returns
    -------
    date: list
        [year,month,day,hour,minute,second]
    """
    return [date.year,date.month,date.day,date.hour,date.minute,date.second]

# PURPOSE: extract parameters from filename
def parse_grace_file(granule):
    """
    Extract dates from GRACE/GRACE-FO files

    Parameters
    ----------
    granule: str
        GRACE/GRACE-FO Level-2 spherical harmonic data file
    """
    # verify that filename is reduced to basename
    file_basename = os.path.basename(granule)
    # compile numerical expression operator for parameters from files
    # UTCSR: The University of Texas at Austin Center for Space Research
    # EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
    # GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
    # JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
    # JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    # GRGS: French Centre National D'Etudes Spatiales (CNES)
    # COSTG: International Combined Time-variable Gravity Fields
    args = r'UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC|GRGS|COSTG'
    regex_pattern = (r'(.*?)-2_(\d{{4}})(\d{{3}})-(\d{{4}})(\d{{3}})_'
        r'(.*?)_({0})_(.*?)_(\d+)(.*?)(\.gz|\.gfc)?$').format(args)
    rx = re.compile(regex_pattern, re.VERBOSE)
    # extract parameters from input filename
    PFX,SY,SD,EY,ED,AUX,PRC,F1,DRL,F2,SFX = rx.findall(file_basename).pop()
    # return the start and end date lists
    return ((SY,SD),(EY,ED))

# PURPOSE: extract dates from GRAZ or Swarm files with regular expressions
def parse_gfc_file(granule, PROC, DSET):
    """
    Extract dates from Gravity Field Coefficient (gfc) files

    Parameters
    ----------
    granule: str
        GRAZ or Swarm spherical harmonic data file
    PROC: str
        GRACE/GRACE-FO Processing Center or Satellite mission

            - ``'GRAZ'``: Institute of Geodesy from GRAZ University of Technology
            - ``'Swarm'``: Time-variable gravity data from Swarm satellites
    DSET: str
        GRACE/GRACE-FO/Swarm dataset

            - ``'GAA'``: non-tidal atmospheric correction
            - ``'GAB'``: non-tidal oceanic correction
            - ``'GAC'``: combined non-tidal atmospheric and oceanic correction
            - ``'GAD'``: ocean bottom pressure product
            - ``'GSM'``: corrected monthly static gravity field product
    """
    # verify that filename is reduced to basename
    file_basename = os.path.basename(granule)
    # extract parameters from input filename
    if (PROC == 'GRAZ'):
        # regular expression operators for ITSG data and models
        itsg_products = []
        itsg_products.append(r'atmosphere')
        itsg_products.append(r'dealiasing')
        itsg_products.append(r'oceanBottomPressure')
        itsg_products.append(r'ocean')
        itsg_products.append(r'Grace2014')
        itsg_products.append(r'Grace2016')
        itsg_products.append(r'Grace2018')
        itsg_products.append(r'Grace_operational')
        regex_pattern=(r'(AOD1B_RL\d+|model|ITSG)[-_]({0})(_n\d+)?_'
            r'(\d+)-(\d+)(\.gfc)').format(r'|'.join(itsg_products))
        # compile regular expression operator for parameters from files
        rx = re.compile(regex_pattern, re.VERBOSE | re.IGNORECASE)
        # extract parameters from input filename
        PFX,PRD,trunc,year,month,SFX = rx.findall(file_basename).pop()
        # number of days in each month for the calendar year
        dpm = calendar_days(int(year))
        # create start and end date lists
        start_date = [int(year),int(month),1,0,0,0]
        end_date = [int(year),int(month),dpm[int(month)-1],23,59,59]
    elif (PROC == 'Swarm') and (DSET == 'GSM'):
        # regular expression operators for Swarm data
        regex_pattern=r'(SW)_(.*?)_(EGF_SHA_2)__(.*?)_(.*?)_(.*?)(\.gfc|\.ZIP)'
        # compile regular expression operator for parameters from files
        rx = re.compile(regex_pattern, re.VERBOSE | re.IGNORECASE)
        # extract parameters from input filename
        SAT,tmp,PROD,starttime,endtime,RL,SFX = rx.findall(file_basename).pop()
        start_date,_ = parse_date_string(starttime)
        end_date,_ = parse_date_string(endtime)
    elif (PROC == 'Swarm') and (DSET != 'GSM'):
        # regular expression operators for Swarm models
        regex_pattern=(r'(GAA|GAB|GAC|GAD)_Swarm_(\d+)_(\d{2})_(\d{4})'
            r'(\.gfc|\.ZIP)')
        # compile regular expression operator for parameters from files
        rx = re.compile(regex_pattern, re.VERBOSE | re.IGNORECASE)
        # extract parameters from input filename
        PROD,trunc,month,year,SFX = rx.findall(file_basename).pop()
        # number of days in each month for the calendar year
        dpm = calendar_days(int(year))
        # create start and end date lists
        start_date = [int(year),int(month),1,0,0,0]
        end_date = [int(year),int(month),dpm[int(month)-1],23,59,59]
    # return the start and end date lists
    return (start_date, end_date)

def reduce_by_date(granules):
    """
    Reduce list of GRACE/GRACE-FO files by date to the newest version

    Parameters
    ----------
    granules: list
        GRACE/GRACE-FO Level-2 spherical harmonic data files
    """
    # list of dates for all input files
    date_list = [parse_grace_file(f) for f in granules]
    unique_list = []
    # compile numerical expression operator for parameters from files
    # UTCSR: The University of Texas at Austin Center for Space Research
    # EIGEN: GFZ German Research Center for Geosciences (RL01-RL05)
    # GFZOP: GFZ German Research Center for Geosciences (RL06+GRACE-FO)
    # JPLEM: NASA Jet Propulsion Laboratory (harmonic solutions)
    # JPLMSC: NASA Jet Propulsion Laboratory (mascon solutions)
    # GRGS: French Centre National D'Etudes Spatiales (CNES)
    # COSTG: International Combined Time-variable Gravity Fields
    args = r'UTCSR|EIGEN|GFZOP|JPLEM|JPLMSC|GRGS|COSTG'
    regex_pattern = (r'(.*?)-2_(\d{{4}})(\d{{3}})-(\d{{4}})(\d{{3}})_(.*?)_'
        r'({0})_(.*?)_(\d{{2}})(\d{{2}})(.*?)(\.gz|\.gfc)?$').format(args)
    rx = re.compile(regex_pattern, re.VERBOSE)
    # for each unique date
    for d in sorted(set(date_list)):
        if (date_list.count(d) == 1):
            i = date_list.index(d)
            unique_list.append(granules[i])
        elif (date_list.count(d) >= 2):
            # if more than 1 file with date use newest version
            indices = [i for i, dt in enumerate(date_list) if (dt == d)]
            # find each version within the file
            versions = []
            for i in indices:
                # verify that filename is reduced to basename
                file_basename = os.path.basename(granules[i])
                # parse filename to get file version
                PFX,SY,SD,EY,ED,AUX,PRC,F1,DRL,VER,F2,SFX = \
                    rx.findall(file_basename).pop()
                # append to list of file versions
                versions.append(int(VER))
            # find file with newest version
            i = versions.index(max(versions))
            unique_list.append(granules[indices[i]])
    # return the sorted list of files with unique dates
    return unique_list

# PURPOSE: Adjust GRACE/GRACE-FO months to fix "Special Cases"
def adjust_months(grace_month):
    """
    Adjust estimated GRACE/GRACE-FO months to fix "Special Cases"

    Parameters
    ----------
    grace_month: int
        GRACE/GRACE-FO months

    Notes
    -----
    The "Special Months" (Nov 2011, Dec 2011 and April 2012) with
    Accelerometer shutoffs make the relation between month number
    and date more complicated as days from other months are used.

    For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)

    For JPL: Dec 2011 (120) is centered in Jan 2012 (121)

    For all: May 2015 (161) is centered in Apr 2015 (160)

    For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)
    """
    # verify dimensions
    grace_month = np.atleast_1d(grace_month)
    # number of months
    nmon = len(grace_month)
    # create temporary months object
    m = np.zeros_like(grace_month)
    # find unique months
    _,i,c = np.unique(grace_month,return_inverse=True,return_counts=True)
    # simple unique months case
    case1, = np.nonzero(c[i] == 1)
    m[case1] = grace_month[case1]
    # Special Months cases
    case2, = np.nonzero(c[i] == 2)
    # for each special case month
    for j in case2:
        # prior month, current month, subsequent 2 months
        mm1 = grace_month[j-1]
        mon = grace_month[j]
        mp1 = grace_month[j+1] if (j < (nmon-1)) else (mon + 1)
        mp2 = grace_month[j+2] if (j < (nmon-2)) else (mp1 + 1)
        # determine the months which meet the criteria need to be adjusted
        if (mon == (mm1 + 1)):
            # case where month is correct
            # but subsequent month needs to be +1
            m[j] = np.copy(grace_month[j])
        elif (mon == mm1) and (mon != m[j-1]):
            # case where prior month needed to be -1
            # but current month is correct
            m[j] = np.copy(grace_month[j])
        elif (mon == mm1):
            # case where month should be +1
            m[j] = grace_month[j] + 1
        elif (mon == mp1) and ((mon == (mm1 + 2)) or (mp2 == (mp1 + 1))):
            # case where month should be -1
            m[j] = grace_month[j] - 1
    # update months and remove singleton dimensions if necessary
    return np.squeeze(m)

# PURPOSE: convert calendar dates to GRACE/GRACE-FO months
def calendar_to_grace(year,month=1,around=np.floor):
    """
    Converts calendar dates to GRACE/GRACE-FO months

    Parameters
    ----------
    year: float
        calendar year
    month: int, default 1
        calendar month
    around: obj, default np.floor
        method of rounding to nearest method

    Returns
    -------
    grace_month: int
        GRACE/GRACE-FO month
    """
    grace_month = around(12.0*(year - 2002.0)) + month
    return np.array(grace_month, dtype=int)

# PURPOSE: convert GRACE/GRACE-FO months to calendar dates
def grace_to_calendar(grace_month):
    """
    Converts GRACE/GRACE-FO months to calendar dates

    Parameters
    ----------
    grace_month: int
        GRACE/GRACE-FO month

    Returns
    -------
    year: int
        calendar year
    month: int
        calendar month
    """
    year = np.array(2002 + (grace_month-1)//12).astype(int)
    month = np.mod(grace_month-1,12) + 1
    return (year, month)

# PURPOSE: convert calendar dates to Julian days
def calendar_to_julian(year_decimal):
    """
    Converts calendar dates to Julian days

    Parameters
    ----------
    year: float
        calendar year

    Returns
    -------
    JD: float
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    """
    # calculate year
    year = np.floor(year_decimal)
    # calculation of day of the year
    dpy = calendar_days(year).sum()
    DofY = dpy*(year_decimal % 1)
    # Calculation of the Julian date from year and DofY
    JD = np.array(367.0*year - np.floor(7.0*year/4.0) -
        np.floor(3.0*(np.floor((7.0*year - 1.0)/700.0) + 1.0)/4.0) +
        DofY + 1721058.5, dtype=np.float64)
    return JD

# days per month in a leap and a standard year
# only difference is February (29 vs. 28)
_dpm_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_dpm_stnd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# PURPOSE: gets the number of days per month for a given year
def calendar_days(year):
    """
    Calculates the number of days per month for a given year

    Parameters
    ----------
    year: int
        calendar year

    Returns
    -------
    dpm: float
        number of days for each month
    """
    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    # find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return np.array(_dpm_leap, dtype=np.float64)
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return np.array(_dpm_stnd, dtype=np.float64)

# PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0):
    """
    Convert delta time from seconds since epoch1 to time since epoch2

    Parameters
    ----------
    delta_time: float
        seconds since epoch1
    epoch1: tuple or NoneType, default None
        epoch for input delta_time
    epoch2: tuple or NoneType, default None
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units
    """
    epoch1 = datetime.datetime(*epoch1)
    epoch2 = datetime.datetime(*epoch2)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    # subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

# PURPOSE: calculate the delta time from calendar date
# http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0,
    epoch=(1992,1,1,0,0,0), scale=1.0):
    """
    Calculate the time in time units since epoch from calendar dates

    Parameters
    ----------
    year: float
        calendar year
    month: float
        month of the year
    day: float
        day of the month
    hour: float, default 0.0
        hour of the day
    minute: float, default 0.0
        minute of the hour
    second: float, default 0.0
        second of the minute
    epoch: tuple, default (1992,1,1,0,0,0)
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units

    Returns
    -------
    delta_time: float
        days since epoch
    """
    # calculate date in Modified Julian Days (MJD) from calendar date
    # MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    epoch1 = datetime.datetime(1858,11,17,0,0,0)
    epoch2 = datetime.datetime(*epoch)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    # return the date in days since epoch (or scaled to units)
    return scale*np.array(MJD - delta_time_epochs/86400.0,dtype=np.float64)

# PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(year, month, day=None, hour=None, minute=None,
    second=None, DofY=None):
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Parameters
    ----------
    year: float
        calendar year
    month: float
        calendar month
    day: float or NoneType, default None
        day of the month
    hour: float or NoneType, default None
        hour of the day
    minute: float or NoneType, default None
        minute of the hour
    second: float or NoneType, default None
        second of the minute
    DofY: float or NoneType, default None
        day of the year (January 1 = 1)

    Returns
    -------
    t_date: float
        date in decimal-year format

    References
    ----------
    .. [Dershowitz2008] Dershowitz, N. and E.M. Reingold.
        *Calendrical Calculations*, (2008).
        Cambridge: Cambridge University Press.
    """

    # number of dates
    n_dates = len(np.atleast_1d(year))

    # create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    # day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    # remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    # create output date variable
    t_date = np.zeros((n_dates))

    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array(_dpm_leap, dtype=np.float64)
    dpm_stnd = np.array(_dpm_stnd, dtype=np.float64)

    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    # find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    # calculate the day of the year
    if DofY is not None:
        # if entered directly as an input
        # remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        # use calendar month and day of the month to calculate day of the year
        # month minus 1: January = 0, February = 1, etc (indice of month)
        # in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=np.int64) - 1

        # day of month
        if day is not None:
            # remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            # if not entering days as an input
            # will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        # create matrix with the lower half = 1
        # this matrix will be used in a matrix multiplication
        # to calculate the total number of days for prior months
        # the -1 will make the diagonal == 0
        # i.e. first row == all zeros and the
        # last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        # using a dot product to calculate total number of days
        # for the months before the input date
        # basically is sum(i*dpm)
        # where i is 1 for all months < the month of interest
        # and i is 0 for all months >= the month of interest
        # month of interest is zero as the exact days will be
        # used to calculate the date

        # calculate the day of the year for leap and standard
        # use total days of all months before date
        # and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    # hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    # minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    # second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    # calculate decimal date
    # convert hours, minutes and seconds into days
    # convert calculated fractional days into decimal fractions of the year
    # Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    # Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

# PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD, **kwargs):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: float
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    astype: str or NoneType, default None
        convert output to variable type
    format: str, default 'dict'
        format of output variables

            - ``'dict'``: dictionary with variable keys
            - ``'tuple'``: tuple in most-to-least-significant order
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    year: float
        calendar year
    month: float
        calendar month
    day: float
        day of the month
    hour: float
        hour of the day
    minute: float
        minute of the hour
    second: float
        second of the minute

    References
    ----------
    .. [Press1988] *Numerical Recipes in C*, William H. Press,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Second Edition, Cambridge University Press, (1988).
    .. [Hatcher1984] Hatcher, D. A., "Simple Formulae for Julian Day Numbers and
        Calendar Dates", *Quarterly Journal of the Royal Astronomical
        Society*, 25(1), (1984).
    """
    # set default keyword arguments
    kwargs.setdefault('astype', None)
    kwargs.setdefault('format', 'dict')
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(ASTYPE='astype', FORMAT='format')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        single_value = True
    else:
        single_value = False

    # verify julian day
    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    # calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    # calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    # calculate day, month, year and hour
    day = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    month = F - 1.0 - 12.0*np.floor(F/14.0)
    year = D - 4715.0 - np.floor((7.0 + month)/10.0)
    hour = np.floor(24.0*(JD + 0.5 - JDO))
    # calculate minute and second
    G = (JD + 0.5 - JDO) - hour/24.0
    minute = np.floor(G*1440.0)
    second = (G - minute/1440.0) * 86400.0

    # convert all variables to output type (from float)
    if kwargs['astype'] is not None:
        year = year.astype(kwargs['astype'])
        month = month.astype(kwargs['astype'])
        day = day.astype(kwargs['astype'])
        hour = hour.astype(kwargs['astype'])
        minute = minute.astype(kwargs['astype'])
        second = second.astype(kwargs['astype'])

    # if only a single value was imported initially: remove singleton dims
    if single_value:
        year = year.item(0)
        month = month.item(0)
        day = day.item(0)
        hour = hour.item(0)
        minute = minute.item(0)
        second = second.item(0)

    # return date variables in output format
    if (kwargs['format'] == 'dict'):
        return dict(year=year, month=month, day=day,
            hour=hour, minute=minute, second=second)
    elif (kwargs['format'] == 'tuple'):
        return (year, month, day, hour, minute, second)
    elif (kwargs['format'] == 'zip'):
        return zip(year, month, day, hour, minute, second)
