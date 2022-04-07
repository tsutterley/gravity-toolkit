#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (04/2022)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

UPDATE HISTORY:
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
import datetime
import numpy as np
import dateutil.parser

#-- PURPOSE: parse a date string into epoch and units scale
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
    #-- try parsing the original date string as a date
    try:
        epoch = dateutil.parser.parse(date_string)
    except ValueError:
        pass
    else:
        #-- return the epoch (as list)
        return (datetime_to_list(epoch),0.0)
    #-- split the date string into units and epoch
    units,epoch = split_date_string(date_string)
    conversion_factors = {'microseconds': 1e-6,'microsecond': 1e-6,
        'microsec': 1e-6,'microsecs': 1e-6,
        'milliseconds': 1e-3,'millisecond': 1e-3,'millisec': 1e-3,
        'millisecs': 1e-3,'msec': 1e-3,'msecs': 1e-3,'ms': 1e-3,
        'seconds': 1.0,'second': 1.0,'sec': 1.0,'secs': 1.0,'s': 1.0,
        'minutes': 60.0,'minute': 60.0,'min': 60.0,'mins': 60.0,
        'hours': 3600.0,'hour': 3600.0,'hr': 3600.0,
        'hrs': 3600.0,'h': 3600.0,
        'day': 86400.0,'days': 86400.0,'d': 86400.0}
    if units not in conversion_factors.keys():
        raise ValueError('Invalid units: {0}'.format(units))
    #-- return the epoch (as list) and the time unit conversion factors
    return (datetime_to_list(epoch),conversion_factors[units])

#-- PURPOSE: split a date string into units and epoch
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
        raise ValueError('Invalid format: {0}'.format(date_string))
    else:
        return (units.lower(),dateutil.parser.parse(epoch))

#-- PURPOSE: convert a datetime object into a list
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

#-- PURPOSE: Adjust GRACE/GRACE-FO months to fix "Special Cases"
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
    #-- verify dimensions
    grace_month = np.atleast_1d(grace_month)
    #-- number of months
    nmon = len(grace_month)
    #-- create temporary months object
    m = np.zeros_like(grace_month)
    #-- find unique months
    _,i,c = np.unique(grace_month,return_inverse=True,return_counts=True)
    #-- simple unique months case
    case1, = np.nonzero(c[i] == 1)
    m[case1] = grace_month[case1]
    #-- Special Months cases
    case2, = np.nonzero(c[i] == 2)
    #-- for each special case month
    for j in case2:
        # prior month, current month, subsequent 2 months
        mm1 = grace_month[j-1]
        mon = grace_month[j]
        mp1 = grace_month[j+1] if (j < (nmon-1)) else (mon + 1)
        mp2 = grace_month[j+2] if (j < (nmon-2)) else (mp1 + 1)
        #-- determine the months which meet the criteria need to be adjusted
        if (mon == (mm1 + 1)):
            #-- case where month is correct
            #-- but subsequent month needs to be +1
            m[j] = np.copy(grace_month[j])
        elif (mon == mm1) and (mon != m[j-1]):
            #-- case where prior month needed to be -1
            #-- but current month is correct
            m[j] = np.copy(grace_month[j])
        elif (mon == mm1):
            #-- case where month should be +1
            m[j] = grace_month[j] + 1
        elif (mon == mp1) and ((mon == (mm1 + 2)) or (mp2 == (mp1 + 1))):
            #-- case where month should be -1
            m[j] = grace_month[j] - 1
    #-- update months and remove singleton dimensions if necessary
    return np.squeeze(m)

#-- PURPOSE: convert calendar dates to GRACE/GRACE-FO months
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
    return np.array(grace_month,dtype=int)

#-- PURPOSE: convert GRACE/GRACE-FO months to calendar dates
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

#-- PURPOSE: convert calendar dates to Julian days
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
    #-- calculate year
    year = np.floor(year_decimal)
    #-- calculation of day of the year
    dpy = calendar_days(year).sum()
    DofY = dpy*(year_decimal % 1)
    #-- Calculation of the Julian date from year and DofY
    JD = np.array(367.0*year - np.floor(7.0*year/4.0) -
        np.floor(3.0*(np.floor((7.0*year - 1.0)/700.0) + 1.0)/4.0) +
        DofY + 1721058.5, dtype=np.float64)
    return JD

#-- PURPOSE: gets the number of days per month for a given year
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
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31],dtype=np.float64)
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float64)
    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    #-- find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return dpm_leap
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return dpm_stnd

#-- PURPOSE: convert times from seconds since epoch1 to time since epoch2
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
    #-- subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

#-- PURPOSE: calculate the delta time from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
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
    #-- calculate date in Modified Julian Days (MJD) from calendar date
    #-- MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    epoch1 = datetime.datetime(1858,11,17,0,0,0)
    epoch2 = datetime.datetime(*epoch)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- return the date in days since epoch (or scaled to units)
    return scale*np.array(MJD - delta_time_epochs/86400.0,dtype=np.float64)

#-- PURPOSE: Converts from calendar dates into decimal years
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

    #-- number of dates
    n_dates = len(np.atleast_1d(year))

    #-- create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    #-- day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    #-- remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    #-- create output date variable
    t_date = np.zeros((n_dates))

    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap=np.array([31,29,31,30,31,30,31,31,30,31,30,31], dtype=np.float64)
    dpm_stnd=np.array([31,28,31,30,31,30,31,31,30,31,30,31], dtype=np.float64)

    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    #-- find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    #-- calculate the day of the year
    if DofY is not None:
        #-- if entered directly as an input
        #-- remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        #-- use calendar month and day of the month to calculate day of the year
        #-- month minus 1: January = 0, February = 1, etc (indice of month)
        #-- in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=np.int64) - 1

        #-- day of month
        if day is not None:
            #-- remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            #-- if not entering days as an input
            #-- will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        #-- create matrix with the lower half = 1
        #-- this matrix will be used in a matrix multiplication
        #-- to calculate the total number of days for prior months
        #-- the -1 will make the diagonal == 0
        #-- i.e. first row == all zeros and the
        #-- last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        #-- using a dot product to calculate total number of days
        #-- for the months before the input date
        #-- basically is sum(i*dpm)
        #-- where i is 1 for all months < the month of interest
        #-- and i is 0 for all months >= the month of interest
        #-- month of interest is zero as the exact days will be
        #-- used to calculate the date

        #-- calculate the day of the year for leap and standard
        #-- use total days of all months before date
        #-- and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    #-- hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    #-- minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    #-- second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    #-- calculate decimal date
    #-- convert hours, minutes and seconds into days
    #-- convert calculated fractional days into decimal fractions of the year
    #-- Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    #-- Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

#-- PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD, ASTYPE=None, FORMAT='dict'):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: float
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    ASTYPE: str or NoneType, default None
        convert output to variable type
    FORMAT: str, default 'dict'
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

    #-- convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        SINGLE_VALUE = True
    else:
        SINGLE_VALUE = False

    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    #-- calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    #-- calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    #-- calculate day, month, year and hour
    DAY = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    MONTH = F - 1.0 - 12.0*np.floor(F/14.0)
    YEAR = D - 4715.0 - np.floor((7.0+MONTH)/10.0)
    HOUR = np.floor(24.0*(JD + 0.5 - JDO))
    #-- calculate minute and second
    G = (JD + 0.5 - JDO) - HOUR/24.0
    MINUTE = np.floor(G*1440.0)
    SECOND = (G - MINUTE/1440.0) * 86400.0

    #-- convert all variables to output type (from float)
    if ASTYPE is not None:
        YEAR = YEAR.astype(ASTYPE)
        MONTH = MONTH.astype(ASTYPE)
        DAY = DAY.astype(ASTYPE)
        HOUR = HOUR.astype(ASTYPE)
        MINUTE = MINUTE.astype(ASTYPE)
        SECOND = SECOND.astype(ASTYPE)

    #-- if only a single value was imported initially: remove singleton dims
    if SINGLE_VALUE:
        YEAR = YEAR.item(0)
        MONTH = MONTH.item(0)
        DAY = DAY.item(0)
        HOUR = HOUR.item(0)
        MINUTE = MINUTE.item(0)
        SECOND = SECOND.item(0)

    #-- return date variables in output format (default python dictionary)
    if (FORMAT == 'dict'):
        return dict(year=YEAR, month=MONTH, day=DAY,
            hour=HOUR, minute=MINUTE, second=SECOND)
    elif (FORMAT == 'tuple'):
        return (YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
    elif (FORMAT == 'zip'):
        return zip(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
