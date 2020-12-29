#!/usr/bin/env python
u"""
test_time.py (09/2020)
Verify time conversion functions
"""
import os
import pytest
import warnings
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.utilities

#-- parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
#-- PURPOSE: verify forward and backwards time conversions
def test_julian(YEAR,MONTH):
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    assert (np.sum(DPM) == gravity_toolkit.time.calendar_days(YEAR).sum())
    #-- calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    MJD = gravity_toolkit.time.convert_calendar_dates(YEAR, MONTH, DAY,
        hour=HOUR, minute=MINUTE, second=SECOND,
        epoch=(1858,11,17,0,0,0))
    #-- convert MJD to calendar date
    JD = np.squeeze(MJD) + 2400000.5
    YY,MM,DD,HH,MN,SS = gravity_toolkit.time.convert_julian(JD,
        FORMAT='tuple', ASTYPE=np.float)
    #-- assert dates
    eps = np.finfo(np.float16).eps
    assert (YY == YEAR)
    assert (MM == MONTH)
    assert (DD == DAY)
    assert (HH == HOUR)
    assert (MN == MINUTE)
    assert (np.abs(SS - SECOND) < eps)

#-- parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
#-- PURPOSE: verify forward and backwards time conversions
def test_decimal_dates(YEAR,MONTH):
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    assert (np.sum(DPM) == gravity_toolkit.time.calendar_days(YEAR).sum())
    #-- calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    #-- calculate year-decimal time
    tdec = gravity_toolkit.time.convert_calendar_decimal(YEAR, MONTH, day=DAY,
        hour=HOUR, minute=MINUTE, second=SECOND)
    #-- day of the year 1 = Jan 1, 365 = Dec 31 (std)
    day_temp = np.mod(tdec, 1)*np.sum(DPM)
    DofY = np.floor(day_temp) + 1
    #-- cumulative sum of the calendar dates
    day_cumulative = np.cumsum(np.concatenate(([0],DPM))) + 1
    #-- finding which month date is in
    i = np.nonzero((DofY >= day_cumulative[0:-1]) & (DofY < day_cumulative[1:]))
    month_range = np.arange(1,13)
    month = month_range[i]
    #-- finding day of the month
    day = (DofY - day_cumulative[i]) + 1
    #-- convert residuals into time (hour, minute and second)
    hour_temp = np.mod(day_temp,1)*24.0
    minute_temp = np.mod(hour_temp,1)*60.0
    second = np.mod(minute_temp,1)*60.0
    #-- assert dates
    eps = np.finfo(np.float16).eps
    assert (np.floor(tdec) == YEAR)
    assert (month == MONTH)
    assert (day == DAY)
    assert (np.floor(hour_temp) == HOUR)
    assert (np.floor(minute_temp) == MINUTE)
    assert (np.abs(second - SECOND) < eps)

#-- PURPOSE: test UNIX time
def test_unix_time():
    #-- ATLAS Standard Data Epoch
    UNIX = gravity_toolkit.utilities.get_unix_time('2018-01-01 00:00:00')
    assert (UNIX == 1514764800)

#-- PURPOSE: test parsing time strings
def test_parse_date_string():
    #-- time string for Modified Julian Days
    time_string = 'days since 1858-11-17T00:00:00'
    epoch,to_secs = gravity_toolkit.time.parse_date_string(time_string)
    #-- check the epoch and the time unit conversion factors
    assert np.all(epoch == [1858,11,17,0,0,0])
    assert (to_secs == 86400.0)
    #-- time string for ATLAS Standard Data Epoch
    time_string = 'seconds since 2018-01-01T00:00:00'
    epoch,to_secs = gravity_toolkit.time.parse_date_string(time_string)
    #-- check the epoch and the time unit conversion factors
    assert np.all(epoch == [2018,1,1,0,0,0])
    assert (to_secs == 1.0)

