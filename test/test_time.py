#!/usr/bin/env python
u"""
test_time.py (01/2023)
Verify time conversion and utility functions

UPDATE HISTORY:
    Updated 01/2023: single implicit import of gravity toolkit
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: added function to test GRACE months adjustments
        test date parser for cases when only a date and no units
    Written 12/2020
"""
import pytest
import warnings
import numpy as np
import gravity_toolkit as gravtk

# parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
# PURPOSE: verify forward and backwards time conversions
def test_julian(YEAR,MONTH):
    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    assert (np.sum(DPM) == gravtk.time.calendar_days(YEAR).sum())
    # calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    MJD = gravtk.time.convert_calendar_dates(YEAR, MONTH, DAY,
        hour=HOUR, minute=MINUTE, second=SECOND,
        epoch=(1858,11,17,0,0,0))
    # convert MJD to calendar date
    JD = np.squeeze(MJD) + 2400000.5
    YY,MM,DD,HH,MN,SS = gravtk.time.convert_julian(JD,
        format='tuple', astype=np.float64)
    # assert dates
    eps = np.finfo(np.float16).eps
    assert (YY == YEAR)
    assert (MM == MONTH)
    assert (DD == DAY)
    assert (HH == HOUR)
    assert (MN == MINUTE)
    assert (np.abs(SS - SECOND) < eps)

# parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
# PURPOSE: verify forward and backwards time conversions
def test_decimal_dates(YEAR,MONTH):
    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    assert (np.sum(DPM) == gravtk.time.calendar_days(YEAR).sum())
    # calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    # calculate year-decimal time
    tdec = gravtk.time.convert_calendar_decimal(YEAR, MONTH, day=DAY,
        hour=HOUR, minute=MINUTE, second=SECOND)
    # day of the year 1 = Jan 1, 365 = Dec 31 (std)
    day_temp = np.mod(tdec, 1)*np.sum(DPM)
    DofY = np.floor(day_temp) + 1
    # cumulative sum of the calendar dates
    day_cumulative = np.cumsum(np.concatenate(([0],DPM))) + 1
    # finding which month date is in
    i = np.nonzero((DofY >= day_cumulative[0:-1]) & (DofY < day_cumulative[1:]))
    month_range = np.arange(1,13)
    month = month_range[i]
    # finding day of the month
    day = (DofY - day_cumulative[i]) + 1
    # convert residuals into time (hour, minute and second)
    hour_temp = np.mod(day_temp,1)*24.0
    minute_temp = np.mod(hour_temp,1)*60.0
    second = np.mod(minute_temp,1)*60.0
    # assert dates
    eps = np.finfo(np.float16).eps
    assert (np.floor(tdec) == YEAR)
    assert (month == MONTH)
    assert (day == DAY)
    assert (np.floor(hour_temp) == HOUR)
    assert (np.floor(minute_temp) == MINUTE)
    assert (np.abs(second - SECOND) < eps)

# PURPOSE: test UNIX time
def test_unix_time():
    # ATLAS Standard Data Epoch
    UNIX = gravtk.utilities.get_unix_time('2018-01-01 00:00:00')
    assert (UNIX == 1514764800)

# PURPOSE: test parsing time strings
def test_parse_date_string():
    # time string for Modified Julian Days
    time_string = 'days since 1858-11-17T00:00:00'
    epoch,to_secs = gravtk.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [1858,11,17,0,0,0])
    assert (to_secs == 86400.0)
    # time string for ATLAS Standard Data Epoch
    time_string = 'seconds since 2018-01-01T00:00:00'
    epoch,to_secs = gravtk.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [2018,1,1,0,0,0])
    assert (to_secs == 1.0)
    # time string for unitless case
    time_string = '2000-01-01T12:00:00'
    epoch,to_secs = gravtk.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [2000,1,1,12,0,0])
    assert (to_secs == 0.0)

# PURPOSE: test months adjustment for special cases
# parameterize calendar dates
@pytest.mark.parametrize("PROC", ['CSR','GFZ','GSFC','JPL'])
def test_adjust_months(PROC):
    # The 'Special Months' (Nov 2011, Dec 2011 and April 2012) with
    # Accelerometer shutoffs make the relation between month number
    # and date more complicated as days from other months are used
    # For CSR and GFZ: Nov 2011 (119) is centered in Oct 2011 (118)
    # For JPL: Dec 2011 (120) is centered in Jan 2012 (121)
    # For all: May 2015 (161) is centered in Apr 2015 (160)
    # For GSFC: Oct 2018 (202) is centered in Nov 2018 (203)

    # dates with special months for each processing center
    center_dates = dict(CSR=[],GFZ=[],GSFC=[],JPL=[])
    # CSR dates to test (year-decimal, GRACE month)
    center_dates['CSR'].append([2011.62465753, 116])
    center_dates['CSR'].append([2011.70821918, 117])
    center_dates['CSR'].append([2011.79178082, 118])
    center_dates['CSR'].append([2011.83287671, 119])
    center_dates['CSR'].append([2011.99041096, 120])
    center_dates['CSR'].append([2012.04371585, 121])
    center_dates['CSR'].append([2012.12568306, 122])
    center_dates['CSR'].append([2015.06027397, 157])
    center_dates['CSR'].append([2015.12465753, 158])
    center_dates['CSR'].append([2015.20547945, 159])
    center_dates['CSR'].append([2015.28904110, 160])
    center_dates['CSR'].append([2015.31917808, 161])
    center_dates['CSR'].append([2015.53698630, 163])
    center_dates['CSR'].append([2015.62465753, 164])
    # GFZ dates to test (year-decimal, GRACE month)
    center_dates['GFZ'].append([2011.62465753, 116])
    center_dates['GFZ'].append([2011.70821918, 117])
    center_dates['GFZ'].append([2011.79178082, 118])
    center_dates['GFZ'].append([2011.83287671, 119])
    center_dates['GFZ'].append([2011.99041096, 120])
    center_dates['GFZ'].append([2012.04371585, 121])
    center_dates['GFZ'].append([2012.12568306, 122])
    center_dates['GFZ'].append([2015.06027397, 157])
    center_dates['GFZ'].append([2015.12465753, 158])
    center_dates['GFZ'].append([2015.20547945, 159])
    center_dates['GFZ'].append([2015.28904110, 160])
    center_dates['GFZ'].append([2015.31917808, 161])
    center_dates['GFZ'].append([2015.53698630, 163])
    center_dates['GFZ'].append([2015.62465753, 164])
    # GSFC dates to test (year-decimal, GRACE month)
    center_dates['GSFC'].append([2018.45479452, 198])
    center_dates['GSFC'].append([2018.53835616, 199])
    center_dates['GSFC'].append([2018.84794520, 202])
    center_dates['GSFC'].append([2018.87397260, 203])
    center_dates['GSFC'].append([2018.95753424, 204])
    center_dates['GSFC'].append([2019.04109589, 205])
    # JPL dates to test (year-decimal, GRACE month)
    center_dates['JPL'].append([2011.62465753, 116])
    center_dates['JPL'].append([2011.70821918, 117])
    center_dates['JPL'].append([2011.79178082, 118])
    center_dates['JPL'].append([2011.83561644, 119])
    center_dates['JPL'].append([2012.00136986, 120])
    center_dates['JPL'].append([2012.04371585, 121])
    center_dates['JPL'].append([2012.12568306, 122])
    center_dates['JPL'].append([2015.06027397, 157])
    center_dates['JPL'].append([2015.12465753, 158])
    center_dates['JPL'].append([2015.20547945, 159])
    center_dates['JPL'].append([2015.28904110, 160])
    center_dates['JPL'].append([2015.31917808, 161])
    center_dates['JPL'].append([2015.53698630, 163])
    center_dates['JPL'].append([2015.62465753, 164])

    # get dates and months for center
    tdec,months = np.transpose(center_dates[PROC])
    # GRACE/GRACE-FO months with duplicates
    temp = np.array(12.0*(tdec-2002.0)+1,dtype='i')
    assert np.any(temp != months.astype('i'))
    # run months adjustment to fix special cases
    temp = gravtk.time.adjust_months(temp)
    assert np.all(temp == months.astype('i'))
