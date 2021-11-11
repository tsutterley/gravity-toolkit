=======
time.py
=======

Utilities for calculating time operations

 - Can convert delta time from seconds since an epoch to time since a different epoch
 - Can calculate the time in days since epoch from calendar dates

Calling Sequence
================

Convert a time from seconds since 1980-01-06T00:00:00 to Modified Julian Days (MJD)

.. code-block:: python

    import gravity_toolkit.time
    MJD = gravity_toolkit.time.convert_delta_time(delta_time, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)

Convert from Julian Days into calendar dates

.. code-block:: python

    import gravity_toolkit.time
    JD = MJD + 2400000.5
    YEAR,MONTH,DAY,HOUR,MINUTE,SECOND = gravity_toolkit.time.convert_julian(JD,
        FORMAT='tuple')

Convert a calendar date into Modified Julian Days (MJD)

.. code-block:: python

    import gravity_toolkit.time
    MJD = gravity_toolkit.time.convert_calendar_dates(YEAR,MONTH,DAY,hour=HOUR,
        minute=MINUTE,second=SECOND,epoch=(1858,11,17,0,0,0))

Convert a calendar date into decimal years

.. code-block:: python

    import gravity_toolkit.time
    t_date = gravity_toolkit.time.convert_calendar_decimal(YEAR,MOMTH,day=DAY,
        hour=HOUR,minute=MINUTE,second=SECOND)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/time.py


General Methods
===============


.. method:: gravity_toolkit.time.parse_date_string(date_string)

    Parse a date string of the form time-units since yyyy-mm-dd hh:mm:ss

    Arguments:

        ``date_string``: time-units since yyyy-mm-dd hh:mm:ss

    Returns:

        ``epoch``: epoch of delta time

        ``scale``: multiplication factor to convert to seconds


.. method:: gravity_toolkit.time.split_date_string(date_string)

    Split a date string into units and epoch

    Arguments:

        ``date_string``: time-units since yyyy-mm-dd hh:mm:ss


.. method:: gravity_toolkit.time.datetime_to_list(date)

    Convert a datetime object into a list

    Arguments:

        ``date``: datetime object


.. method:: gravity_toolkit.time.adjust_months(grace_months)

    Adjust estimated GRACE/GRACE-FO months to fix "Special Cases"

    Arguments:

        ``grace_months``: GRACE/GRACE-FO months


.. method:: gravity_toolkit.time.calendar_to_grace(year,month=1,around=np.floor)

    Converts calendar dates to GRACE/GRACE-FO months

    Arguments:

        ``year``: calendar year

    Keyword arguments:

        ``month``: calendar month

        ``around``: method of rounding to nearest method

    Returns:

        ``grace_month``: GRACE/GRACE-FO month


.. method:: gravity_toolkit.time.grace_to_calendar(grace_month)

    Converts GRACE/GRACE-FO months to calendar dates

    Arguments:

        ``grace_month``: GRACE/GRACE-FO month

    Returns:

        ``year``: calendar year

        ``month``: calendar month


.. method:: gravity_toolkit.time.calendar_to_julian(year_decimal)

    Converts calendar dates to Julian days

    Arguments

        ``year_decimal``: calendar year

    Returns

        ``JD``: Julian Day (days since 01-01-4713 BCE at 12:00:00)


.. method:: gravity_toolkit.time.calendar_days(year)

    Calculates the number of days per month for a given year

    Arguments:

        ``year``: calendar year

    Returns:

        ``dpm``: number of days for each month


.. method:: gravity_toolkit.time.convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0)

    Convert delta time from seconds since epoch1 to time since epoch2

    Arguments:

        ``delta_time``: seconds since epoch1

    Keyword arguments:

        ``epoch1``: epoch for input delta_time

        ``epoch2``: epoch for output delta_time

        ``scale``: scaling factor for converting time to output units


.. method:: gravity_toolkit.time.convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0, epoch=None, scale=1.0)

    Calculate the time in time units since epoch from calendar dates

    Arguments:

        ``year``: calendar month

        ``month``: month of the year

        ``day``: day of the month

    Keyword arguments:

        ``hour``: hour of the day

        ``minute``: minute of the hour

        ``second``: second of the minute

        ``epoch``: epoch for output delta_time

        ``scale``: scaling factor for converting time to output units


.. method:: gravity_toolkit.time.convert_calendar_decimal(year, month, day=None, hour=None, minute=None, second=None, DofY=None)

    Converts from calendar date into decimal years taking into account leap years

    Arguments:

        ``year``: calendar year

        ``month``: calendar month

    Keyword arguments:

        ``day``: Number of day of the month

        ``hour``: hour of the day

        ``minute``: minute of the hour

        ``second``: second (and fractions of a second) of the minute

        ``DofY``: day of the year

    Returns:

        ``t_date`` date in decimal-year format


.. method:: gravity_toolkit.time.convert_julian(JD, ASTYPE=None, FORMAT=None)

    Converts from Julian day to calendar date and time

    Arguments:

        ``JD``: Julian Day (days since 01-01-4713 BCE at 12:00:00)

    Keyword arguments:

        ``ASTYPE``: convert output to variable type

        ``FORMAT``: format of output variables

            ``'dict'``: dictionary with variable keys

            ``'tuple'``: tuple with variable order year,month,day,hour,minute,second

            ``'zip'``: aggregated variable sets

    Returns:

        ``year``: Calendar year

        ``month``: Calendar month

        ``day``: Calendar day of the month

        ``hour``: hour of the day

        ``minute``: minute of the hour

        ``second``: second (and fractions of a second) of the minute
