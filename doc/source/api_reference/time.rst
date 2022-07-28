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

.. autofunction:: gravity_toolkit.time.parse_date_string

.. autofunction:: gravity_toolkit.time.split_date_string

.. autofunction:: gravity_toolkit.time.datetime_to_list

.. autofunction:: gravity_toolkit.time.adjust_months

.. autofunction:: gravity_toolkit.time.calendar_to_grace

.. autofunction:: gravity_toolkit.time.grace_to_calendar

.. autofunction:: gravity_toolkit.time.calendar_to_julian

.. autofunction:: gravity_toolkit.time.calendar_days

.. autofunction:: gravity_toolkit.time.convert_delta_time

.. autofunction:: gravity_toolkit.time.convert_calendar_dates

.. autofunction:: gravity_toolkit.time.convert_calendar_decimal

.. autofunction:: gravity_toolkit.time.convert_julian

