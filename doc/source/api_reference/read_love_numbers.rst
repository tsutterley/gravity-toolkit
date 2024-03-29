=================
read_love_numbers
=================

- Reads sets of load Love/Shida numbers computed using outputs from the Preliminary Reference Earth Model (PREM) or other Earth models
- Linearly interpolates load Love/Shida numbers for missing degrees
- Linearly extrapolates load Love/Shida numbers beyond maximum degree of dataset
- Applies isomorphic parameters for different reference frames

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.utilities import get_data_path
    from gravity_toolkit.read_love_numbers import read_love_numbers
    love_numbers_file = get_data_path(['data','love_numbers'])
    hl,kl,ll = read_love_numbers(love_numbers_file, FORMAT='tuple', REFERENCE='CF')

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/read_love_numbers.py

.. autofunction:: gravity_toolkit.read_love_numbers

.. autofunction:: gravity_toolkit.read_love_numbers.love_number_formatter

.. autofunction:: gravity_toolkit.read_love_numbers.extract_love_numbers

.. autofunction:: gravity_toolkit.read_love_numbers.load_love_numbers

.. autoclass:: gravity_toolkit.love_numbers
   :members:
