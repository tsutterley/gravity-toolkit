read_love_numbers.py
====================

 - Reads sets of load Love numbers computed using outputs from the Preliminary Reference Earth Model (PREM) as described by [Han and Wahr (1995)](https://10.1111/j.1365-246X.1995.tb01819.x)

#### Calling Sequence
```
from gravity_toolkit.read_love_numbers import read_love_numbers
hl,kl,ll = read_love_numbers(love_numbers_file, FORMAT='tuple')
```

#### Inputs
 - `love_numbers_file`: Elastic load Love numbers file

#### Options
 - `HEADER`: file contains header text to be skipped (default: True)
 - `FORMAT`: format of output variables
     - `'dict'`: dictionary with variable keys as listed above
     - `'tuple'`: tuple with variable order hl,kl,ll
     - `'zip'`: aggregated variable sets

#### Outputs
 - `kl`: Love number of Gravitational Potential
 - `hl`: Love number of Vertical Displacement
 - `ll`: Love number of Horizontal Displacement
