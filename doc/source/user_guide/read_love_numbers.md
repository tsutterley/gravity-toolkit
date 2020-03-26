read_love_numbers.py
====================

 - Reads sets of load Love numbers computed using outputs from the Preliminary Reference Earth Model (PREM) as described by [Han and Wahr (1995)](https://doi.org/10.1111/j.1365-246X.1995.tb01819.x)
 - Applies isomorphic parameters for different reference frames following [Blewitt (2003)](https://doi.org/10.1029/2002JB002082)

#### Calling Sequence
```python
from gravity_toolkit.read_love_numbers import read_love_numbers
hl,kl,ll = read_love_numbers(love_numbers_file, FORMAT='tuple', REFERENCE='CF')
```

#### Inputs
 - `love_numbers_file`: Elastic load Love numbers file

#### Options
 - `HEADER`: file contains header text to be skipped (default: True)
 - `REFERENCE`: Reference frame for calculating degree 1 love numbers
     * `'CF'`: Center of Surface Figure
     * `'CL'`: Center of Surface Lateral Figure
     * `'CH'`: Center of Surface Height Figure
     * `'CM'`: Center of Mass of Earth System
     * `'CE'`: Center of Mass of Solid Earth (default)
 - `FORMAT`: format of output variables
     * `'dict'`: dictionary with variable keys as listed above
     * `'tuple'`: tuple with variable order hl,kl,ll
     * `'zip'`: aggregated variable sets

#### Outputs
 - `kl`: Love number of Gravitational Potential
 - `hl`: Love number of Vertical Displacement
 - `ll`: Love number of Horizontal Displacement
