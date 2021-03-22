read_love_numbers.py
====================

- Reads sets of load Love numbers computed using outputs from the Preliminary Reference Earth Model (PREM) ([Dziewonski and Anderson, 1981](https://doi.org/10.1016/0031-9201(81)90046-7))
- Applies isomorphic parameters for different reference frames following [Blewitt (2003)](https://doi.org/10.1029/2002JB002082)
- Can read load Love numbers from:
    * [Han and Wahr (1995)](https://doi.org/10.1111/j.1365-246X.1995.tb01819.x)
    * [Gegout (2005)](http://gemini.gsfc.nasa.gov/aplo/)
    * [Wang et al. (2012)](https://doi.org/10.1016/j.cageo.2012.06.022)

#### Calling Sequence
```python
from gravity_toolkit.utilities import get_data_path
from gravity_toolkit.read_love_numbers import read_love_numbers
love_numbers_file = get_data_path(['data','love_numbers'])
hl,kl,ll = read_love_numbers(love_numbers_file, FORMAT='tuple', REFERENCE='CF')
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_love_numbers.py)

#### Arguments
- `love_numbers_file`: Elastic load Love numbers file

#### Keyword arguments
- `LMAX`: truncate or interpolate to maximum spherical harmonic degree
- `HEADER`: number of header lines to be skipped
- `COLUMNS`: column names of ascii file
    * `'l'`: spherical harmonic degree
    * `'hl'`: vertical displacement
    * `'kl'`: gravitational potential
    * `'ll'`: horizontal displacement
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

#### Returns
- `hl`: Love number of Vertical Displacement
- `kl`: Love number of Gravitational Potential
- `ll`: Love number of Horizontal Displacement
