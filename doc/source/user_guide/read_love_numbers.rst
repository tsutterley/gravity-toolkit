====================
read_love_numbers.py
====================

- Reads sets of load Love numbers computed using outputs from the Preliminary Reference Earth Model (PREM) [Dziewonski1981]_
- Applies isomorphic parameters for different reference frames following [Blewett2003]_
- Can read load Love numbers from [Han1995]_, [Gegout2010]_, and [Wang2012]_

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.utilities import get_data_path
    from gravity_toolkit.read_love_numbers import read_love_numbers
    love_numbers_file = get_data_path(['data','love_numbers'])
    hl,kl,ll = read_love_numbers(love_numbers_file, FORMAT='tuple', REFERENCE='CF')

`Source code`__

    .. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/read_love_numbers.py

Arguments
#########

- ``love_numbers_file``: Elastic load Love numbers file

Keyword arguments
#################

- ``LMAX``: truncate or interpolate to maximum spherical harmonic degree
- ``HEADER``: number of header lines to be skipped
- ``COLUMNS``: column names of ascii file

    * ``'l'``: spherical harmonic degree
    * ``'hl'``: vertical displacement
    * ``'kl'``: gravitational potential
    * ``'ll'``: horizontal displacement
- ``REFERENCE``: Reference frame for calculating degree 1 love numbers

    * ``'CF'``: Center of Surface Figure
    * ``'CL'``: Center of Surface Lateral Figure
    * ``'CH'``: Center of Surface Height Figure
    * ``'CM'``: Center of Mass of Earth System
    * ``'CE'``: Center of Mass of Solid Earth (default)
- ``FORMAT``: format of output variables

    * ``'dict'``: dictionary with variable keys as listed above
    * ``'tuple'``: tuple with variable order ``hl``, ``kl``, ``ll``
    * ``'zip'``: aggregated variable sets

Returns
#######

- ``hl``: Love number of Vertical Displacement
- ``kl``: Love number of Gravitational Potential
- ``ll``: Love number of Horizontal Displacement

References
##########

.. [Blewett2003] G. Blewitt, "Self‐consistency in reference frames, geocenter definition, and surface loading of the solid Earth", *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003). `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

.. [Dziewonski1981] A. M. Dziewonski and D. L. Anderson, "Preliminary reference Earth model", *Physics of the Earth and Planetary Interiors*, 25(4), 297--356, (1981). `doi: 10.1016/0031-9201(81)90046-7 <https://doi.org/10.1016/0031-9201(81)90046-7>`_

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
