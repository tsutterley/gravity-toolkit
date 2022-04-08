=========================
run_sea_level_equation.py
=========================

- Solves the sea level equation with the option of including polar motion feedback [Farrell1976]_ [Kendall2005]_ [Mitrovica2003]_
- Uses a Clenshaw summation to calculate the spherical harmonic summation [Holmes2002]_ [Tscherning1982]_

Calling Sequence
################

.. code-block:: bash

     python run_sea_level_equation --format netCDF4 --lmax 512 input_file output_file

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/run_sea_level_equation.py

Command Line Options
####################

- ``--mask X``: input land-sea function (netCDF4) with variable ``LSMASK`` as mask
- ``-l X``, ``--lmax X``: Maximum spherical harmonic degree
- ``-n X``, ``--love X``: Treatment of the Love Love numbers

     * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
     * ``1``: Gegout (2005) values from PREM [Gegout2010]_
     * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``-b X``, ``--body X``: Treatment of the body tide Love number

        - ``0``: [Wahr1981]_ and [Wahr1985]_ values from PREM
        - ``1``: [Farrell1972]_ values from Gutenberg-Bullen oceanic mantle model
- ``-f X``, ``--fluid X:`` Treatment of the fluid Love number

        - ``0``: [Han1989]_ fluid love number
        - ``1``: [Munk1960]_ secular love number
        - ``2``: [Munk1960]_ fluid love number
        - ``3``: [Lambeck1980]_ fluid love number
- ``--polar-feedback:`` Include polar feedback [Wahr1985]_
- ``--reference X``: Reference frame for load love numbers [Blewett2003]_

     * ``'CF'``: Center of Surface Figure (default)
     * ``'CM'``: Center of Mass of Earth System
     * ``'CE'``: Center of Mass of Solid Earth
- ``-I X``, ``--iterations X``: maximum number of iterations for the solver
- ``-F X``, ``--format X``: input and output data format

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created
-

References
##########

.. [Blewett2003] G. Blewitt, "Self‐consistency in reference frames, geocenter definition, and surface loading of the solid Earth", *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003). `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

.. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads", *Reviews of Geophysics*, 10(3), 761--797, (1972). `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_

.. [Farrell1976] W. E. Farrell and J. A. Clark, "On Postglacial Sea Level", *Geophysical Journal of the Royal Astronomical Society*, 46(3), 647--667, (1976). `doi: 10.1111/j.1365-246X.1976.tb01252.x <https://doi.org/10.1111/j.1365-246X.1976.tb01252.x>`_

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1989] D. Han and J. Wahr, "Post-Glacial Rebound Analysis for a Rotating Earth", *Slow Deformation and Transmission of Stress in the Earth*, 49, (1989). `doi: 10.1029/GM049p0001 <https://doi.org/10.1029/GM049p0001>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions", *Journal of Geodesy*, 76, 279--299, (2002). `doi: 10.1007/s00190-002-0216-2 <https://doi.org/10.1007/s00190-002-0216-2>`_

.. [Kendall2005] R. A. Kendall, J. X. Mitrovica, and G. A. Milne, "On post-glacial sea level -- II. Numerical formulation and comparative results on spherically symmetric models", *Geophysical Journal International*, 161(3), 679--706, (2005). `doi: 10.1111/j.1365-246X.2005.02553.x <https://doi.org/10.1111/j.1365-246X.2005.02553.x>`_

.. [Lambeck1980] K. Lambeck, *The Earth's Variable Rotation: Geophysical Causes and Consequences*, First Edition, (1980).

.. [Mitrovica2003] J. X. Mitrovica and G. A. Milne, "On post-glacial sea level: I. General theory", *Geophysical Journal International*, 154(2), 253--267, (2003). `doi: 10.1046/j.1365-246X.2003.01942.x <https://doi.org/10.1046/j.1365-246X.2003.01942.x>`_

.. [Munk1960] W. H. Munk and G. J. F. MacDonald, *The Rotation of the Earth: A Geophysical Discussion*, First Edition, (1960).

.. [Tscherning1982] C. C. Tscherning and K. Poder, "Some Geodetic Applications of Clenshaw Summation", *Bollettino di Geodesia e Scienze*, 4, 349--375, (1982).

.. [Wahr1981] J. M. Wahr, "Body tides on an elliptical, rotating, elastic and oceanless Earth", *Geophysical Journal of the Royal Astronomical Society*, 64(3), 677--703, (1981). `doi: 10.1111/j.1365-246X.1981.tb02690.x <https://doi.org/10.1111/j.1365-246X.1981.tb02690.x>`_

.. [Wahr1985] J. M. Wahr, "Deformation induced by polar motion", *Journal of Geophysical Research: Solid Earth*, 90(B11), 9363--9368, (1985). `doi: 10.1029/JB090iB11p09363 <https://doi.org/10.1029/JB090iB11p09363>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
