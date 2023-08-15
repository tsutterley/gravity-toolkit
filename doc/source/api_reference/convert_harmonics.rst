====================
convert_harmonics.py
====================

- Converts a file from the spatial domain into the spherical harmonic domain [Wahr1998]_

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/convert_harmonics.py

Calling Sequence
################

.. argparse::
    :filename: ../scripts/convert_harmonics.py
    :func: arguments
    :prog: convert_harmonics.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
        * ``1``: Gegout (2005) values from PREM [Gegout2010]_
        * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
        * ``3``: Wang et al. (2012) values from PREM with hard sediment [Wang2012]_
        * ``4``: Wang et al. (2012) values from PREM with soft sediment [Wang2012]_

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth

    --units -U : @after
        * ``1``: cm of water thickness (cm.w.e., g/cm\ :sup:`2`)
        * ``2``: Gigatonnes (Gt)
        * ``3``: mm of water thickness (kg/m\ :sup:`2`)

    --interval -I : @after
        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
