============
Spatial Maps
============

The GRACE/GRACE-FO spherical harmonic products can be converted into sets of spatial maps
if we assume that the mass redistributions are concentrated within a thin layer
(thickness |mlt| horizontal resolution) [Wahr1998]_.
To calculate accurate maps of spatial variability, several processing steps need to be accounted for
to convert the data into the proper reference frame, reduce the impact of noisy data,
remove unwanted sources of gravitational variability, and convert to appropriate units.

The ``grace_spatial_maps.py`` program will output spatial files in ascii, netCDF4 or HDF5 format
for each GRACE/GRACE-FO month.

Load Love Numbers
#################

A variation in mass at the Earth's surface will load and deform the solid Earth,
which will induce density anomalies at depth [Wahr1998]_.
To accurately assess the surface load from time-variable gravity,
we need to compensate for the Earth's elastic yielding [Wahr1998]_.
The elastic deformation of the solid Earth induced by a change in surface load
can be estimated using load Love numbers.
Using load Love numbers to calculate the elastic yielding assumes that
all other time-variable solid Earth contributions have been independently
removed from the spherical harmonic coefficients [Wahr1998]_.
Here, we use load Love numbers with parameters calculated from
the Preliminary Reference Earth model (PREM) [Farrell1972]_ [Dziewonski1981]_.
In order to help estimate the uncertainty in elastic deformation,

``grace_spatial_maps.py`` can use different sets of load Love numbers by adjusting the
``--love`` command line option.

Reference Frames
################

Measurements of time-variable gravity from the Gravity Recovery and Climate Experiment (GRACE)
and the GRACE Follow-On (GRACE-FO) missions are set in a center of mass (CM) reference frame,
in which the total degree one variations are inherently zero.
The individual contributions to degree one variations in the CM reference frame,
such as from oceanic processes or terrestrial water storage change, are not necessarily zero [Wahr1998]_.
Applications set in a center of figure (CF) reference frame,
such as the recovery of mass variations of the oceans, hydrosphere and cryosphere,
require the inclusion of degree one terms to be fully accurate [Swenson2008]_.

``grace_spatial_maps.py`` has geocenter options to select the degree one product to
include with the GRACE/GRACE-FO derived harmonics.
There are options for using measurements from satellite laser ranging [Cheng2013]_ and
calculations from time-variable gravity and ocean model outputs [Swenson2008]_ [Sutterley2019]_.
If including degree one harmonics and changing the reference frame,
the reference frame for the load Love numbers needs to be updated accordingly [Blewett2003]_.
In ``grace_spatial_maps.py`` and other GRACE/GRACE-FO programs, the reference frame for the load Love numbers
is adjusted by setting the ``--reference`` command line option to ``'CF'``.

Low-Degree Harmonics
####################

For both GRACE and GRACE-FO, there have been operational issues that have affected the
quality of the time-variable gravity fields.
During the late stages of the GRACE mission, procedures were enacted to preserve the
battery life of the GRACE satellites and extend the mission lifetime.
This included turning off the accelerometer and microwave ranging instrument (MWI) to
reduce the battery load on the spacecraft during periods of low |beta| angle when solar
input to the spacecraft is lowest.
In late 2016, the accelerometer onboard GRACE-B was permanently powered down to help
maintain the operation of the ranging instrument.
For the Follow-On mission, an anomaly shortly after the launch of GRACE-FO caused a
malfunction and the accelerometer onboard GRACE-FO 2 (GF2) has been operating in a
less optimal Large-Range-Mode.
For both of these cases, the GRACE/GRACE-FO processing centers have developed
independent methods to spatiotemporally transplant the accelerometer data retrieved
from GRACE-A to GRACE-B and from GF1 to GF2 [Bandikova2019]_.
These single-accelerometer months for both GRACE and GRACE-FO contain significantly
more noise, particularly the low-degree zonal harmonics
(predominantly :math:`C_{20}` and :math:`C_{30}` but possibly :math:`C_{40}` and :math:`C_{50}`).
:math:`C_{20}` has also been difficult for GRACE and GRACE-FO to independently measure
throughout both missions.
The figure axis harmonics (:math:`C_{21}` and :math:`S_{21}`) may also be contaminated
by noise during the single-accelerometer months in the GFZ products [Dahle2019]_.
Measurements from satellite laser ranging (SLR) can provide an independent assessment
for some low degree and order spherical harmonics.
``grace_spatial_maps.py`` has options for replacing
:math:`C_{20}`,
:math:`C_{21}`,
:math:`S_{21}`,
:math:`C_{22}`,
:math:`S_{22}`,
:math:`C_{30}`,
:math:`C_{40}`,
and :math:`C_{50}` with
SLR low degree harmonic products [Cheng2011]_ [Cheng2018]_ [Koenig2019]_ [Loomis2019]_ [Loomis2020]_.

Corrections
###########

Prior to GRACE/GRACE-FO Release-6, corrections needed to be applied to compensate
for long-period signals in the pole tide that were contaminating the
:math:`C_{21}` and :math:`S_{21}` harmonics [Wahr2015]_,
as well as for discontinuities in the atmospheric de-aliasing product that were
introduced with upgrades in the ECMWF weather prediction model [Fagiolini2015]_.
The Pole Tide and Atmospheric corrections do not need to be applied to the Release-6 data.

Geophysical Leakage
###################

Gravity measurements from GRACE and GRACE-FO are global, near-monthly and
are directly related to changes in mass.
Several mass transport processes can occur concurrently for a given region,
which means that the total time-dependent geopotential from GRACE/GRACE-FO
can relate to multiple time-varying components [Wahr1998]_.
These mass transport processes include but are not limited to terrestrial water storage,
glacier and ice sheet mass, atmospheric and oceanic circulation and geodynamic processes.
In order to isolate the mass change of a single process, each of the other processes
needs to be independently estimated and removed from the GRACE/GRACE-FO data.
Uncertainties in the components removed from the GRACE/GRACE-FO data will directly
impact the precision of the final mass balance estimate.

Filtering
#########

The GRACE/GRACE-FO coefficients are impacted by random spherical harmonic errors
that increase as a function of spherical harmonic degree [Wahr1998]_ [Swenson2002]_.
The truncation of the spherical harmonics series also results
in spurious ringing artifacts from Gibbs phenomenon.
The impact of these errors can be reduced using Gaussian averaging functions
as described in [Jekeli1981]_ [Swenson2002]_.
GRACE/GRACE-FO coefficients are also impacted by correlated north/south "striping" errors,
which can be spectrally filtered following [Swenson2006]_.

References
##########

.. [Bandikova2019] T. Bandikova, C. McCullough, G. L. Kruizinga, H. Save, and B. Christophe, "GRACE accelerometer data transplant", *Advances in Space Research*, 64(3), 623--644, (2019). `doi: 10.1016/j.asr.2019.05.021 <10.1016/j.asr.2019.05.021>`_

.. [Blewett2003] G. Blewitt, "Self‐consistency in reference frames, geocenter definition, and surface loading of the solid Earth", *Journal of Geophysical Research: Solid Earth*, 108(B2), 2103, (2003). `doi: 10.1029/2002JB002082 <https://doi.org/10.1029/2002JB002082>`_

.. [Cheng2011] M. Cheng, J. C. Ries, and B. D. Tapley, "Variations of the Earth's figure axis from satellite laser ranging and GRACE", *Journal of Geophysical Research: Solid Earth*, 116, B01409, (2011). `doi: 10.1029/2010JB000850 <https://doi.org/10.1029/2010JB000850>`_

.. [Cheng2013] M. Cheng, "Geocenter Variations from Analysis of SLR Data", *Reference Frames for Applications in Geosciences*, 19--25, (2013). `doi: 10.1007/978-3-642-32998-2_4 <https://doi.org/10.1007/978-3-642-32998-2_4>`_

.. [Cheng2018] M. Cheng and J. C. Ries, "Decadal variation in Earth's oblateness (J2) from satellite laser ranging data", *Geophysical Journal International*, 212(2), 1218--1224 (2018). `doi: 10.1093/gji/ggx483 <https://doi.org/10.1093/gji/ggx483>`_

.. [Dahle2019] C. Dahle et al. "The GFZ GRACE RL06 Monthly Gravity Field Time Series: Processing Details, and Quality Assessment", *Remote Sensing*, 11(18), 2116, (2019). `doi: 10.3390/rs11182116 <https://doi.org/10.3390/rs11182116>`_

.. [Dziewonski1981] A. M. Dziewonski and D. L. Anderson, "Preliminary reference Earth model", *Physics of the Earth and Planetary Interiors*, 25(4), 297--356, (1981). `doi: 10.1016/0031-9201(81)90046-7 <https://doi.org/10.1016/0031-9201(81)90046-7>`_

.. [Fagiolini2015] E. Fagiolini, F. Flechtner, M. Horwath, and H. Dobslaw, "Correction of inconsistencies in ECMWF's operational analysis data during de-aliasing of GRACE gravity models", *Geophysical Journal International*, 202(3), 2150--2158, (2015). `doi: 10.1093/gji/ggv276 <https://doi.org/10.1093/gji/ggv276>`_

.. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads", *Reviews of Geophysics*, 10(3), 761--797, (1972). `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Koenig2019] R. Koenig, P. Schreiner, and C. Dahle, "Monthly estimates of C(2,0) generated by GFZ from SLR satellites based on GFZ GRACE/GRACE-FO RL06 background models", V. 1.0. GFZ Data Services, (2019). `doi: 10.5880/GFZ.GRAVIS_06_C20_SLR <http://doi.org/10.5880/GFZ.GRAVIS_06_C20_SLR>`_

.. [Loomis2019] B. D. Loomis, K. E. Rachlin, and S. B. Luthcke, "Improved Earth oblateness rate reveals increased ice sheet losses and mass‐driven sea level rise". *Geophysical Research Letters*, 46, 6910--6917, (2019). `doi: 10.1029/2019GL082929 <https://doi.org/10.1029/2019GL082929>`_

.. [Loomis2020] B. D. Loomis, K. E. Rachlin, D. N. Wiese, F. W. Landerer, and S. B. Luthcke, "Replacing GRACE/GRACE‐FO *C*\ :sub:`30` with satellite laser ranging: Impacts on Antarctic Ice Sheet mass change". *Geophysical Research Letters*, 47, (2020). `doi: 10.1029/2019GL085488 <https://doi.org/10.1029/2019GL085488>`_

.. [Sutterley2019] T. C. Sutterley and I. Velicogna, "Improved Estimates of Geocenter Variability from Time-Variable Gravity and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

.. [Swenson2002] S. Swenson and J. Wahr, "Methods for inferring regional surface‐mass anomalies from Gravity Recovery and Climate Experiment (GRACE) measurements of time‐variable gravity", *Journal of Geophysical Research: Solid Earth*, 107(B9), 2193, (2002). `doi: 10.1029/2001JB000576 <https://doi.org/10.1029/2001JB000576>`_

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wahr2015] J. Wahr, R. S. Nerem, and S. V. Bettadpur, "The pole tide and its effect on GRACE time‐variable gravity measurements: Implications for estimates of surface mass variations". *Journal of Geophysical Research: Solid Earth*, 120(6), 4597--4615, (2015). `doi: 10.1002/2015JB011986 <https://doi.org/10.1002/2015JB011986>`_

.. |beta|    unicode:: U+03B2 .. GREEK SMALL LETTER BETA
.. |mu|      unicode:: U+03BC .. GREEK SMALL LETTER MU
.. |mlt|     unicode:: U+226A .. MUCH LESS-THAN