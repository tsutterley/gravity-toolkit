============
Spatial Maps
============

The GRACE/GRACE-FO spherical harmonic products can be converted into sets of spatial maps
if we assume that the mass redistributions are concentrated within a thin layer
(thickness |mlt| horizontal resolution) :cite:p:`Wahr:1998hy`.
To calculate accurate maps of spatial variability, several processing steps need to be accounted for
to convert the data into the proper reference frame, reduce the impact of noisy data,
remove unwanted sources of gravitational variability, and convert to appropriate units.

The ``grace_spatial_maps.py`` program will output spatial files in ascii, netCDF4 or HDF5 format
for each GRACE/GRACE-FO month.

Load Love Numbers
#################

A variation in mass at the Earth's surface will load and deform the solid Earth,
which will induce density anomalies at depth :cite:p:`Wahr:1998hy`.
To accurately assess the surface load from time-variable gravity,
we need to compensate for the Earth's elastic yielding :cite:p:`Wahr:1998hy`.
The elastic deformation of the solid Earth induced by a change in surface load
can be estimated using load Love numbers.
Using load Love numbers to calculate the elastic yielding assumes that
all other time-variable solid Earth contributions have been independently
removed from the spherical harmonic coefficients :cite:p:`Wahr:1998hy`.
Here, we use load Love and Shida numbers with parameters calculated from
the Preliminary Reference Earth model (PREM) :cite:p:`Farrell:1972cm,Dziewonski:1981bz`.
In order to help estimate the uncertainty in elastic deformation,
``grace_spatial_maps.py`` can use different sets of load Love numbers by adjusting the
``--love`` command line option.

Reference Frames
################

Measurements of time-variable gravity from the Gravity Recovery and Climate Experiment (GRACE)
and the GRACE Follow-On (GRACE-FO) missions are set in a center of mass (CM) reference frame,
in which the total degree one variations are inherently zero.
The individual contributions to degree one variations in the CM reference frame,
such as from oceanic processes or terrestrial water storage change, are not necessarily zero :cite:p:`Wahr:1998hy`.
Applications set in a center of figure (CF) reference frame,
such as the recovery of mass variations of the oceans, hydrosphere and cryosphere,
require the inclusion of degree one terms to be fully accurate :cite:p:`Swenson:2008cr`.

``grace_spatial_maps.py`` has geocenter options to select the degree one product to
include with the GRACE/GRACE-FO derived harmonics.
There are options for using measurements from satellite laser ranging :cite:p:`Cheng:2013tz` and
calculations from time-variable gravity and ocean model outputs :cite:p:`Swenson:2008cr,Sutterley:2019bx`.
If including degree one harmonics and changing the reference frame,
the reference frame for the load Love numbers needs to be updated accordingly :cite:p:`Blewitt:2003bz`.
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
from GRACE-A to GRACE-B and from GF1 to GF2 :cite:p:`Bandikova:2019ji`.
These single-accelerometer months for both GRACE and GRACE-FO contain significantly
more noise, particularly the low-degree zonal harmonics
(predominantly :math:`C_{20}` and :math:`C_{30}` but possibly :math:`C_{40}` and :math:`C_{50}`).
:math:`C_{20}` has also been difficult for GRACE and GRACE-FO to independently measure
throughout both missions.
The figure axis harmonics (:math:`C_{21}` and :math:`S_{21}`) may also be contaminated
by noise during the single-accelerometer months in the GFZ products :cite:p:`Dahle:2019jf`.
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
SLR low degree harmonic products :cite:p:`Cheng:2011hh,Cheng:2018jz,Koenig:2019gs,Loomis:2019dc,Loomis:2020bq`.

Corrections
###########

Prior to GRACE/GRACE-FO Release-6, corrections needed to be applied to compensate
for long-period signals in the pole tide that were contaminating the
:math:`C_{21}` and :math:`S_{21}` harmonics :cite:p:`Wahr:2015dg`,
as well as for discontinuities in the atmospheric de-aliasing product that were
introduced with upgrades in the ECMWF weather prediction model :cite:p:`Fagiolini:2015kc`.
The Pole Tide and Atmospheric corrections do not need to be applied to the Release-6 data.

Geophysical Leakage
###################

Gravity measurements from GRACE and GRACE-FO are global, near-monthly and
are directly related to changes in mass.
Several mass transport processes can occur concurrently for a given region,
which means that the total time-dependent geopotential from GRACE/GRACE-FO
can relate to multiple time-varying components :cite:p:`Wahr:1998hy`.
These mass transport processes include but are not limited to terrestrial water storage,
glacier and ice sheet mass, atmospheric and oceanic circulation and geodynamic processes.
In order to isolate the mass change of a single process, each of the other processes
needs to be independently estimated and removed from the GRACE/GRACE-FO data.
Uncertainties in the components removed from the GRACE/GRACE-FO data will directly
impact the precision of the final mass balance estimate.

Filtering
#########

The GRACE/GRACE-FO coefficients are impacted by random spherical harmonic errors
that increase as a function of spherical harmonic degree :cite:p:`Wahr:1998hy,Swenson:2002hs`.
The truncation of the spherical harmonics series also results
in spurious ringing artifacts from Gibbs phenomenon.
The impact of these errors can be reduced using Gaussian averaging functions
as described in :cite:p:`Jekeli:1981vj,Swenson:2002hs`.
GRACE/GRACE-FO coefficients are also impacted by correlated north/south "striping" errors,
which can be spectrally filtered following :cite:t:`Swenson:2006hu`.

.. |beta|    unicode:: U+03B2 .. GREEK SMALL LETTER BETA

.. |mu|      unicode:: U+03BC .. GREEK SMALL LETTER MU

.. |mlt|     unicode:: U+226A .. MUCH LESS-THAN