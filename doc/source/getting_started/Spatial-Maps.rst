============
Spatial Maps
============

The GRACE/GRACE-FO spherical harmonic products can be converted into sets of spatial maps
if we assume that the mass redistributions are concentrated within a thin layer
(thickness |mlt| horizontal resolution) (`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
To calculate accurate maps of spatial variability, several processing steps need to be accounted for
to convert the data into the proper reference frame, reduce the impact of noisy data,
remove unwanted sources of gravitational variability, and convert to appropriate units.

The `grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
program will output individual files for each GRACE/GRACE-FO month.
The monthly spatial files can be in ascii, netCDF4 or HDF5 formats.

Load Love Numbers
#################

A variation in mass at the Earth's surface will load and deform the solid Earth,
which will induce density anomalies at depth (`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
To accurately assess the surface load from time-variable gravity,
we need to compensate for the Earth's elastic yielding (`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
This program accounts for the elastic deformation of the solid Earth using load Love numbers
with parameters from the Preliminary reference Earth model (PREM)
(`Farrell, 1972 <http://dx.doi.org/10.1029/RG010i003p00761>`_;
`Dziewonski and Anderson, 1981 <https://doi.org/10.1016/0031-9201(81)90046-7>`_).
In order to help estimate the uncertainty in elastic deformation,
`grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
can use different sets of load Love numbers by adjusting the ``--love`` command line option.

Reference Frames
################

Measurements of time-variable gravity from the Gravity Recovery and Climate Experiment (GRACE)
and the GRACE Follow-On (GRACE-FO) missions are set in a center of mass (CM) reference frame,
in which the total degree one variations are inherently zero.
The individual contributions to degree one variations in the CM reference frame,
such as from oceanic processes or terrestrial water storage change, are not necessarily zero
(`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
Applications set in a center of figure (CF) reference frame,
such as the recovery of mass variations of the oceans, hydrosphere and cryosphere,
require the inclusion of degree one terms to be fully accurate
(`Swenson et al., 2008 <https://doi.org/10.1029/2007JB005338>`_).
`grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
has geocenter options to select the degree one product to include with the GRACE/GRACE-FO derived harmonics.
There are options for using measurements from
satellite laser ranging (`Cheng, 2013 <https://doi.org/10.1007/978-3-642-32998-2_4>`_) and
calculations from time-variable gravity and ocean model outputs
(`Swenson et al., 2008 <https://doi.org/10.1029/2007JB005338>`_;
`Sutterley and Velicogna, 2019 <https://doi.org/10.3390/rs11182108>`_).
If including degree one harmonics and changing the reference frame,
the reference frame for the load Love numbers needs to be updated accordingly
(`Blewett, 2003 <https://doi.org/10.1029/2002JB002082>`_).
In `grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
and other GRACE/GRACE-FO programs, the reference frame for the load Love numbers
is adjusted by setting the ``--reference`` command line option to ``'CF'``.

Low-Degree Zonals
#################

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
from GRACE-A to GRACE-B and from GF1 to GF2 (`Bandikova et al., 2019 <https://doi.org/10.1016/j.asr.2019.05.021>`_).
These single-accelerometer months for both GRACE and GRACE-FO contain significantly
more noise, particularly the low-degree zonal harmonics
(particularly *C*\ :sub:`20` and *C*\ :sub:`30` but possibly *C*\ :sub:`40` and *C*\ :sub:`50`).
*C*\ :sub:`20` has also been difficult for GRACE and GRACE-FO to independently measure throughout both missions.
Measurements from satellite laser ranging (SLR) can provide an independent assessment
for some low degree and order spherical harmonics
(`Cheng and Ries, 2018 <https://doi.org/10.1093/gji/ggx483>`_;
`Loomis et al., 2019 <https://doi.org/10.1029/2019GL082929>`_;
`Loomis et al., 2020 <https://doi.org/10.1029/2019GL085488>`_).
`grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
has options for replacing both *C*\ :sub:`20` and *C*\ :sub:`30` with SLR low degree zonal products.

Corrections
###########

Prior to GRACE/GRACE-FO Release-6, corrections needed to be applied to compensate
for long-period signals in the pole tide that were contaminating the
*C*\ :sub:`21` and *S*\ :sub:`21` harmonics
(`Wahr et al., 2015 <https://doi.org/10.1002/2015JB011986>`_),
as well as for discontinuities in the atmospheric de-aliasing product that were
introduced with upgrades in the ECMWF weather prediction model
(`Fagiolini et al., 2015 <https://doi.org/10.1093/gji/ggv276>`_).
The Pole Tide and Atmospheric corrections do not need to be applied to the Release-6 data.

Geophysical Leakage
###################

Gravity measurements from GRACE and GRACE-FO are global, near-monthly and
are directly related to changes in mass.
Several mass transport processes can occur concurrently for a given region,
which means that the total time-dependent geopotential from GRACE/GRACE-FO
can relate to multiple time-varying components
(`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
These mass transport processes include but are not limited to terrestrial water storage,
glacier and ice sheet mass, atmospheric and oceanic circulation and geodynamic processes.
In order to isolate the mass change of a single process, each of the other processes
needs to be independently estimated and removed from the GRACE/GRACE-FO data.
Uncertainties in the components removed from the GRACE/GRACE-FO data will directly
impact the precision of the final mass balance estimate.

Filtering
#########

The GRACE/GRACE-FO coefficients are impacted by random spherical harmonic errors
that increase as a function of spherical harmonic degree
(`Wahr et al., 1998 <https://doi.org/10.1029/98JB02844>`_).
The impact of these errors can be reduced using Gaussian averaging functions as described in
`Jekeli, (1981) <http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf>`_.
GRACE/GRACE-FO coefficients are also impacted by correlated north/south "striping" errors,
which can be spectrally filtered following `Swenson and Wahr (2006) <https://doi.org/10.1029/2005GL025285>`_.

Parameter Files
###############

The `grace_spatial_maps.py <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/scripts/grace_spatial_maps.py>`_
program works by accepting parameter file system arguments (``sys.argv``) listed after the program.
These parameter files can be ran individually or in a series.

.. code-block:: bash

     python grace_spatial_maps.py inp1 inp2 inp3
     python grace_spatial_maps.py --np 2 inp1 inp2 inp3
     python grace_spatial_maps.py -P 2 inp1 inp2 inp3


The program can be run in series (default) or in parallel with the python
multiprocessing package using options ``--np`` and ``-P``.
The number after the ``--np`` and ``-P`` options declares the number of processes to run in parallel.

The parameter files are read line by line to fill a python dictionary variable
mapping specific parameter names with their values.
The parameter names are the first column in the file and the parameter values are the second column.
For parameters consisting of lists (e.g. GRACE/GRACE-FO missing months),
the parameter values are separated by commas.

- Column 1: parameter name (such as ``LMAX``)
- Column 2: parameter value (e.g. ``60``)
- Column 3: comments (which are discarded)

Dataset Parameters
##################

- ``PROC``: GRACE Processing Center (CSR, GFZ, JPL, CNES)
- ``DREL``: GRACE data release for given processing center
- ``DSET``: GRACE data product (see `GRACE Data File Formats <./GRACE-Data-File-Formats.html>`_)
- ``LMIN``: minimum spherical harmonic degree (lower bound of truncation)
- ``LMAX``: maximum spherical harmonic degree (upper bound of truncation)
- ``MMAX``: maximum spherical harmonic order (None if LMAX)
- ``START``: first month to be analyzed
- ``END``: last month to be analyzed
- ``MISSING``: GRACE/GRACE-FO months that are not be analyzed (see available GRACE/GRACE-FO months)
- ``RAD``: Gaussian smoothing radius in km (`Jekeli, 1981 <http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf>`_)
- ``DESTRIPE``: filter coefficients using destriping procedure (`Swenson et al., 2006 <https://doi.org/10.1029/2005GL025285>`_)
- ``SLR_C20``: replace *C*\ :sub:`20` coefficients with values from Satellite Laser Ranging (SLR)

     * `None`: use original values
     * ``'CSR'``: use values from CSR (TN-07, TN-09, TN-11)
     * ``'GSFC'``: use values from GSFC (TN-14)

- ``SLR_C30``: replace *C*\ :sub:`30` coefficients with values from Satellite Laser Ranging (SLR)

     * `None`: use original values
     * ``'CSR'``: use values from CSR (5x5 with 6,1)
     * ``'GSFC'``: use values from GSFC (TN-14)
     * ``'LARES'``: use filtered values from CSR (John Ries)

- ``DEG1``: account for variations in geocenter with specified values

     * `None`
     * ``'Tellus'``: GRACE/GRACE-FO TN-13 coefficients from PO.DAAC
     * ``'SLR'``: satellite laser ranging coefficients from CSR
     * ``'SLF'``: Sutterley and Velicogna coefficients, Remote Sensing (2019)
     * ``'Swenson'``: GRACE-derived coefficients from Sean Swenson

- ``MODEL_DEG1``: use a least-squares regression model to predict geocenter values where unavailable
- ``GIA``: GIA model type

     * `None`
     * ``'IJ05-R2'``: `Ivins R2 GIA Models <https://doi.org/10.1002/jgrb.50208>`_
     * ``'W12a'``: `Whitehouse GIA Models <https://doi.org/10.1111/j.1365-246X.2012.05557.x>`_
     * ``'SM09'``: `Simpson/Milne GIA Models <https://doi.org/10.1029/2010JB007776>`_
     * ``'ICE6G'``: `ICE-6G GIA Models <https://doi.org/10.1002/2014JB011176>`_
     * ``'Wu10'``: `Wu (2010) GIA Correction <https://doi.org/10.1038/ngeo938>`_
     * ``'AW13-ICE6G'``: `Geruo A ICE-6G GIA Models <https://doi.org/10.1093/gji/ggs030>`_
     * ``'Caron'``: `Caron JPL GIA Assimilation <https://doi.org/10.1002/2017GL076644>`_
     * ``'ICE6G-D'``: `ICE-6G Version-D GIA Models <https://doi.org/10.1002/2016JB013844>`_
     * ``'netCDF4'``: reformatted GIA in netCDF4 format
     * ``'HDF5'``: reformatted GIA in HDF5 format

- ``GIA_FILE``: path to specific GIA file to be read
- ``DATAFORM``: input data format and output data format

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``

- ``MEAN``: Remove a mean field to isolate the time-variable gravity field
- ``MEAN_FILE``: use a file to remove as static field (default: mean of imported month)
- ``MEANFORM``: Data format for input ``MEAN_FILE``

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
     * ``'gfc'``

- ``DIRECTORY``: Directory to output data
- ``REMOVE_FILE``: Remove sets of spherical harmonics (can be multiple files)
- ``REMOVEFORM``: Data format for input ``REMOVE_FILE`` (can be a single value for a uniform type or values for each file)

     * ``'ascii'``
     * ``'netCDF4'``
     * ``'HDF5'``
     * ``'index'``: index file containing monthly files in ``DATAFORM``

- ``REDISTRIBUTE_REMOVED``: Redistribute total mass of removed harmonics over the ocean
- ``POLE_TIDE``: correct GSM *C*\ :sub:`21` and *S*\ :sub:`21` for pole tide (`Wahr et al., 2015 <https://doi.org/10.1002/2015JB011986>`_)
- ``ATM``: correct Atmosphere with `ECMWF "jump" corrections <https://doi.org/10.1093/gji/ggv276>`_
- ``UNITS``: Output units of the spatial fields

     * ``1``: Equivalent Water Thickness (cm)
     * ``2``: Geoid Height (mm)
     * ``3``: Elastic Crustal Uplift (mm)
     * ``4``: Gravitational Undulation (\ |mu|\ Gal)
     * ``5``: Equivalent surface pressure (millibar)

- ``DDEG``: spatial longitude and latitude degree spacing
- ``INTERVAL``: determines the spatial field degree interval

     * ``1``: (90:-90,0:360)
     * ``2``: (degree interval/2)

- ``FILENAME``: Start of the output filename

.. |beta|    unicode:: U+03B2 .. GREEK SMALL LETTER BETA
.. |mu|      unicode:: U+03BC .. GREEK SMALL LETTER MU
.. |mlt|     unicode:: U+226A .. MUCH LESS-THAN