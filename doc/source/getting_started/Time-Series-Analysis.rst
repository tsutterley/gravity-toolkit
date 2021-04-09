====================
Time Series Analysis
====================

Least-squares mascons are a method of extracting a regional signal from the
GRACE/GRACE-FO spherical harmonic data.
The procedure was outlined in procedure outlined in
[Tiwari2009]_ and [Jacob2012]_.
Least-squares mascons can be considered a post-processing technique for
analyzing the GRACE/GRACE-FO data.
The technique calculates the scaling factor between an input kernel and the
GRACE/GRACE-FO data for a given month.
For example, for a uniform kernel equivalent to 1 cm w.e.,
if GRACE/GRACE-FO measures 6 cm w.e. in the region, then the scaling factor would be 6.

Ideally we would want as small of kernels as possible to get as geophysically
relevant as possible (perfect fit to region shape).
However, GRACE/GRACE-FO data is limited to a specific spherical harmonic range
with noise increasing with higher degree and order (data-to-noise ratio decreases).
If the kernels are too small, there will be a higher degree of ringing when
truncated (Gibbs phenomenon) and require more information at the higher degrees
and orders to distinguish from the adjacent kernels.
Thus there is a balance between wanting kernels as geophysically relevant as
possible with wanting idealistic kernels which minimize ringing and are resolvable.
Getting the kernels "just right" in order to isolate regions of interest takes some time.

The set of least-squares mascon programs have been used in [Velicogna2014]_
and other publications for regional time series analysis.
The ``calc_mascon.py`` program additionally calculates the GRACE/GRACE-FO error
harmonics following [Wahr2006]_.

The ``calc_mascon.py`` program will output a text file of the time series for each mascon
(format: GRACE/GRACE-FO month, mid-month date in decimal-year format,
estimated monthly mass anomaly [Gt], estimated monthly error [Gt],
mascon area [km\ :sup:`2`]).

Parameter Files
###############

The ``calc_mascon.py`` program works by accepting parameter file system arguments
(``sys.argv``) listed after the program.
These parameter files can be ran individually or in a series.

.. code-block:: bash

	python calc_mascon.py inp1 inp2 inp3
	python calc_mascon.py --np=2 inp1 inp2 inp3
	python calc_mascon.py -P 2 inp1 inp2 inp3

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
- ``MMAX``: maximum spherical harmonic order (None if ``LMAX``)
- ``START``: first month to be analyzed
- ``END``: last month to be analyzed
- ``MISSING``: GRACE/GRACE-FO months that are not be analyzed (see available GRACE/GRACE-FO months)
- ``RAD``: Gaussian smoothing radius in km [Jekeli1981]_
- ``DESTRIPE``: filter coefficients using destriping procedure [Swenson2006]_
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

- ``DATAFORM``: input data format for mascon files and files to be removed from the GRACE/GRACE-FO data

	* ``'ascii'``
	* ``'netCDF4'``
	* ``'HDF5'``

- ``DIRECTORY``: Directory to output data (will create directory if non-existent)
- ``MASCON_INDEX``: file index listing the full path to each mascon file to fit to the GRACE data
- ``FIT_METHOD``: method of fitting mascons coefficients

	* ``1``: convert coefficients to mass
	* ``2``: keep coefficients as normalized geoid

- ``MEAN``: Remove a mean field to isolate the time-variable gravity field
- ``MEAN_FILE``: use a file to remove as static field (default: mean of imported month)
- ``MEANFORM``: Data format for input ``MEAN_FILE``

	* ``'ascii'``
	* ``'netCDF4'``
	* ``'HDF5'``
	* ``'gfc'``

- ``REMOVE_FILE``: Remove sets of spherical harmonics (can be multiple files)
- ``REMOVEFORM``: Data format for input ``REMOVE_FILE`` (can be a single value for a uniform type or values for each file)

	* ``'ascii'``
	* ``'netCDF4'``
	* ``'HDF5'``
	* ``'index'``: index file containing monthly files in ``DATAFORM``

- ``REDISTRIBUTE_REMOVED``: Redistribute total mass of removed harmonics over the ocean
- ``MASCON_OCEAN``: remove uniformly distributed mascon mass over ocean
- ``RECONSTRUCT``: remove the reconstructed time series for a region to get the statistical leakage
- ``POLE_TIDE``: correct GSM *C*\ :sub:`21` and *S*\ :sub:`21` for pole tide [Wahr2015]_
- ``ATM``: correct Atmosphere with ECMWF "jump" corrections [Fagiolini2015]_

References
##########

.. [Fagiolini2015] E. Fagiolini, F. Flechtner, M. Horwath, and H. Dobslaw, "Correction of inconsistencies in ECMWF's operational analysis data during de-aliasing of GRACE gravity models", *Geophysical Journal International*, 202(3), 2150--2158, (2015). `doi: 10.1093/gji/ggv276 <https://doi.org/10.1093/gji/ggv276>`_

.. [Jacob2012] T. Jacob, J. Wahr, W. T. Pfeffer, and S. Swenson, "Recent contributions of glaciers and ice caps to sea level rise", *Nature*, 482, 514--518, (2012). `doi: 10.1038/nature10847 <https://doi.org/10.1038/nature10847>`_

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Tiwari2009] V. M. Tiwari, J. Wahr, and S. Swenson, "Dwindling groundwater resources in northern India, from satellite gravity observations", *Geophysical Research Letters*, 36(L18401), (2009). `doi: 10.1029/2009GL039401 <https://doi.org/10.1029/2009GL039401>`_

.. [Velicogna2014] I. Velicogna, T. C. Sutterley, and M. R. van den Broeke, "Regional acceleration in ice mass loss from Greenland and Antarctica using GRACE time‐variable gravity data", *Geophysical Research Letters*, 119, 8130--8137, (2014). `doi: 10.1002/2014GL061052 <https://doi.org/10.1002/2014GL061052>`_

.. [Wahr2006] J. Wahr, S. Swenson, and I. Velicogna, "Accuracy of GRACE mass estimates", Geophysical Research Letters, 33(L06401), (2006). `doi: 10.1029/2005GL025305 <https://doi.org/10.1029/2005GL025305>`_

.. [Wahr2015] J. Wahr, R. S. Nerem, and S. V. Bettadpur, "The pole tide and its effect on GRACE time‐variable gravity measurements: Implications for estimates of surface mass variations". *Journal of Geophysical Research: Solid Earth*, 120, 4597--4615. `doi: 10.1002/2015JB011986 <https://doi.org/10.1002/2015JB011986>`_
