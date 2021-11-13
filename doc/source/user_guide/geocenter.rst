============
geocenter.py
============

Data class for reading and processing geocenter data

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/geocenter.py

General Attributes and Methods
==============================

.. class:: geocenter(object)

    .. attribute:: object.C10

        cosine spherical harmonics of degree 1 and order 0

    .. attribute:: object.C11

        cosine spherical harmonics of degree 1 and order 1

    .. attribute:: object.S11

        sine spherical harmonics of degree 1 and order 1

    .. attribute:: object.X

        X-component of Cartesian geocenter coordinates

    .. attribute:: object.Y

        Y-component of Cartesian geocenter coordinates

    .. attribute:: object.Z

        Z-component of Cartesian geocenter coordinates

    .. attribute:: object.time

        time variable of the spherical harmonics

    .. attribute:: object.month

        GRACE/GRACE-FO months variable of the spherical harmonics

    .. attribute:: object.radius

        Average Radius of the Earth [mm]


    .. method:: object.case_insensitive_filename(filename)

        Searches a directory for a filename without case dependence


    .. method:: object.from_AOD1B(release, calendar_year, calendar_month)

        Reads monthly non-tidal ocean and atmospheric variation geocenter files

        Arguments:

            ``release``: GRACE/GRACE-FO/Swarm data release for dealiasing product

            ``calendar_year``: calendar year of data
            
            ``calendar_month``: calendar month of data


    .. method:: object.from_gravis(geocenter_file, **kwargs)

        - Reads `monthly geocenter coefficients <ftp://isdcftp.gfz-potsdam.de/grace/GravIS/GFZ/Level-2B/aux_data/GRAVIS-2B_GFZOP_GEOCENTER_0002.dat>`_ from GFZ GravIS [Dahle2019]_

        - Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_ [Sun2016]_

        Arguments:

            ``geocenter_file``: degree 1 file

        Keyword arguments:

            ``header``: file contains header text to be skipped


    .. method:: object.from_SLR(geocenter_file, **kwargs)

        - Reads monthly geocenter spherical harmonic data files from [satellite laser ranging (SLR)](ftp://ftp.csr.utexas.edu/pub/slr/geocenter/) provided by CSR

        - Can be corrected for non-tidal ocean and atmospheric variation

        Arguments:

            ``geocenter_file``: Satellite Laser Ranging file

        Keyword arguments:

            ``AOD``: remove Atmospheric and Oceanic Dealiasing products

            ``release``: GRACE/GRACE-FO/Swarm data release for AOD

            ``header``: rows of data to skip when importing data

            ``columns``: column names of ascii file


    .. method:: object.from_UCI(geocenter_file, **kwargs)

        - Reads geocenter file and extracts dates and spherical harmonic data [Sutterley2019]_

        - Estimates calculated using GRACE/GRACE-FO measurements and Ocean Models of Degree 1 [Swenson2008]_ [Sutterley2019]_

        Arguments:

            ``geocenter_file``: degree 1 file


    .. method:: object.from_swenson(geocenter_file, **kwargs)

        - Reads `monthly geocenter coefficients <https://github.com/swensosc/GRACE_Tiles/blob/master/ancillary_data/gad_gsm.rl05.txt>`_ in units mm w.e

        - Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_

        Arguments:

            ``geocenter_file``: degree 1 file

        Keyword arguments:

            ``header``: file contains header text to be skipped


    .. method:: object.from_tellus(geocenter_file, **kwargs)

        - Reads monthly geocenter spherical harmonic data files from `GRACE Tellus Technical Notes <https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L2/degree_1>`_
  
        - Estimates calculated using GRACE measurements and Ocean Models of Degree 1 [Swenson2008]_ [Sun2016]_

        Arguments:

            ``geocenter_file``: degree 1 file
            
                * CSR: ``TN-13_GEOC_CSR_RL06.txt``

                * GFZ: ``TN-13_GEOC_GFZ_RL06.txt``
                
                * JPL: ``TN-13_GEOC_JPL_RL06.txt``

        Keyword arguments:

            ``header``: file contains header text to be skipped

            ``JPL``: use JPL TN-13 geocenter files calculated following [Sun2016]_


    .. method:: object.copy(**kwargs)

        Copy a geocenter object to a new geocenter object

        Keyword arguments:

            ``fields``: default keys in geocenter object


    .. method:: object.from_dict(temp, **kwargs):

        Convert a dictionary object to a geocenter object

        Arguments:

            ``temp``: dictionary object to be converted

        Keyword arguments:

            ``fields``: default keys in dictionary


    .. method:: object.from_harmonics(temp, **kwargs):

        Convert a harmonics object to a geocenter object

        Arguments:

            ``temp``: harmonics object to be converted

        Keyword arguments:

            ``fields``: default keys in harmonics object


    .. method:: object.from_matrix(clm, slm):

        Converts spherical harmonic matrices to a geocenter object

        Arguments:

            ``clm``: cosine spherical harmonics of degree 1

            ``slm``: sine spherical harmonics of degree 1


    .. method:: object.to_dict():

        Convert a geocenter object to a dictionary object

        Keyword arguments:

            ``fields``: default attributes in geocenter object


    .. method:: object.to_matrix():

        Converts a geocenter object to spherical harmonic matrices


    .. method:: object.to_cartesian(kl=0.0):

        Converts normalized spherical harmonics to Cartesian geocenter variations

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1


    .. method:: object.to_cmwe(kl=0.0):

        Converts normalized spherical harmonics to centimeters water equivalent

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1     


    .. method:: object.to_mmwe(kl=0.0):

        Converts normalized spherical harmonics to millimeters water equivalent

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1


    .. method:: object.from_cartesian(kl=0.0):

        Converts Cartesian geocenter variations to normalized spherical harmonics

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1


    .. method:: object.from_cmwe(kl=0.0):

        Normalizes spherical harmonics from centimeters water equivalent (cmwe)

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1
       
        
    .. method:: object.from_mmwe(kl=0.0):

        Normalizes spherical harmonics from millimeters water equivalent (mmwe)

        Keyword arguments:

            ``kl``: gravitational load love number of degree 1


    .. method:: object.mean(apply=False, indices=Ellipsis)

        Compute mean gravitational field and remove from data if specified

        Keyword arguments:

            ``apply``: remove the mean field from the input harmonics

            ``indices``: of input harmonics object to compute mean


    .. method:: object.add(self, temp):

        Add two geocenter objects

        Arguments:

            ``temp``: geocenter object to be added

    .. method:: object.subtract(self, temp):

        Subtract one geocenter object from another

        Arguments:
        
            ``temp``: geocenter object to be subtracted

    .. method:: object.multiply(self, temp):

        Multiply two geocenter objects

        Arguments:
        
            ``temp``: geocenter object to be multiplied


    .. method:: object.divide(self, temp):

        Divide one geocenter object from another

        Arguments:
        
            ``temp``: geocenter object to be divided


    .. method:: object.scale(self, var):

        Multiply a geocenter object by a constant

        Arguments:
        
            ``var``: scalar value to which the geocenter object will be multiplied


    .. method:: object.power(self, power):

        Raise a geocenter object to a power

        Arguments:
        
            ``power``: power to which the geocenter object will be raised


References
##########

.. [Dahle2019] Dahle and Murboeck, "Post-processed GRACE/GRACE-FO Geopotential GSM Coefficients GFZ RL06 (Level-2B Product)." V. 0002. *GFZ Data Services*, (2019). `doi: 10.5880/GFZ.GRAVIS_06_L2B <http://doi.org/10.5880/GFZ.GRAVIS_06_L2B>`_

.. [Sun2016] Y. Sun, P. Ditmar, and R. Riva, "Observed changes in the Earth's dynamic oblateness from GRACE data and geophysical models", *Journal of Geodesy*, 90(1), 81--89, (2016). `doi: 10.1007/s00190-015-0852-y <https://doi.org/10.1007/s00190-015-0852-y>`_

.. [Sutterley2019] T. C. Sutterley and I. Velicogna, "Improved Estimates of Geocenter Variability from Time-Variable Gravity and Ocean Model Outputs", *Remote Sensing*, 11(18), 2108, (2019). `doi: 10.3390/rs11182108 <https://doi.org/10.3390/rs11182108>`_

.. [Swenson2008] S. Swenson, D. Chambers, and J. Wahr, "Estimating geocenter variations from a combination of GRACE and ocean model output", *Journal of Geophysical Research: Solid Earth*, 113(B08410), (2008). `doi: 10.1029/2007JB005338 <https://doi.org/10.1029/2007JB005338>`_
