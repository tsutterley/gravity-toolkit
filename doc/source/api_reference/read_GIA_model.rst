==============
read_GIA_model
==============

- Reads Glacial Isostatic Adjustment (GIA) files for given modeling group formats
- Outputs spherical harmonics for the GIA rates and the GIA model parameters
- Can also output fully normalized harmonics to netCDF4 or HDF5 formats

Calling Sequence
################

.. code-block:: python

   from gravity_toolkit.read_GIA_model import read_GIA_model
   GIA_Ylms = read_GIA_model('Stokes.R2_65_.2_1.5_L120',GIA='IJ05-R2',LMAX=60)

`Source code`__

   .. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/read_GIA_model.py

.. autofunction:: gravity_toolkit.read_GIA_model

.. autoclass:: gravity_toolkit.read_GIA_model.gia
   :members:
