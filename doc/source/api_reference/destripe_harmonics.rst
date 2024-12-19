==================
destripe_harmonics
==================

- Filters spherical harmonic coefficients for correlated "striping" errors following :cite:p:`Swenson:2006hu`

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.destripe_harmonics import destripe_harmonics
    Ylms = destripe_harmonics(clm,slm,LMAX=60)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/destripe_harmonics.py

.. autofunction:: gravity_toolkit.destripe_harmonics
