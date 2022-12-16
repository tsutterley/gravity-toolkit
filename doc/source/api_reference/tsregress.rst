=========
tsregress
=========

- Fits a synthetic signal to data over a time period by ordinary or weighted least-squares

Calling Sequence
################

.. code-block:: python

   from gravity_toolkit.tsregress import tsregress
   tsbeta = tsregress(t_in, d_in, ORDER=1, CYCLES=[0.5,1.0], CONF=0.95)

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/tsregress.py

.. autofunction:: gravity_toolkit.tsregress
