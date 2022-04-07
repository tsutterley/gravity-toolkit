====================
piecewise_regress.py
====================

- Fits a synthetic signal to data over a time period by ordinary or weighted least-squares for breakpoint analysis

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.piecewise_regress import piecewise_regress
    tsbeta = piecewise_regress(t_in, d_in, BREAKPOINT=len(t_in)//2, CYCLES=[0.5,1.0])

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/piecewise_regress.py

.. autofunction:: gravity_toolkit.piecewise_regress
