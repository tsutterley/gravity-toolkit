====================
gen_averaging_kernel
====================

- Generates averaging kernel coefficients which minimize the total error

Calling Sequence
################

.. code-block:: python

    from gravity_toolkit.gen_averaging_kernel import gen_averaging_kernel
    Wlms = gen_averaging_kernel(gclm,gslm,eclm,eslm,sigma,hw,UNITS=0,LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_averaging_kernel.py

.. autofunction:: gravity_toolkit.gen_averaging_kernel
