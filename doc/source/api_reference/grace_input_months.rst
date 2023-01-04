==================
grace_input_months
==================

- Reads GRACE/GRACE-FO/Swarm files for a specified spherical harmonic degree and order and for a specified date range

Calling Sequence
################

.. code-block:: python

      from gravity_toolkit.grace_input_months import grace_input_months
      GRACE_Ylms = grace_input_months(base_dir, PROC, DREL, DSET, LMAX,
        start_mon, end_mon, missing, SLR_C20, DEG1, SLR_C30=SLR_C30)

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/gravity_toolkit/grace_input_months.py

.. autofunction:: gravity_toolkit.grace_input_months

.. autofunction:: gravity_toolkit.grace_input_months.read_ecmwf_corrections

.. autofunction:: gravity_toolkit.grace_input_months.regress_model