==========================
calc_sensitivity_kernel.py
==========================

- Calculates spatial sensitivity kernels through a least-squares mascon procedure following :cite:p:`Tiwari:2009bx` :cite:p:`Jacob:2012gv`

`Source code`__

.. __: https://github.com/tsutterley/gravity-toolkit/blob/main/scripts/calc_sensitivity_kernel.py

Calling Sequence
################

.. argparse::
    :filename: calc_sensitivity_kernel.py
    :func: arguments
    :prog: calc_sensitivity_kernel.py
    :nodescription:
    :nodefault:

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM :cite:p:`Han:1995go`
        * ``1``: Gegout (2005) values from PREM :cite:p:`Gegout:2010gc`
        * ``2``: Wang et al. (2012) values from PREM :cite:p:`Wang:2012gc`
        * ``3``: Wang et al. (2012) values from PREM with hard sediment :cite:p:`Wang:2012gc`
        * ``4``: Wang et al. (2012) values from PREM with soft sediment :cite:p:`Wang:2012gc`

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth

    --fit-method : @after
        * ``1``: mass coefficients
        * ``2``: geoid coefficients

    --solver -s : @replace
        Least squares solver for sensitivity kernels

        * ``'inv'``: matrix inversion
        * ``'lstsq'``: least squares solution
        * ``'gelsy'``: complete orthogonal factorization solution
        * ``'gelss'``: singular value decomposition (SVD) solution
        * ``'gelsd'``: singular value decomposition (SVD) solution with a divide and conquer method
