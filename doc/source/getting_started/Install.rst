============
Installation
============

``gravity-toolkit`` is available for download from the `GitHub repository <https://github.com/tsutterley/gravity-toolkit>`_,
the `Python Package Index (pypi) <https://pypi.org/project/gravity-toolkit/>`_,
and from `conda-forge <https://anaconda.org/conda-forge/gravity-toolkit>`_.


The simplest installation for most users will likely be using ``conda`` or ``mamba``:

.. code-block:: bash

    conda install -c conda-forge gravity-toolkit

``conda`` installed versions of ``gravity-toolkit`` can be upgraded to the latest stable release:

.. code-block:: bash

    conda update gravity-toolkit

To use the development repository, please fork ``gravity-toolkit`` into your own account and then clone onto your system:

.. code-block:: bash

    git clone https://github.com/tsutterley/gravity-toolkit.git

``gravity-toolkit`` can then be installed within the package directory using ``pip``:

.. code-block:: bash

    python3 -m pip install --user .

The development version of ``gravity-toolkit`` can also be installed directly from GitHub using ``pip``:

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/gravity-toolkit.git
