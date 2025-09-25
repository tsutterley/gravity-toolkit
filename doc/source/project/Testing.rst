=======
Testing
=======

``gravity-toolkit`` uses the ``pytest`` framework to run tests and verify outputs.

Running the Test Suite
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    pytest --directory <path_to_tide_models> test/

The test suite is run in verbose mode as a default.

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^
We use `GitHub Actions <https://github.com/tsutterley/gravity-toolkit/actions>`_ continuous integration (CI) services to build and test the project on Linux (``ubuntu-latest``), Mac (``macos-latest``) and Windows (``windows-latest``) Operating Systems.
The configuration files for this service are in the `GitHub workflows <https://github.com/tsutterley/gravity-toolkit/tree/main/.github/workflows>`_ directory.
The workflows use ``mamba`` and the `environment.yml <https://github.com/tsutterley/gravity-toolkit/blob/main/environment.yml>`_ file to install the required dependencies and build the environment.

The GitHub Actions jobs include:

* Running `flake8 <https://flake8.pycqa.org/en/latest/>`_ to check the code for style and compilation errors
* Running the test suite on multiple combinations of OS and Python version
* Uploading source and wheel distributions to `PyPI <https://pypi.org/project/gravity-toolkit/>`_ (on releases)
