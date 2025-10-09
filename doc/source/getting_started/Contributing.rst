=======================
Contribution Guidelines
=======================

``gravity-toolkit`` is an open source project.
We welcome any help in maintaining and developing the software and documentation.
Anyone at *any career stage and with any level of coding experience* can contribute towards the development of ``gravity-toolkit``.
Please read our `code of conduct <./Code-of-Conduct.html>`_ before contributing to ``gravity-toolkit`` development.
You will be recognized for your work by being listed as one of the `project contributors <../project/Contributors.html>`_.

.. note::

    If you have found a problem in ``gravity-toolkit``, or you would like to suggest an improvement or modification,
    please submit a `GitHub issue <https://github.com/tsutterley/gravity-toolkit/issues>`_ and we will get back to you.

Ways to Contribute
------------------

1) Fixing typographical or coding errors
2) Submitting bug reports or feature requests through the use of `GitHub issues <https://github.com/tsutterley/gravity-toolkit/issues>`_
3) Improving documentation and testing
4) Sharing use cases and examples (such as `Jupyter Notebooks <../user_guide/Examples.html>`_)
5) Providing code for everyone to use

Requesting a Feature
--------------------
Check the `project issues tab <https://github.com/tsutterley/gravity-toolkit/issues>`_ to see if the feature has already been suggested.
If not, please submit a new issue describing your requested feature or enhancement .
Please give your feature request both a clear title and description.
Let us know if this is something you would like to contribute to ``gravity-toolkit`` in your description as well.

Reporting a Bug
---------------
Check the `project issues tab <https://github.com/tsutterley/gravity-toolkit/issues>`_ to see if the problem has already been reported.
If not, *please* submit a new issue so that we are made aware of the problem.
Please provide as much detail as possible when writing the description of your bug report.
Providing information and examples will help us resolve issues faster.

Contributing Code
-----------------
We follow a standard Forking Workflow for code changes and additions.
Submitted code goes through the pull request process for `continuous integration (CI) testing <../project/Testing.html#continuous-integration>`_ and comments.

General Guidelines
^^^^^^^^^^^^^^^^^^

- Make each pull request as small and simple as possible
- `Commit messages should be clear and describe the changes <./Contributing.html#semantic-commit-messages>`_
- Larger changes should be broken down into their basic components and integrated separately
- Bug fixes should be their own pull requests with an associated `GitHub issue <https://github.com/tsutterley/gravity-toolkit/issues>`_
- Write a descriptive pull request message with a clear title
- Be patient as reviews of pull requests take time

Steps to Contribute
^^^^^^^^^^^^^^^^^^^

1) Fork the repository to your personal GitHub account by clicking the "Fork" button on the project `main page <https://github.com/tsutterley/gravity-toolkit>`_.  This creates your own server-side copy of the repository.
2) Either by cloning to your local system or working in `GitHub Codespaces <https://github.com/features/codespaces>`_, create a work environment to make your changes.
3) Add your fork as the ``origin`` remote and the original project repository as the ``upstream`` remote.  While this step isn't a necessary, it allows you to keep your fork up to date in the future.
4) Create a new branch to do your work.
5) Make your changes on the new branch and add yourself to the list of project `contributors <https://github.com/tsutterley/gravity-toolkit/blob/main/CONTRIBUTORS.md>`_.
6) Push your work to GitHub under your fork of the project.
7) Submit a `Pull Request <https://github.com/tsutterley/gravity-toolkit/pulls>`_ from your forked branch to the project repository.

Adding Examples
^^^^^^^^^^^^^^^
Examples may be in the form of executable scripts or interactive `Jupyter Notebooks <../user_guide/Examples.html>`_.
Fully working (but unrendered) examples should be submitted with the same steps as above.

Semantic Commit Messages
^^^^^^^^^^^^^^^^^^^^^^^^

Please follow the `Conventional Commits <https://www.conventionalcommits.org/>`_ specification for your commit messages to help organize the pull requests:

.. code-block:: bash

    <type>: <subject>

    [optional message body]

where ``<type>`` is one of the following:

- ``feat``: adding new features or programs
- ``fix``: fixing bugs or problems
- ``docs``: changing the documentation
- ``style``: changing the line order or adding comments
- ``refactor``: changing the names of variables or programs
- ``ci``: changing the `continuous integration <../project/Testing.html#continuous-integration>`_ configuration files or scripts
- ``test``: adding or updating `continuous integration tests <../project/Testing.html#continuous-integration>`_
