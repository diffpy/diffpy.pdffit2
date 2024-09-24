|Icon| |title|_
===============

.. |title| replace:: diffpy.pdffit2
.. _title: https://diffpy.github.io/diffpy.pdffit2

.. |Icon| image:: https://avatars.githubusercontent.com/diffpy
        :target: https://diffpy.github.io/diffpy.pdffit2
        :height: 100px

|PyPi| |Forge| |PythonVersion| |PR|

|CI| |Codecov| |Black| |Tracking|

.. |Black| image:: https://img.shields.io/badge/code_style-black-black
        :target: https://github.com/psf/black

.. |CI| image:: https://github.com/diffpy/diffpy.pdffit2/actions/workflows/matrix-and-codecov-on-merge-to-main.yml/badge.svg
        :target: https://github.com/diffpy/diffpy.pdffit2/actions/workflows/matrix-and-codecov-on-merge-to-main.yml

.. |Codecov| image:: https://codecov.io/gh/diffpy/diffpy.pdffit2/branch/main/graph/badge.svg
        :target: https://codecov.io/gh/diffpy/diffpy.pdffit2

.. |Forge| image:: https://img.shields.io/conda/vn/conda-forge/diffpy.pdffit2
        :target: https://anaconda.org/conda-forge/diffpy.pdffit2

.. |PR| image:: https://img.shields.io/badge/PR-Welcome-29ab47ff

.. |PyPi| image:: https://img.shields.io/pypi/v/diffpy.pdffit2
        :target: https://pypi.org/project/diffpy.pdffit2/

.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/diffpy.pdffit2
        :target: https://pypi.org/project/diffpy.pdffit2/

.. |Tracking| image:: https://img.shields.io/badge/issue_tracking-github-blue
        :target: https://github.com/diffpy/diffpy.pdffit2/issues

PDFfit2 - space structure refinement to atomic pair distribution function

The diffpy.pdffit2 package provides functions for calculation and
refinement of atomic Pair Distribution Function (PDF) from crystal
structure model.  It is used as a computational engine by PDFgui.  All
refinements possible in PDFgui can be done with diffpy.pdffit2,
although less conveniently and with a fair knowledge of Python.
The package includes an extension for the interactive `IPython
<http://ipython.org>`__ shell, which tries to mimic the old PDFFIT
program.  To start IPython with this extension and also with plotting
functions enabled, use ::

   ipython --ext=diffpy.pdffit2.ipy_ext --pylab

The IPython extension is suitable for interactive use, however
refinement scripts should be preferably written as a standard
Python code.  This is more reliable and needs only a few extra
statements.

To learn more about diffpy.pdffit2 library, see the examples directory
included in this distribution or the API documentation at
http://www.diffpy.org/doc/pdffit2.

For more information about the diffpy.pdffit2 library, please consult our `online documentation <https://diffpy.github.io/diffpy.pdffit2>`_.

Citation
--------

If you use diffpy.pdffit2 in a scientific publication, we would like you to cite the following paper:

        C\. L. Farrow, P. Juhás, J. W. Liu, D. Bryndin, E. S. Božin, J. Bloch, Th. Proffen
        and S. J. L. Billinge, PDFfit2 and PDFgui: computer programs for studying nanostructure
        in crystals (https://stacks.iop.org/0953-8984/19/335219), *J. Phys.: Condens. Matter*, 19, 335219 (2007)

Installation
------------

diffpy.pdffit2 requires Python 3.10 or later and
the following external software:

* ``setuptools`` - software distribution tools for Python
* ``python-dev`` - header files for interfacing Python with C
* ``GSL`` - GNU Scientific Library for C
* ``g++`` - GNU C++ compiler
* ``diffpy.structure`` - simple storage and manipulation of atomic
  structures, https://github.com/diffpy/diffpy.structure

----

The preferred method is to use `Miniconda Python
<https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>`_
and install from the "conda-forge" channel of Conda packages.

To add "conda-forge" to the conda channels, run the following in a terminal. ::

        conda config --add channels conda-forge

We want to install our packages in a suitable conda environment.
The following creates and activates a new environment named ``diffpy.pdffit2_env`` ::

        conda create -n diffpy.pdffit2_env python=3
        conda activate diffpy.pdffit2_env

Then, to fully install ``diffpy.pdffit2`` in our active environment, run ::

        conda install diffpy.pdffit2

Another option is to use ``pip`` to download and install the latest release from
`Python Package Index <https://pypi.python.org>`_.
To install using ``pip`` into your ``diffpy.pdffit2_env`` environment, type ::

        pip install diffpy.pdffit2

If you prefer to install from sources, after installing the dependencies, obtain the source archive from
`GitHub <https://github.com/diffpy/diffpy.pdffit2/>`_. Once installed, ``cd`` into your ``diffpy.pdffit2`` directory
and run the following ::

        pip install .

Support and Contribute
----------------------

PDFfit2 is not developed anymore and is only maintained due to its
status of a sole computational engine for PDFgui.  We don't expect any
major developments to the code beyond simple bug fixes and compatibility
features.  The source code to PDFfit2 is available in a git repository
at https://github.com/diffpy/diffpy.pdffit2.

For an actively developed codes for PDF simulations see the
DiffPy-CMI framework at http://www.diffpy.org.

----

`Diffpy user group <https://groups.google.com/g/diffpy-users>`_ is the discussion forum for general questions and discussions about the use of diffpy.pdffit2. Please join the diffpy.pdffit2 users community by joining the Google group. The diffpy.pdffit2 project welcomes your expertise and enthusiasm!

If you see a bug or want to request a feature, please `report it as an issue <https://github.com/diffpy/diffpy.pdffit2/issues>`_ and/or `submit a fix as a PR <https://github.com/diffpy/diffpy.pdffit2/pulls>`_. You can also post it to the `Diffpy user group <https://groups.google.com/g/diffpy-users>`_.

Feel free to fork the project and contribute. To install diffpy.pdffit2
in a development mode, with its sources being directly used by Python
rather than copied to a package directory, use the following in the root
directory ::

        pip install -e .

To ensure code quality and to prevent accidental commits into the default branch, please set up the use of our pre-commit
hooks.

1. Install pre-commit in your working environment by running ``conda install pre-commit``.

2. Initialize pre-commit (one time only) ``pre-commit install``.

Thereafter your code will be linted by black and isort and checked against flake8 before you can commit.
If it fails by black or isort, just rerun and it should pass (black and isort will modify the files so should
pass after they are modified). If the flake8 test fails please see the error messages and fix them manually before
trying to commit again.

Improvements and fixes are always appreciated.

Before contribuing, please read our `Code of Conduct <https://github.com/diffpy/diffpy.pdffit2/blob/main/CODE_OF_CONDUCT.rst>`_.

Contact
-------

For more information on diffpy.pdffit2 please visit the project `web-page <https://diffpy.github.io/>`_ or email Prof. Simon Billinge at sb2896@columbia.edu.
