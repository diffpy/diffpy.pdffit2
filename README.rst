.. image:: https://travis-ci.org/diffpy/diffpy.pdffit2.svg?branch=master
   :target: https://travis-ci.org/diffpy/diffpy.pdffit2

.. image:: https://codecov.io/gh/diffpy/diffpy.pdffit2/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/diffpy/diffpy.pdffit2


PDFfit2
========================================================================

Real space structure refinement to atomic pair distribution function.

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


REQUIREMENTS
------------------------------------------------------------------------

diffpy.pdffit2 requires Python 3.7 or later or 2.7 and
the following external software:

* ``setuptools`` - software distribution tools for Python
* ``python-dev`` - header files for interfacing Python with C
* ``GSL`` - GNU Scientific Library for C
* ``g++`` - GNU C++ compiler
* ``diffpy.structure`` - simple storage and manipulation of atomic
  structures, https://github.com/diffpy/diffpy.structure

We recommend to use `Anaconda Python <https://www.anaconda.com/distribution>`_
as it allows to install all software dependencies together with
PDFfit2.  For other Python distributions it is necessary to
install the required software separately.  As an example, on Ubuntu
Linux some of the required software can be installed using ::

   sudo apt-get install \
      python-setuptools python-dev libgsl0-dev build-essential


INSTALLATION
------------------------------------------------------------------------

The preferred method is to use Anaconda Python and install from the
"conda-forge" channel of Anaconda packages ::

   conda install -c conda-forge diffpy.pdffit2

If you don't use Anaconda or prefer to install from sources, make
sure the required software is in place and run ::

   python setup.py install

By default the files get installed to standard system directories,
which may require the use of ``sudo`` for write permissions.  If
administrator (root) access is not available, consult the output from
``python setup.py install --help`` for options to install as a regular
user to other locations.  The installation integrity can be
verified by changing to the HOME directory and running ::

   python -m diffpy.pdffit2.tests.rundeps

Anaconda Python allows to later update PDFfit2 using ::

   conda update diffpy.pdffit2

With other Python distributions use the easy_install program to
upgraded to the latest version ::

   easy_install --upgrade diffpy.pdffit2


DEVELOPMENT
------------------------------------------------------------------------

PDFfit2 is not developed anymore and is only maintained due to its
status of a sole computational engine for PDFgui.  We don't expect any
major developments to the code beyond simple bug fixes and compatibility
features.  The source code to PDFfit2 is available in a git repository
at https://github.com/diffpy/diffpy.pdffit2.

For an actively developed codes for PDF simulations see the
DiffPy-CMI framework at http://www.diffpy.org.


INSTALLATION AND PACKAGE BUILDING VERIFICATION
-----------------------------------------------
To verify pdffit2 installed correctly and that the package is buildable across OS’s from Conda-forge, Github, and PyPi use the following procedure:
- Create a new conda env on your local machine and preset the python version to match python that the current version pdffit2 is running on.
   - Conda-forge: **Recommended**, follow the steps above for the installation process.
   - Github: Download the tar file for the newest version under Release Tab on the GitHub page, then run the setup file (same command as above):
   ::
   python setup.py install
   - PyPi: Clone the repository locally, then install the package using pip install pdffit2==<version_of_pdffit2>
      - Further documentation and installation notes found ‘here <https://pypi.org/project/diffpy.pdffit2/>’
- Note: When installing the package or running the local test, there might be an error about a specific package (i.e. numpy, etc.). Just install that package to your local env and rerun the command.
- To verify the package’s installation (irrespective of the installation method used) locally run the dependencies test (same command as provided above in installation):
   ::
   python -m diffpy.pdffit2.tests.rundeps
- Some helpful installation notes
   - If you a gsl fatal error when trying to install on Linux, make sure you have installed GNU scientific library. If not it is not instaled use:
   ::
      sudo apt-get install libgsl-dev

CONTACTS
------------------------------------------------------------------------

For more information on diffpy.pdffit2 please visit the project web-page:

http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu.
