diffpy.pdffit2 - real space structure refinement to atomic PDF

The diffpy.pdffit2 package enables calculation and refinement of atomic
Pair Distribution Function (PDF) from crystal structure model.  It is
used as a computational engine by PDFgui.  All refinements possible
in PDFgui can be done with diffpy.pdffit2, it is a lower level library
and requires fair Python knowledge.  The package includes a pdffit2
command-line application, which tries to mimic the old PDFFIT program. 
The pdffit2 program launches Python interactive session with several
libraries and functions preloaded for convenience.  The pdffit2 program
is suitable for interactive use, however refinement scripts should be
rather written as normal Python; this is more reliable and needs only
few extra lines of code.

To learn more about diffpy.pdffit2 library, see the examples directory
included in this distribution or the API documentation at

    http://docs.danse.us/diffraction/diffpy.pdffit2


REQUIREMENTS

diffpy.pdffit2 requires Python2.5 and the following external software:

    setuptools  -- software distribution tools for Python
    python-dev  -- header files for interfacing Python with C
    GSL         -- GNU Scientific Library for C
    g++         -- GNU C++ compiler

On Ubuntu Linux the required software can be easily installed using
the system package manager:

    sudo aptitude install \
        python-setuptools python-dev libgsl0-dev g++

For other Linux distributions use their respective package manager; note
the packages may have slightly different names.  diffpy.pdffit2 should work on
other Unix-like operating systems and on Mac as well.  Please, search the web
for instructions how to install external dependencies on your particular
system.


INSTALLATION

Once the requirements are satisfied, diffpy.pdffit2 can be installed with

    python setup.py install

This installs the library and the pdffit2 program to the default system
directories.  See the usage info "./setup.py install --help" for options
to install under different locations.  Note that installation to non-standard
directories you may require adjustment to the PATH and PYTHONPATH environment
variables.

The Python setuptools library provides "easy_install" script, which can
be used to update diffpy.pdffit2 installation or even to perform a new
install without explicit need to download and unzip the code:

    easy_install -U diffpy.pdffit2

This checks the package repository at http://www.diffpy.org/packages/
for any newer releases of diffpy.pdffit2 and if present, it updates the
installation.  The easy_install can be also used to get in sync with the
latest development sources in the subversion repository:

    easy_install -U \
        svn://svn@danse.us/diffraction/diffraction/diffpy.pdffit2/trunk


CONTACTS

For more information on diffpy.pdffit2 please visit the project web-page:

    http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu

Last modified $Date$
