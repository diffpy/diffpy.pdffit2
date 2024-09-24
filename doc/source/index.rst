#######
|title|
#######

.. |title| replace:: diffpy.pdffit2 documentation

diffpy.pdffit2 - PDFfit2 - real space structure refinement program.

| Software version |release|.
| Last updated |today|.

The diffpy.pdffit2 package provides functions for calculation and
refinement of atomic Pair Distribution Function (PDF) from crystal
structure model.  It is used as a computational engine by PDFgui. All
refinements possible in PDFgui can be done with diffpy.pdffit2,
although less conveniently and with a fair knowledge of Python.
The package includes an extension for the interactive `IPython
<http://ipython.org>`_ shell, which tries to mimic the old PDFFIT
program.  To start IPython with this extension and also with plotting
functions enabled, use ::

   ipython --ext=diffpy.pdffit2.ipy_ext --pylab

The IPython extension is suitable for interactive use, however
refinement scripts should be preferably written as a standard
Python code.  This is more reliable and needs only a few extra
statements.

=======
Authors
=======

This code was derived from the first PDFFIT program by Thomas Proffen.
The sources were converted to C++ by Jacques Bloch and then extensively hacked,
extended and purged from most glaring bugs by Chris Farrow and Pavol Juhas.
This code is currently maintained as part of the DiffPy project to create
python modules for structure investigations from diffraction data.

The DiffPy team is located in the Billinge-group at the Applied Physics
and Applied Mathematics Department of the Columbia University in New York.
Previous significant contributors to this code were

   Pavol Juhas,  Chris Farrow, Jacques Bloch, Wenduo Zhou

For a detailed list of contributors see
https://github.com/diffpy/diffpy.pdffit2/graphs/contributors.


=========
Reference
=========

If you use this program for a scientific research that leads to publication,
we ask that you acknowledge use of the program by citing the following paper
in your publication:

   C. L. Farrow, P. Juhás, J. W. Liu, D. Bryndin, E. S. Božin, J. Bloch, Th. Proffen
   and S. J. L. Billinge, PDFfit2 and PDFgui: computer programs for studying nanostructure
   in crystals (https://stacks.iop.org/0953-8984/19/335219), *J. Phys.: Condens. Matter*, 19, 335219 (2007)

============
Installation
============

See the `README <https://github.com/diffpy/diffpy.pdffit2#installation>`_
file included with the distribution.

=================
Table of contents
=================
.. toctree::
   :titlesonly:

   license
   release
   examples
   Package API <api/diffpy.pdffit2>

=======
Indices
=======

* :ref:`genindex`
* :ref:`search`
