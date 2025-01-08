#!/usr/bin/env python
##############################################################################
#
# (c) 2006 trustees of the Michigan State University.
# All rights reserved.
# (c) 2024 The Trustees of Columbia University in the City of New York.
# All rights reserved.
#
# File coded by: Billinge Group members and community contributors.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.pdffit2/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""PDFfit2 - real space structure refinement program."""

# WARNING: Do NOT remove the isort: off/on comments in this file.
# These tags are used to prevent isort from reordering imports in this file.
# __version__ must be initialized before importing C++ extensions.

# isort: off
# Import the package version before C++ extensions are loaded.
from diffpy.pdffit2.output import redirect_stdout
from diffpy.pdffit2.version import __date__, __version__

# Import C++ related modules since the __version__ attribute is used.
from diffpy.pdffit2.pdffit import PdfFit
from diffpy.pdffit2.pdffit2 import is_element

# isort: on

# Ensure all necessary components are imported and initialized
assert __version__ or True
assert __date__ or True
assert all((PdfFit, redirect_stdout, is_element))

# End of file
