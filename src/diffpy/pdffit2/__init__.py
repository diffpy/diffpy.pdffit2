#!/usr/bin/env python
##############################################################################
#
# (c) 2006 trustees of the Michigan State University.
# All rights reserved.
# (c) 2025 The Trustees of Columbia University in the City of New York.
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

from diffpy.pdffit2.output import redirect_stdout  # noqa

# package version
from diffpy.pdffit2.version import __date__, __version__  # noqa

# silence the pyflakes syntax checker
assert __version__ or True

# End of file
