#!/usr/bin/env python
##############################################################################
#
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

"""Definition of __version__."""

#  We do not use the other three variables, but can be added back if needed.
__all__ = ["__date__", "__git_commit__", "__timestamp__", "__version__"]

# obtain version information
from importlib.metadata import version

__version__ = version("diffpy.pdffit2")
__date__ = __all__[0]
__git_commit__ = __all__[1]
__timestamp__ = __all__[2]

# End of file