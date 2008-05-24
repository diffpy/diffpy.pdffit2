########################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""PDFfit2 - real space structure refinement program.
Classes:
    PdfFit
Routines:
    redirect_stdout
"""

__id__ = "$Id$"

from PdfFit import PdfFit
from output import redirect_stdout
from pdffit2 import is_element

# obtain version information
from pkg_resources import get_distribution
__version__ = get_distribution(__name__).version
__date__ = __id__.split()[3]

# cleanup what should not get imported
del get_distribution

#  End of file
