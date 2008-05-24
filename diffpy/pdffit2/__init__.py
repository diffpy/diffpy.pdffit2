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
from version import __version__, __date__

#  End of file
