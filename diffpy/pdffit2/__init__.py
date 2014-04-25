#!/usr/bin/env python
##############################################################################
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
##############################################################################

"""PDFfit2 - real space structure refinement program.
Classes:
    PdfFit
Routines:
    redirect_stdout
"""


from diffpy.pdffit2.version import __version__, __date__
from diffpy.pdffit2.pdffit import PdfFit
from diffpy.pdffit2.output import redirect_stdout
from diffpy.pdffit2.pdffit2 import is_element


# unit tests
def test():
    '''Execute all unit tests for the diffpy.pdffit2 package.
    Return a unittest TestResult object.
    '''
    from diffpy.pdffit2.tests import test
    return test()

# End of file
