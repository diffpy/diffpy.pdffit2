#!/usr/bin/env python
##############################################################################
#
# diffpy.pdffit2    by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2012 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Convenience module for executing unit tests for all pdffit2 dependencies

python -m diffpy.pdffit2.tests.rundeps
"""


if __name__ == '__main__':
    import sys
    from diffpy.pdffit2.tests import testdeps
    # produce zero exit code for a successful test
    sys.exit(not testdeps().wasSuccessful())

# End of file
