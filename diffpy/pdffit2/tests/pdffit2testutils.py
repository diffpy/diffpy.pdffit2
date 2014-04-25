#!/usr/bin/env python
##############################################################################
#
# diffpy.srreal     by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Helper routines for running other unit tests.
Import of this module suppresses the chatty output from the C++ extension.
"""


import os.path
import diffpy.pdffit2

# silence the C++ engine output
diffpy.pdffit2.redirect_stdout(open(os.path.devnull, 'w'))

# path variables
thisfile = locals().get('__file__', 'file.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

def datafile(filename):
    """prepend testdata_dir to filename.
    """
    return os.path.join(testdata_dir, filename)

# End of file
