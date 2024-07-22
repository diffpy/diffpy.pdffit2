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

"""Convenience module for executing all unit tests with

python -m diffpy.pdffit2.tests.run
"""


if __name__ == "__main__":
    import sys

    # show warnings by default
    if not sys.warnoptions:
        import os
        import warnings

        warnings.simplefilter("default")
        # also affect subprocesses
        os.environ["PYTHONWARNINGS"] = "default"
    from diffpy.pdffit2.tests import test

    # produce zero exit code for a successful test
    sys.exit(not test().wasSuccessful())

# End of file
