#!/usr/bin/env python
##############################################################################
#
# (c) 2012 Trustees of the Columbia University in the City of New York.
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


# Consider upgrading to pytest
# import sys
#
# import pytest
#
# if __name__ == "__main__":
#     # show output results from every test function
#     args = ["-v"]
#     # show the message output for skipped and expected failure tests
#     if len(sys.argv) > 1:
#         args.extend(sys.argv[1:])
#     print("pytest arguments: {}".format(args))
#     # call pytest and exit with the return code from pytest
#     exit_res = pytest.main(args)
#     sys.exit(exit_res)

# End of file
