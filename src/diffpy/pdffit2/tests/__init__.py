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

"""Unit tests for the diffpy.pdffit2 package.
"""

import unittest


def testsuite(pattern=""):
    """Create a unit tests suite for the diffpy.pdffit2 package.

    Parameters
    ----------
    pattern : str, optional
        Regular expression pattern for selecting test cases.
        Select all tests when empty.  Ignore the pattern when
        any of unit test modules fails to import.

    Returns
    -------
    suite : `unittest.TestSuite`
        The TestSuite object containing the matching tests.
    """
    import re
    from itertools import chain
    from os.path import dirname

    from pkg_resources import resource_filename

    loader = unittest.defaultTestLoader
    thisdir = resource_filename(__name__, "")
    depth = __name__.count(".") + 1
    topdir = thisdir
    for i in range(depth):
        topdir = dirname(topdir)
    suite_all = loader.discover(thisdir, pattern="*Test*.py", top_level_dir=topdir)
    # always filter the suite by pattern to test-cover the selection code.
    suite = unittest.TestSuite()
    rx = re.compile(pattern)
    tsuites = list(chain.from_iterable(suite_all))
    tsok = all(isinstance(ts, unittest.TestSuite) for ts in tsuites)
    if not tsok:  # pragma: no cover
        return suite_all
    tcases = chain.from_iterable(tsuites)
    for tc in tcases:
        tcwords = tc.id().split(".")
        shortname = ".".join(tcwords[-3:])
        if rx.search(shortname):
            suite.addTest(tc)
    # verify all tests are found for an empty pattern.
    assert pattern or suite_all.countTestCases() == suite.countTestCases()
    return suite


def test():
    """Execute all unit tests for the diffpy.pdffit2 package.

    Returns
    -------
    result : `unittest.TestResult`
    """
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


def testdeps():
    """Execute all unit tests for diffpy.pdffit2 and its dependencies.

    Returns
    -------
    result : `unittest.TestResult`
    """
    from importlib import import_module

    modulenames = """
        diffpy.pdffit2.tests
        diffpy.structure.tests
    """.split()
    suite = unittest.TestSuite()
    for mname in modulenames:
        mod = import_module(mname)
        suite.addTests(mod.testsuite())
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file
