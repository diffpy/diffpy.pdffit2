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

"""Unit tests for the diffpy.pdffit2 package.
"""


def testsuite():
    '''Build a unit tests suite for the diffpy.pdffit2 package.

    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        diffpy.pdffit2.tests.ExceptionsTest
        diffpy.pdffit2.tests.TestPdfFit
        diffpy.pdffit2.tests.TestPhaseFractions
        diffpy.pdffit2.tests.TestShapeFactors
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    mobj = None
    for mname in modulenames:
        exec ('import %s as mobj' % mname)
        suite.addTests(loader.loadTestsFromModule(mobj))
    return suite


def test():
    '''Execute all unit tests for the diffpy.pdffit2 package.
    Return a unittest TestResult object.
    '''
    import unittest
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


def testdeps():
    '''Execute all unit tests for diffpy.pdffit2 and its dependencies.

    Return a unittest TestResult object.
    '''
    import unittest
    modulenames = '''
        diffpy.pdffit2.tests
        diffpy.Structure.tests
    '''.split()
    suite = unittest.TestSuite()
    t = None
    for mname in modulenames:
        exec ('from %s import testsuite as t' % mname)
        suite.addTests(t())
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file
