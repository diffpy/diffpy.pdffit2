#!/usr/bin/env python

"""Unit tests for particle shape envelope factors.
"""

# version
__id__ = '$Id$'

import os
import unittest
import numpy

# useful variables
thisfile = locals().get('__file__', 'file.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

from diffpy.pdffit2 import PdfFit
from diffpy.pdffit2 import pdffit2

def testdata(filename):
    """prepend testdata_dir to filename.
    """
    return os.path.join(testdata_dir, filename)

##############################################################################
class TestSphereEnvelope(unittest.TestCase):

    places = 6

    def setUp(self):
        self.P = PdfFit()
        return

    def tearDown(self):
        del self.P
        return

    def test_calculation(self):
        """check calculation of sphere envelope factor
        """
        self.P.read_struct(testdata('Ni.stru'))
        self.P.alloc('X', 0.0, 0.05, 0.1, 10, 200)
        self.P.calc()
        d = 8.0
        r = numpy.array(self.P.getR())
        G0 = numpy.array(self.P.getpdf_fit())
        self.P.setvar('spdiameter', d)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        fsph = 1.0 - 1.5*r/d + 0.5*(r/d)**3
        fsph[r > d] = 0.0
        dG = (G0*fsph - G1)
        msd = numpy.dot(dG, dG)/len(r)
        self.assertAlmostEqual(0.0, numpy.sqrt(msd), self.places)
        return

    def test_refinement(self):
        """check refinement of sphere envelope factor
        """
        from numpy.random import rand
        dcheck = 8.0
        dstart = 12.0
        self.P.read_struct(testdata('Ni.stru'))
        self.P.alloc('X', 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar('spdiameter', dcheck)
        self.P.calc()
        r = numpy.array(self.P.getR())
        Gd8 = numpy.array(self.P.getpdf_fit())
        Gd8noise = Gd8
        Gd8noise[::2] += 0.01
        Gd8noise[1::2] -= 0.01
        self.P.reset()
        self.P.read_struct(testdata('Ni.stru'))
        self.P.read_data_lists('X', 0.0, 0.05, list(r), list(Gd8noise))
        self.P.constrain('spdiameter', '@8')
        self.P.setpar(8, dstart)
        self.P.refine()
        dfinal = self.P.getvar('spdiameter')
        self.assertAlmostEqual(dcheck, dfinal, 3)
        return

# End of class TestSphereEnvelope

if __name__ == '__main__':
    unittest.main()

# End of file
