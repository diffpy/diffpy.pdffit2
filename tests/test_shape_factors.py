#!/usr/bin/env python

"""Unit tests for particle shape envelope factors."""


import unittest

import numpy
import pytest

from diffpy.pdffit2 import PdfFit, pdffit2


def spherefactor(r, d):
    """Calculate spherical envelope correction.

    r -- PDF radius
    d -- diameter of spherical particle

    Return numpy array of shape correction envelope.
    """
    r1 = numpy.array(r)
    fsph = 1.0 - 1.5 * r1 / d + 0.5 * (r1 / d) ** 3
    fsph[r1 > d] = 0.0
    return fsph


##############################################################################
class TestSphereEnvelope(unittest.TestCase):

    places = 6

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        return

    def tearDown(self):
        self.P = None
        return

    def test_calculation(self):
        """Check calculation of sphere envelope factor."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.calc()
        d = 8.0
        r = numpy.array(self.P.getR())
        G0 = numpy.array(self.P.getpdf_fit())
        self.P.setvar("spdiameter", d)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        dG = G0 * spherefactor(r, d) - G1
        msd = numpy.dot(dG, dG) / len(r)
        self.assertAlmostEqual(0.0, numpy.sqrt(msd), self.places)
        return

    def test_refinement(self):
        """Check refinement of sphere envelope factor."""
        dcheck = 8.0
        dstart = 12.0
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("spdiameter", dcheck)
        self.P.calc()
        r = numpy.array(self.P.getR())
        Gd8 = numpy.array(self.P.getpdf_fit())
        Gd8noise = Gd8
        Gd8noise[::2] += 0.01
        Gd8noise[1::2] -= 0.01
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data_lists("X", 0.0, 0.05, list(r), list(Gd8noise))
        self.P.constrain("spdiameter", "@8")
        self.P.setpar(8, dstart)
        self.P.refine()
        dfinal = self.P.getvar("spdiameter")
        self.assertAlmostEqual(dcheck, dfinal, 3)
        return

    def test_twophase_calculation(self):
        """Check PDF calculation for 2 phases with different spdiameters."""
        d1 = 6
        d2 = 9
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("spdiameter", d1)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        self.P.reset()
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("spdiameter", d2)
        self.P.calc()
        G2 = numpy.array(self.P.getpdf_fit())
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setphase(1)
        self.P.setvar("spdiameter", d1)
        self.P.setphase(2)
        self.P.setvar("spdiameter", d2)
        self.P.calc()
        Gtot = numpy.array(self.P.getpdf_fit())
        dG = G1 + G2 - Gtot
        r = numpy.array(self.P.getR())
        msd = numpy.dot(dG, dG) / len(r)
        self.assertAlmostEqual(0.0, numpy.sqrt(msd), self.places)
        return

    def test_twophase_refinement(self):
        """Check PDF refinement of 2 phases that have different spdiameter."""
        dcheck1 = 8.0
        dstart1 = 8.2
        dcheck2 = 6.0
        dstart2 = 5.5
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("spdiameter", dcheck1)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        self.P.reset()
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("spdiameter", dcheck2)
        self.P.calc()
        G2 = numpy.array(self.P.getpdf_fit())
        r = numpy.array(self.P.getR())
        Gnoise = G1 + G2
        Gnoise[::2] += 0.01
        Gnoise[1::2] -= 0.01
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.read_data_lists("X", 0.0, 0.05, list(r), list(Gnoise))
        self.P.setphase(1)
        self.P.constrain("spdiameter", "@11")
        self.P.setphase(2)
        self.P.constrain("spdiameter", "@12")
        self.P.setpar(11, dstart1)
        self.P.setpar(12, dstart2)
        self.P.refine()
        dfinal2 = self.P.getvar("spdiameter")
        self.P.setphase(1)
        dfinal1 = self.P.getvar("spdiameter")
        self.assertAlmostEqual(dcheck1, dfinal1, 3)
        self.assertAlmostEqual(dcheck2, dfinal2, 3)
        return

    def test_spdiameter_io(self):
        """Check reading and writing of spdiameter from structure file."""
        import re

        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertEqual(0.0, self.P.getvar("spdiameter"))
        # engine should not write shape factor when not defined
        spdnone = self.P.save_struct_string(1)
        self.assertTrue(not re.search("(?m)^shape +sphere,", spdnone))
        self.P.setvar("spdiameter", 7)
        spd7 = self.P.save_struct_string(1)
        # spd7 should contain shape factor data
        self.assertTrue(re.search("(?m)^shape +sphere,", spd7))
        self.P.reset()
        self.P.read_struct_string(spd7)
        self.assertEqual(7.0, self.P.getvar("spdiameter"))
        # try to read without comma
        spd14 = re.sub("(?m)^shape +sphere.*$", "shape sphere 14.00", spd7)
        self.P.read_struct_string(spd14)
        self.assertEqual(14.0, self.P.getvar("spdiameter"))
        # try to read invalid shape data
        sinvalid = re.sub("(?m)^shape .*", "shape invalid, 1", spd7)
        self.assertRaises(pdffit2.structureError, self.P.read_struct_string, sinvalid)
        return


# End of class TestSphereEnvelope


##############################################################################
class TestStepCutEnvelope(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    places = 6

    def setUp(self):
        self.P = PdfFit()
        return

    def tearDown(self):
        self.P = None
        return

    def test_stepcut_calculation(self):
        """Check calculation of sphere envelope factor."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.calc()
        stepcut = 8.0
        r = numpy.array(self.P.getR())
        G0 = numpy.array(self.P.getpdf_fit())
        G0[r > stepcut] = 0.0
        self.P.setvar("stepcut", stepcut)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        dG = G0 - G1
        msd = numpy.dot(dG, dG) / len(r)
        self.assertAlmostEqual(0.0, numpy.sqrt(msd), self.places)
        return

    def test_twophase_stepcut_calculation(self):
        """Check PDF calculation for 2 phases with different spdiameters."""
        d1 = 6
        d2 = 9
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("stepcut", d1)
        self.P.calc()
        G1 = numpy.array(self.P.getpdf_fit())
        self.P.reset()
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setvar("stepcut", d2)
        self.P.calc()
        G2 = numpy.array(self.P.getpdf_fit())
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.setphase(1)
        self.P.setvar("stepcut", d1)
        self.P.setphase(2)
        self.P.setvar("stepcut", d2)
        self.P.calc()
        Gtot = numpy.array(self.P.getpdf_fit())
        dG = G1 + G2 - Gtot
        r = numpy.array(self.P.getR())
        msd = numpy.dot(dG, dG) / len(r)
        self.assertAlmostEqual(0.0, numpy.sqrt(msd), self.places)
        # G after step should be zero
        self.assertTrue(numpy.all(0 == Gtot[r > max(d1, d2)]))
        return

    def test_stepcut_io(self):
        """Check reading and writing of stepcut from structure file."""
        import re

        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertEqual(0.0, self.P.getvar("stepcut"))
        # engine should not write shape factor when not defined
        sscnone = self.P.save_struct_string(1)
        self.assertTrue(not re.search("(?m)^shape +stepcut,", sscnone))
        self.P.setvar("stepcut", 7)
        ssc7 = self.P.save_struct_string(1)
        # ssc7 should contain shape factor data
        self.assertTrue(re.search("(?m)^shape +stepcut,", ssc7))
        self.P.reset()
        self.P.read_struct_string(ssc7)
        self.assertEqual(7.0, self.P.getvar("stepcut"))
        # try to read without comma
        ssc14 = re.sub("(?m)^shape +stepcut.*$", "shape stepcut 14.00", ssc7)
        self.P.read_struct_string(ssc14)
        self.assertEqual(14.0, self.P.getvar("stepcut"))
        # try to read invalid shape data
        sinvalid = re.sub("(?m)^shape .*", "shape invalid, 1", ssc7)
        self.assertRaises(pdffit2.structureError, self.P.read_struct_string, sinvalid)
        return


# End of class TestStepCutEnvelope

if __name__ == "__main__":
    unittest.main()

# End of file
