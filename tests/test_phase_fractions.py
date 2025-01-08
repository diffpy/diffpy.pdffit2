#!/usr/bin/env python

"""Unit tests for phase fraction calculations."""

import unittest

import pytest

from diffpy.pdffit2 import PdfFit


##############################################################################
class TestPhaseFractions(unittest.TestCase):
    places = 4

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.alloc("X", 0.0, 0.05, 0.1, 10, 200)
        self.P.alloc("N", 0.0, 0.05, 0.1, 10, 200)
        return

    def tearDown(self):
        del self.P
        return

    def test_xray_fractions(self):
        """test_xray_fractions -- check phase fractions in x-ray dataset."""
        self.P.setdata(1)
        ph = self.P.phase_fractions()
        bb1 = 28**2
        bb2 = ((8 * 82 + 24 * 8 + 4 * 21 + 2 * 74 + 2 * 22) / 40.0) ** 2
        self.assertAlmostEqual(1.0, sum(ph["atom"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph["cell"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph["mass"]), self.places)
        self.assertAlmostEqual(bb2 / bb1, ph["atom"][0] / ph["atom"][1], self.places)
        self.assertAlmostEqual(bb2 / bb1 * 40.0 / 4.0, ph["cell"][0] / ph["cell"][1], self.places)
        mavg1 = 58.69
        mavg2 = (8 * 207.19 + 24 * 15.994 + 4 * 44.956 + 2 * 183.85 + 2 * 47.90) / 40.0
        self.assertAlmostEqual(bb2 / bb1 * mavg1 / mavg2, ph["mass"][0] / ph["mass"][1], self.places)
        self.assertEqual(0.0, sum(ph["stdatom"]))
        self.assertEqual(0.0, sum(ph["stdcell"]))
        self.assertEqual(0.0, sum(ph["stdmass"]))
        self.P.setphase(1)
        self.P.setvar("pscale", 2.0)
        ph2 = self.P.phase_fractions()
        self.assertAlmostEqual(1.0, sum(ph2["atom"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph2["cell"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph2["mass"]), self.places)
        self.assertAlmostEqual(2.0, ph2["atom"][0] / ph2["atom"][1] / (ph["atom"][0] / ph["atom"][1]), self.places)
        self.assertAlmostEqual(2.0, ph2["cell"][0] / ph2["cell"][1] / (ph["cell"][0] / ph["cell"][1]), self.places)
        self.assertAlmostEqual(2.0, ph2["mass"][0] / ph2["mass"][1] / (ph["mass"][0] / ph["mass"][1]), self.places)
        return

    def test_neutron_fractions(self):
        """test_neutron_fractions -- check phase fractions in neutron dataset."""
        self.P.setdata(2)
        ph = self.P.phase_fractions()
        bb1 = 10.31**2
        bPb = 9.4012
        bO = 5.8054
        bSc = 12.11
        bW = 4.75518
        bTi = -3.37013
        bb2 = ((8 * bPb + 24 * bO + 4 * bSc + 2 * bW + 2 * bTi) / 40.0) ** 2
        self.assertAlmostEqual(1.0, sum(ph["atom"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph["cell"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph["mass"]), self.places)
        self.assertAlmostEqual(bb2 / bb1, ph["atom"][0] / ph["atom"][1], self.places)
        self.assertAlmostEqual(bb2 / bb1 * 40.0 / 4.0, ph["cell"][0] / ph["cell"][1], self.places)
        mavg1 = 58.69
        mavg2 = (8 * 207.19 + 24 * 15.994 + 4 * 44.956 + 2 * 183.85 + 2 * 47.90) / 40.0
        self.assertAlmostEqual(bb2 / bb1 * mavg1 / mavg2, ph["mass"][0] / ph["mass"][1], self.places)
        self.assertEqual(0.0, sum(ph["stdatom"]))
        self.assertEqual(0.0, sum(ph["stdcell"]))
        self.assertEqual(0.0, sum(ph["stdmass"]))
        self.P.setphase(1)
        self.P.setvar("pscale", 2.0)
        ph2 = self.P.phase_fractions()
        self.assertAlmostEqual(1.0, sum(ph2["atom"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph2["cell"]), self.places)
        self.assertAlmostEqual(1.0, sum(ph2["mass"]), self.places)
        self.assertAlmostEqual(2.0, ph2["atom"][0] / ph2["atom"][1] / (ph["atom"][0] / ph["atom"][1]), self.places)
        self.assertAlmostEqual(2.0, ph2["cell"][0] / ph2["cell"][1] / (ph["cell"][0] / ph["cell"][1]), self.places)
        self.assertAlmostEqual(2.0, ph2["mass"][0] / ph2["mass"][1] / (ph["mass"][0] / ph["mass"][1]), self.places)
        return


# End of class TestSphereEnvelope

if __name__ == "__main__":
    unittest.main()

# End of file
