#!/usr/bin/env python

"""Unit tests for PdfFit.py."""

import unittest

import pytest

from diffpy.pdffit2 import PdfFit, pdffit2
from diffpy.structure import loadStructure

# ----------------------------------------------------------------------------


class TestPdfFit(unittest.TestCase):
    places = 6

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile, capture_output):
        self.datafile = datafile
        self.capture_output = capture_output

    def setUp(self):
        self.P = PdfFit()
        return

    def tearDown(self):
        del self.P
        return

    def test__exportAll(self):
        "check PdfFit._exportAll()"
        ns = {}
        self.P._exportAll(ns)
        self.assertEqual("ALL", ns["ALL"])
        self.assertEqual("FSQR", ns["FSQR"])
        self.assertEqual("N", ns["N"])
        self.assertIs("N", ns["N"])
        self.assertIs(self.P.lat, ns["lat"])
        self.assertEqual(self.P.reset, ns["reset"])
        return

    #   def test_intro(self):
    #       """check PdfFit.intro()
    #       """
    #       return

    def test_add_structure(self):
        """Check PdfFit.add_structure()"""
        ni = loadStructure(self.datafile("Ni.stru"))
        self.P.add_structure(ni)
        self.assertEqual(4, self.P.num_atoms())
        return

    #   def test_read_struct(self):
    #       """check PdfFit.read_struct()
    #       """
    #       return
    #
    #   def test_read_struct_string(self):
    #       """check PdfFit.read_struct_string()
    #       """
    #       return
    #
    #   def test_read_data(self):
    #       """check PdfFit.read_data()
    #       """
    #       return

    def test_read_data_string(self):
        """Check PdfFit.read_data_string()"""
        pf = self.P
        with open(self.datafile("300K.gr")) as fp:
            s = fp.read()
        self.assertEqual([], pf.data_files)
        pf.read_data_string(s, "N", 32, 0.03, "lmo")
        self.assertEqual(1, len(pf.data_files))
        gobs = pf.getpdf_obs()
        self.assertEqual(2000, len(gobs))
        self.assertEqual(0.384, gobs[-1])
        self.assertEqual(0.03, pf.getvar("qdamp"))
        return

    #   def test_read_data_lists(self):
    #       """check PdfFit.read_data_lists()
    #       """
    #       return
    #
    #   def test_pdfrange(self):
    #       """check PdfFit.pdfrange()
    #       """
    #       return
    #
    #   def test_reset(self):
    #       """check PdfFit.reset()
    #       """
    #       return

    def test_alloc(self):
        """Check PdfFit.alloc()"""
        # alloc and read_struct can be called in any order.
        self.P.alloc("X", 25, 0.0, 0.01, 10, 1000)
        # without a structure calculated PDF is all zero
        self.P.calc()
        Gzero = self.P.getpdf_fit()
        self.assertEqual(1000 * [0.0], Gzero)
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.calc()
        # check r-values
        r = self.P.getR()
        self.assertEqual(1000, len(r))
        for i in range(1000):
            self.assertAlmostEqual(0.01 * (i + 1), r[i], self.places)
        Gfit_alloc_read = self.P.getpdf_fit()
        # now try the other order
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 25, 0.0, 0.01, 10, 1000)
        self.P.calc()
        Gfit_read_alloc = self.P.getpdf_fit()
        # and they should be the same
        self.assertEqual(Gfit_read_alloc, Gfit_alloc_read)
        return

    #   def test_calc(self):
    #       """check PdfFit.calc()
    #       """
    #       return
    #
    #   def test_refine(self):
    #       """check PdfFit.refine()
    #       """
    #       return
    #
    #   def test_refine_step(self):
    #       """check PdfFit.refine_step()
    #       """
    #       return
    #
    #   def test_save_pdf(self):
    #       """check PdfFit.save_pdf()
    #       """
    #       return
    #
    #   def test_save_pdf_string(self):
    #       """check PdfFit.save_pdf_string()
    #       """
    #       return
    #
    #   def test_save_dif(self):
    #       """check PdfFit.save_dif()
    #       """
    #       return
    #
    #   def test_save_dif_string(self):
    #       """check PdfFit.save_dif_string()
    #       """
    #       return
    #
    #   def test_save_res(self):
    #       """check PdfFit.save_res()
    #       """
    #       return
    #
    #   def test_save_res_string(self):
    #       """check PdfFit.save_res_string()
    #       """
    #       return

    def test_get_structure(self):
        """Check PdfFit.get_structure()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        stru1 = self.P.get_structure(1)
        self.assertEqual(4, len(stru1))
        self.assertEqual("Ni", stru1[0].element)
        stru2 = self.P.get_structure(2)
        self.assertEqual(56, len(stru2))
        self.assertEqual("Ti", stru2[-1].element)
        return

    #   def test_save_struct(self):
    #       """check PdfFit.save_struct()
    #       """
    #       return
    #
    #   def test_save_struct_string(self):
    #       """check PdfFit.save_struct_string()
    #       """
    #       return
    #
    #   def test_show_struct(self):
    #       """check PdfFit.show_struct()
    #       """
    #       return
    #
    #   def test_constrain(self):
    #       """check PdfFit.constrain()
    #       """
    #       return

    def test_setpar(self):
        """Check PdfFit.setpar()"""
        pf = self.P
        pf.read_struct(self.datafile("Ni.stru"))
        pf.setpar(1, "lat(1)")
        self.assertEqual(3.52, pf.getpar(1))
        pf.setpar(1, 4.0)
        self.assertEqual(4, pf.getpar(1))
        pf.setpar(1, pf.lat("a"))
        self.assertEqual(3.52, pf.getpar(1))
        return

    def test_setvar(self):
        """Check PdfFit.setvar()"""
        pf = self.P
        pf.read_struct(self.datafile("Ni.stru"))
        pf.setvar(pf.delta1, 1.2)
        self.assertEqual(1.2, pf.getvar(pf.delta1))
        pf.setvar("delta1", 1.7)
        self.assertEqual(1.7, pf.getvar("delta1"))
        return

    #   def test_getvar(self):
    #       """check PdfFit.getvar()
    #       """
    #       return
    #
    #   def test_getrw(self):
    #       """check PdfFit.getrw()
    #       """
    #       return
    #
    #   def test_getR(self):
    #       """check PdfFit.getR()
    #       """
    #       return
    #
    #   def test_getpdf_fit(self):
    #       """check PdfFit.getpdf_fit()
    #       """
    #       return
    #
    #   def test_getpdf_obs(self):
    #       """check PdfFit.getpdf_obs()
    #       """
    #       return
    #
    #   def test_getpdf_diff(self):
    #       """check PdfFit.getpdf_diff()
    #       """
    #       return

    def test_get_atoms(self):
        """Check PdfFit.get_atoms()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.setphase(1)
        a1 = self.P.get_atoms()
        a2 = self.P.get_atoms(2)
        self.assertEqual(4 * ["NI"], a1)
        self.assertEqual(8 * ["PB"] + 24 * ["O"] + 8 * ["SC"] + 8 * ["W"] + 8 * ["TI"], a2)
        return

    def test_get_atom_types(self):
        """Check PdfFit.get_atom_types()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.setphase(1)
        atp1 = self.P.get_atom_types()
        atp2 = self.P.get_atom_types(2)
        self.assertEqual(["NI"], atp1)
        self.assertEqual(["PB", "O", "SC", "W", "TI"], atp2)
        return

    def test_num_phases(self):
        """Check PdfFit.num_phases()"""
        self.assertEqual(0, self.P.num_phases())
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertEqual(1, self.P.num_phases())
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.assertEqual(2, self.P.num_phases())
        self.P.reset()
        self.assertEqual(0, self.P.num_phases())
        return

    def test_num_datasets(self):
        """Check PdfFit.num_datasets()"""
        self.assertEqual(0, self.P.num_datasets())
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.5)
        self.assertEqual(1, self.P.num_datasets())
        # failed data should not increase num_datasets
        try:
            self.P.read_data(self.datafile("badNi.dat"))
        except (RuntimeError, TypeError, NameError, ValueError, IOError):
            pass
        self.assertEqual(1, self.P.num_datasets())
        # alloc should increase number of datasets
        # alloc requires a loaded structure
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 30.0, 0.05, 2, 10, 100)
        self.assertEqual(2, self.P.num_datasets())
        self.P.reset()
        self.assertEqual(0, self.P.num_datasets())
        return

    def test_getcrw(self):
        """Check PdfFit.getcrw()"""
        import numpy

        self.assertEqual(0, self.P.num_datasets())
        # Setting qmax=0 so that partial crw are not disturbed by
        # termination ripples.
        self.P.read_data(self.datafile("Ni.dat"), "X", 0.0, 0.0)
        # crw is empty before data refinement
        self.assertEqual([], self.P.getcrw())
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.pdfrange(1, 2, 19)
        self.P.refine()
        crw19 = numpy.array(self.P.getcrw())
        self.assertTrue(numpy.all(crw19 >= 0.0))
        # check that crw19 is non decreasing
        self.assertTrue(numpy.all(numpy.diff(crw19) >= 0.0))
        # check that crw19 and getrw give the same value
        rw19 = crw19[-1]
        self.assertAlmostEqual(self.P.getrw(), rw19, self.places)
        # renormalize cumulative Rw and compare with Rw at r=15
        Gobs19 = numpy.array(self.P.getpdf_obs())
        Gnorm19 = numpy.sqrt(numpy.sum(Gobs19**2))
        r = numpy.array(self.P.getR())
        idx = numpy.nonzero(r <= 15)[0]
        Gnorm15 = numpy.sqrt(numpy.sum(Gobs19[idx] ** 2))
        i15 = idx[-1]
        rw15 = crw19[i15] * Gnorm19 / Gnorm15
        self.P.pdfrange(1, 2, r[i15] + 1e-5)
        self.P.refine()
        self.assertAlmostEqual(self.P.getrw(), rw15, self.places)
        return

    def test_getcrw_two_datasets(self):
        """Check that getcrw() and getrw() are consistent for two datasets."""
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.pdfrange(1, 2, 8)
        self.P.read_data(self.datafile("300K.gr"), "N", 32.0, 0.0)
        self.P.pdfrange(2, 1, 11)
        self.P.read_struct(self.datafile("Ni.stru"))
        # mess lattice parameters to have comparable Rw contributions
        self.P.setvar("lat(1)", 3)
        self.P.setvar("lat(2)", 3)
        self.P.setvar("lat(3)", 3)
        self.P.refine()
        rwtot = self.P.getrw()
        self.assertTrue(rwtot > 0.0)
        self.P.setdata(1)
        rw1 = self.P.getcrw()[-1]
        self.P.setdata(2)
        rw2 = self.P.getcrw()[-1]
        self.assertAlmostEqual(rwtot**2, rw1**2 + rw2**2, self.places)
        return

    #   def test_getpar(self):
    #       """check PdfFit.getpar()
    #       """
    #       return

    def test_fixpar(self):
        """Check PdfFit.fixpar()"""
        self.P.fixpar("all")
        self.assertRaises(TypeError, self.P.fixpar, "x")
        return

    def test_freepar(self):
        """Check PdfFit.freepar()"""
        self.P.freepar("all")
        self.assertRaises(TypeError, self.P.freepar, "x")
        return

    #   def test_setphase(self):
    #       """check PdfFit.setphase()
    #       """
    #       return
    #
    #   def test_setdata(self):
    #       """check PdfFit.setdata()
    #       """
    #       return
    #
    def test_psel(self):
        """Check PdfFit.psel()"""

        def doalloc():
            self.P.alloc("X", 30.0, 0.05, 2, 10, 100)
            return

        self.assertRaises(pdffit2.unassignedError, self.P.psel, 0)
        self.assertRaises(pdffit2.unassignedError, self.P.psel, 1)
        self.P.read_struct(self.datafile("Ni.stru"))
        doalloc()
        self.P.calc()
        G1 = self.P.getpdf_fit()
        self.P.reset()
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        doalloc()
        self.P.calc()
        G2 = self.P.getpdf_fit()
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        doalloc()
        self.P.pdesel("ALL")
        self.P.psel(1)
        self.P.calc()
        self.assertEqual(G1, self.P.getpdf_fit())
        self.P.pdesel("ALL")
        self.P.psel(2)
        self.P.calc()
        self.assertEqual(G2, self.P.getpdf_fit())
        self.P.psel("ALL")
        self.P.calc()
        Gall = self.P.getpdf_fit()
        dGmax = max([abs(g1 + g2 - gall) for g1, g2, gall in zip(G1, G2, Gall)])
        self.assertAlmostEqual(0, dGmax, self.places)
        self.assertRaises(pdffit2.unassignedError, self.P.psel, 10)
        self.assertRaises(pdffit2.unassignedError, self.P.psel, 0)
        self.assertRaises(pdffit2.unassignedError, self.P.psel, -100)
        return

    def test_pdesel(self):
        """Check PdfFit.pdesel()"""

        def doalloc():
            self.P.alloc("X", 30.0, 0.05, 2, 10, 100)
            return

        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, 0)
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, 1)
        self.P.read_struct(self.datafile("Ni.stru"))
        doalloc()
        self.P.calc()
        G1 = self.P.getpdf_fit()
        self.P.reset()
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        doalloc()
        self.P.calc()
        G2 = self.P.getpdf_fit()
        self.P.reset()
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        doalloc()
        self.P.psel("ALL")
        self.P.pdesel(2)
        self.P.calc()
        self.assertEqual(G1, self.P.getpdf_fit())
        self.P.psel("ALL")
        self.P.pdesel(1)
        self.P.calc()
        self.assertEqual(G2, self.P.getpdf_fit())
        self.P.pdesel("ALL")
        self.P.calc()
        G0 = self.P.getpdf_fit()
        self.assertEqual([0.0] * len(G0), G0)
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, 10)
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, 0)
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, -100)
        return

    #
    #   def test_selectAtomType(self):
    #       """check PdfFit.selectAtomType()
    #       """
    #       return
    #
    #   def test_selectAtomIndex(self):
    #       """check PdfFit.selectAtomIndex()
    #       """
    #       return
    #
    #   def test_selectAll(self):
    #       """check PdfFit.selectAll()
    #       """
    #       return
    #
    #   def test_selectNone(self):
    #       """check PdfFit.selectNone()
    #       """
    #       return

    def test_bond_angle(self):
        """Check PdfFit.bond_angle()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        a, e = self.P.bond_angle(1, 2, 3)
        self.assertAlmostEqual(60.0, a, self.places)
        self.assertRaises(ValueError, self.P.bond_angle, 0, 1, 2)
        self.assertRaises(ValueError, self.P.bond_angle, 1, 2, 7)
        return

    def test_bang(self):
        "check PdfFit.bang() function"
        self.P.read_struct(self.datafile("Ni.stru"))
        out = self.capture_output(self.P.bang, 1, 2, 3).strip()
        self.assertTrue(out.endswith("60 degrees"))
        self.assertTrue(out.startswith("NI (#1) - NI (#2) - NI (#3)"))
        return

    def test_bond_length_atoms(self):
        """Check PdfFit.bond_length_atoms()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        dij, ddij = self.P.bond_length_atoms(1, 5)
        self.assertAlmostEqual(4.03635, dij, self.places)
        self.P.setphase(1)
        self.assertRaises(ValueError, self.P.bond_length_atoms, 1, 5)
        return

    def test_bond_length_types(self):
        """Check PdfFit.bond_length_types()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        dPbO = self.P.bond_length_types("Pb", "O", 0.1, 3.0)
        # check if keys are present
        self.assertTrue("dij" in dPbO)
        self.assertTrue("ddij" in dPbO)
        self.assertTrue("ij0" in dPbO)
        self.assertTrue("ij1" in dPbO)
        # check if they have the same length
        npts = len(dPbO["dij"])
        self.assertEqual(npts, len(dPbO["ddij"]))
        self.assertEqual(npts, len(dPbO["ij0"]))
        self.assertEqual(npts, len(dPbO["ij1"]))
        # 8 Pb atoms have coordination 12 in perovskite structure
        self.assertEqual(8 * 12, len(dPbO["dij"]))
        self.P.setphase(1)
        dfcc = self.P.bond_length_types("ALL", "ALL", 0.1, 2.6)
        # 4 Ni atoms with coordination 12
        self.assertEqual(4 * 12, len(dfcc["dij"]))
        # invalid element
        self.assertRaises(ValueError, self.P.bond_length_types, "Ni", "Nix", 0.1, 5.0)
        # check indices ij0
        allij0 = sum(dfcc["ij0"], tuple())
        self.assertEqual(0, min(allij0))
        self.assertEqual(3, max(allij0))
        # check indices ij1
        allij1 = sum(dfcc["ij1"], tuple())
        self.assertEqual(1, min(allij1))
        self.assertEqual(4, max(allij1))
        # check index values
        ij0check = [(i1 - 1, j1 - 1) for i1, j1 in dfcc["ij1"]]
        self.assertEqual(ij0check, dfcc["ij0"])
        # test valid element which is not present in the structure
        dnone = self.P.bond_length_types("Ni", "Au", 0.1, 5.0)
        self.assertEqual(0, len(dnone["dij"]))
        self.assertEqual(0, len(dnone["ddij"]))
        self.assertEqual(0, len(dnone["ij0"]))
        self.assertEqual(0, len(dnone["ij1"]))
        return

    def test_blen(self):
        """Check PdfFit.blen()"""
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        blen = self.P.blen
        o = self.capture_output(blen, 1, 5).strip()
        self.assertTrue(o.endswith("4.03635 A"))
        self.assertTrue("PB (#1)" in o)
        self.assertTrue("PB (#5)" in o)
        self.assertRaises(ValueError, blen, 1, 99)
        self.assertRaises(ValueError, blen, 0, 1)
        o1 = self.capture_output(blen, 1, 1, 0.1, 1)
        self.assertTrue("No pairs found" in o1)
        o2 = self.capture_output(blen, 1, 50, 0.1, 1)
        self.assertEqual("", o2)
        o3 = self.capture_output(blen, "Sc", "O", 0.5, 2.3).strip()
        self.assertEqual(1 + 48, len(o3.split("\n")))
        self.assertEqual(6, o3.count("SC (#33)"))
        self.assertEqual(2, o3.count("O (#9)"))
        self.assertRaises(TypeError, blen, "Sc", "O", 0.5)
        return

    #   def test_show_scat(self):
    #       """check PdfFit.show_scat()
    #       """
    #       return
    #
    #   def test_get_scat_string(self):
    #       """check PdfFit.get_scat_string()
    #       """
    #       return

    def test_get_scat(self):
        """Check PdfFit.get_scat()"""
        # x-ray scattering factors
        fPb = self.P.get_scat("X", "Pb")
        self.assertEqual(82.0, fPb)
        fTi = self.P.get_scat("X", "tI")
        self.assertEqual(22.0, fTi)
        # neutron scattering lengths
        bPb = self.P.get_scat("N", "PB")
        self.assertAlmostEqual(9.401, bPb, 3)
        bTi = self.P.get_scat("N", "ti")
        self.assertAlmostEqual(-3.370, bTi, 3)
        # exceptions
        self.assertRaises(ValueError, self.P.get_scat, "N", "zz")
        self.assertRaises(ValueError, self.P.get_scat, "Z", "Ti")
        return

    def test_set_scat(self):
        """Check PdfFit.set_scat()"""
        # raises exception when no phase exists
        self.assertRaises(pdffit2.unassignedError, self.P.set_scat, "N", "Ti", -11)
        # check if it is local to phase
        fPb = self.P.get_scat("X", "Pb")
        bPb = self.P.get_scat("N", "Pb")
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.set_scat("X", "Pb", 142)
        self.assertEqual(142, self.P.get_scat("X", "Pb"))
        self.assertEqual(bPb, self.P.get_scat("N", "Pb"))
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.assertEqual(fPb, self.P.get_scat("X", "Pb"))
        self.P.setphase(1)
        self.assertEqual(142, self.P.get_scat("X", "Pb"))
        self.P.setphase(2)
        self.assertEqual(fPb, self.P.get_scat("X", "Pb"))
        # check exception for invalid inputs
        self.assertRaises(ValueError, self.P.set_scat, "Z", "C", 123)
        self.assertRaises(ValueError, self.P.set_scat, "X", "ZZ", 123)
        return

    def test_reset_scat(self):
        """Check PdfFit.reset_scat()"""
        # raises exception when no phase exists
        self.assertRaises(pdffit2.unassignedError, self.P.reset_scat, "Ti")
        # check if it is local to phase
        fPb = self.P.get_scat("X", "Pb")
        bPb = self.P.get_scat("N", "Pb")
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.set_scat("X", "Pb", 142)
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.P.set_scat("N", "Pb", -17)
        self.P.setphase(1)
        self.assertNotEqual(fPb, self.P.get_scat("X", "Pb"))
        self.P.reset_scat("Pb")
        self.assertEqual(fPb, self.P.get_scat("X", "Pb"))
        self.P.setphase(2)
        self.assertNotEqual(bPb, self.P.get_scat("N", "Pb"))
        self.P.reset_scat("Pb")
        self.assertEqual(bPb, self.P.get_scat("N", "Pb"))
        # check exception for invalid inputs
        self.assertRaises(ValueError, self.P.reset_scat, "Zz")
        return

    def test_num_atoms(self):
        """Check PdfFit.num_atoms()"""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertEqual(4, self.P.num_atoms())
        self.P.read_struct(self.datafile("PbScW25TiO3.stru"))
        self.assertEqual(56, self.P.num_atoms())
        self.P.setphase(1)
        self.assertEqual(4, self.P.num_atoms())
        self.P.setphase(2)
        self.assertEqual(56, self.P.num_atoms())
        return

    def test_lat(self):
        """Check PdfFit.lat()"""
        pf = self.P
        pf.read_struct(self.datafile("Ni.stru"))
        for i in ("a", "b", "c", 1, 2, 3):
            self.assertEqual(3.52, pf.getvar(pf.lat(i)))
        for i in ("alpha", "beta", "gamma", 4, 5, 6):
            self.assertEqual(90, pf.getvar(pf.lat(i)))
        return

    def test_xyz(self):
        """Check PdfFit.x() PdfFit.y(), PdfFit.z()"""
        pf = self.P
        pf.read_struct(self.datafile("Ni.stru"))
        self.assertEqual(0.5, pf.getvar(pf.x(3)))
        self.assertEqual(0, pf.getvar(pf.y(3)))
        self.assertEqual(0.5, pf.getvar(pf.z(3)))
        return

    def test_uij(self):
        """Check PdfFit.uij()"""
        ni = loadStructure(self.datafile("Ni.stru"))
        ni[2].anisotropy = True
        ni[2].U11, ni[2].U22, ni[2].U33 = 1, 2, 3
        ni[2].U12, ni[2].U13, ni[2].U23 = 4, 5, 6
        pf = self.P
        pf.add_structure(ni)
        self.assertEqual(1, pf.getvar(pf.u11(3)))
        self.assertEqual(2, pf.getvar(pf.u22(3)))
        self.assertEqual(3, pf.getvar(pf.u33(3)))
        self.assertEqual(4, pf.getvar(pf.u12(3)))
        self.assertEqual(5, pf.getvar(pf.u13(3)))
        self.assertEqual(6, pf.getvar(pf.u23(3)))
        return

    def test_occ(self):
        """Check PdfFit.occ()"""
        pf = self.P
        pf.read_struct(self.datafile("Ni.stru"))
        for i in range(1, 5):
            self.assertEqual(1, pf.getvar(pf.occ(i)))
        return

    #   def test_pscale(self):
    #       """check PdfFit.pscale()
    #       """
    #       return
    #
    #   def test_pscale(self):
    #       """check PdfFit.pscale()
    #       """
    #       return
    #
    #   def test_sratio(self):
    #       """check PdfFit.sratio()
    #       """
    #       return
    #
    #   def test_delta1(self):
    #       """check PdfFit.delta1()
    #       """
    #       return
    #
    #   def test_delta2(self):
    #       """check PdfFit.delta2()
    #       """
    #       return
    #
    #   def test_dscale(self):
    #       """check PdfFit.dscale()
    #       """
    #       return
    #
    #   def test_qdamp(self):
    #       """check PdfFit.qdamp()
    #       """
    #       return
    #
    #   def test_qbroad(self):
    #       """check PdfFit.qbroad()
    #       """
    #       return
    #
    #   def test_rcut(self):
    #       """check PdfFit.rcut()
    #       """
    #       return
    #

    def test___init__(self):
        """Check PdfFit.__init__()"""
        output_true = self.capture_output(PdfFit, create_intro=True).strip()
        output_false = self.capture_output(PdfFit, create_intro=False).strip()

        self.assertGreater(len(output_true), 0)
        self.assertEqual(len(output_false), 0)

        return


#
#   def test__PdfFit__getRef(self):
#       """check PdfFit._PdfFit__getRef()
#       """
#       return

# End of class TestPdfFit

if __name__ == "__main__":
    unittest.main()

# End of file
