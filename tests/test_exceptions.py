#!/usr/bin/env python
##############################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################


import unittest

import pytest

from diffpy.pdffit2 import PdfFit, pdffit2


class read_structExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when structure file does not exist."""
        self.assertRaises(IOError, self.P.read_struct, "Nofile.stru")

    def test_structureError(self):
        """Raise pdffit2.structureError when structure is malformed."""
        self.assertRaises(pdffit2.structureError, self.P.read_struct, self.datafile("badNi.stru"))

    def test_structureErrorZeroVolume(self):
        """Raise pdffit2.structureError when unit cell volume is negative."""
        # I don't know how to test for this, but it's in the library code
        self.assertRaises(pdffit2.structureError, self.P.read_struct, self.datafile("badNiZeroVolume.stru"))


class read_dataExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when data file does not exist."""
        self.assertRaises(IOError, self.P.read_data, "Nofile.dat", "X", 25.0, 0.5)

    def test_dataError(self):
        """Raise pdffit2.dataError when data has improper spacing."""
        self.assertRaises(pdffit2.dataError, self.P.read_data, self.datafile("badNi.dat"), "X", 25.0, 0.5)


class read_data_listsExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.r_data = [0.1, 0.2]
        self.Gr_data = [1, 2, 3]
        self.qmax = 10
        self.qdamp = 0.5

    def tearDown(self):
        del self.P

    def test_ValueError1(self):
        """Raise ValueError when lists are of different length."""
        self.assertRaises(
            ValueError, self.P.read_data_lists, "X", self.qmax, self.qdamp, self.r_data, self.Gr_data
        )

    def test_ValueError2(self):
        """Raise ValueError when qmax < 0."""
        self.assertRaises(
            ValueError, self.P.read_data_lists, "X", -self.qmax, self.qdamp, self.r_data, self.Gr_data
        )

    def test_ValueError3(self):
        """Raise ValueError when qdamp < 0."""
        self.assertRaises(
            ValueError, self.P.read_data_lists, "X", self.qmax, -self.qdamp, self.r_data, self.Gr_data
        )

    def test_dataError(self):
        """Raise pdffit2.dataError when data has improper spacing."""
        r_data = [0.1, 0.52, 0.2]
        self.assertRaises(
            pdffit2.dataError, self.P.read_data_lists, "X", self.qmax, self.qdamp, r_data, self.Gr_data
        )


class pdfrangeExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.iset = 1
        self.rmin = 4.0
        self.rmax = 10.0

    def tearDown(self):
        del self.P

    def test_ValueError1(self):
        """Raise ValueError when iset does not exist."""
        self.assertRaises(ValueError, self.P.pdfrange, self.iset, self.rmin, self.rmax)

    def test_ValueError2(self):
        """Raise ValueError when rmax < rmin."""
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.5)
        self.assertRaises(ValueError, self.P.pdfrange, self.iset, self.rmax, self.rmin)

    def test_ValueError3(self):
        """Raise ValueError when range outside of data."""
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.5)
        self.assertRaises(ValueError, self.P.pdfrange, self.iset, -self.rmin, self.rmax)


class allocExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.qmax = 25
        self.qdamp = 0.5
        self.rmin = 4.0
        self.rmax = 10.0
        self.bin = 100

    def tearDown(self):
        del self.P

    def test_ValueError1(self):
        """Raise ValueError when qmax < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", -self.qmax, self.qdamp, self.rmin, self.rmax, self.bin)

    def test_ValueError2(self):
        """Raise ValueError when qdamp < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, -self.qdamp, self.rmin, self.rmax, self.bin)

    def test_ValueError3(self):
        """Raise ValueError when rmin < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, self.qdamp, -self.rmin, self.rmax, self.bin)

    def test_ValueError4(self):
        """Raise ValueError when rmax < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, self.qdamp, self.rmin, -self.rmax, self.bin)

    def test_ValueError5(self):
        """Raise ValueError when bin < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, self.qdamp, self.rmin, self.rmax, -self.bin)

    def test_ValueError6(self):
        """Raise ValueError when rmax < rmin."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, self.qdamp, self.rmax, self.rmin, self.bin)

    def test_ValueError7(self):
        """Raise ValueError when qdamp < 0."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.alloc, "X", self.qmax, self.qdamp, self.rmin, self.rmax, -self.bin)


class calcExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.P.read_struct(self.datafile("Ni.stru"))

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when no space has been allocated."""
        self.assertRaises(pdffit2.unassignedError, self.P.calc)


# PJ 2006-03-06
#
# test_calculationError raised exception, because for Qmax=0.5, rmax would
# increase to 4010A and this would throw exception when the size of bnd array
# would exceed MAXBND limit.  However, bnd vector can now grow, thus rmax
# is not limited and the following test would hang indefinitely.

#   def test_calculationError(self):
#       """raise pdffit2.calculationError when calculation cannot be done"""
#       self.P.alloc('X', 0.01, 0.5, 2, 10, 100)
#       self.assertRaises(pdffit2.calculationError, self.P.calc)


# class refineExceptions(unittest.TestCase):
# I'm not sure how to test these

#    def setUp(self):
#        self.P = PdfFit()
#
#    def test_calculationError(self):
#        """raise pdffit2.calculationError when model pdf can't be calculated"""
#        #self.assertRaises(pdffit2.calculationError, self.P.calc)
#
#    def test_constraintError(self):
#        """raise pdffit2.constraintError for bad constraint(s)"""
#        #self.assertRaises(pdffit2.constraintError, self.P.calc)


# class refine_stepExceptions(unittest.TestCase):
# I'm not sure how to test these

#    def setUp(self):
#        self.P = PdfFit()
#
#    def test_calculationError(self):
#        """raise pdffit2.calculationError when model pdf can't be calculated"""
#        #self.assertRaises(pdffit2.calculationError, self.P.calc)
#
#    def test_constraintError(self):
#        """raise pdffit2.constraintError for bad constraint(s)"""
#        #self.assertRaises(pdffit2.constraintError, self.P.calc)


class save_pdfExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.strufile = "temp.pdf"

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when structure cannot be saved."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 30.0, 0.05, 2, 10, 100)
        self.P.calc()
        self.assertRaises(IOError, self.P.save_pdf, 1, "nodir183160/" + self.strufile)

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when structure is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.save_pdf, 1, self.strufile)


class save_difExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.strufile = "temp.dif"

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when dif cannot be saved."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.alloc("X", 30.0, 0.05, 2, 10, 100)
        self.P.calc()
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.5)
        self.assertRaises(IOError, self.P.save_dif, 1, "nodir183160/" + self.strufile)

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when structure is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.save_dif, 1, self.strufile)


class save_resExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.resfile = "temp.res"

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when residual file cannot be saved."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 30.0, 0.0)
        self.P.constrain(self.P.lat(1), 1)
        self.P.setpar(1, 3.0)
        self.P.pdfrange(1, 2.0, 10.0)
        self.P.refine_step()
        self.assertRaises(IOError, self.P.save_res, "nodir183160/" + self.resfile)

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when structure is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.save_res, self.resfile)


class save_structExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    # Same code as show_struct

    def setUp(self):
        self.P = PdfFit()
        self.strufile = "temp.stru"

    def tearDown(self):
        del self.P

    def test_IOError(self):
        """Raise IOError when structure cannot be saved."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(IOError, self.P.save_struct, 1, "nodir183160/" + self.strufile)

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when structure is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.save_struct, 1, self.strufile)


class constrainExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.par = 1
        return

    def tearDown(self):
        self.P = None
        return

    def test_constraintError(self):
        """Raise constraintError when constraint is bad."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.constrain("x(1)", "junk+@1")
        self.P.setpar(1, 0.01)
        self.assertRaises(pdffit2.constraintError, self.P.calc)
        self.assertRaises(pdffit2.constraintError, self.P.refine)
        return

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when variable is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.constrain, self.P.x(1), self.par)
        return

    def test_ValueError(self):
        """Raise ValueError when a variable index does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.constrain, self.P.x(6), self.par)
        return

    def test_constrainNonRefVar(self):
        "raise constraintError when attempting to constrain non-refinables"
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(pdffit2.constraintError, self.P.constrain, "rcut", "@7")
        self.assertRaises(pdffit2.constraintError, self.P.constrain, "rcut", 13)
        self.assertRaises(pdffit2.constraintError, self.P.constrain, "stepcut", "@17")
        return


class setvarExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.val = 3.0

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when variable is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.setvar, self.P.lat(1), self.val)

    def test_ValueError(self):
        """Raise ValueError when a variable index does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.setvar, self.P.lat(7), self.val)


class getvarExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when variable is undefined."""
        self.assertRaises(pdffit2.unassignedError, self.P.getvar, self.P.pscale())

    def test_ValueError(self):
        """Raise ValueError when a variable index does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.getvar, self.P.lat(7))


class getRExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getR)


class getpdf_fitExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getpdf_fit)


class getpdf_obsExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getpdf_obs)


class getpdf_diffExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getpdf_diff)


class get_atomsExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.get_atoms)


class getparExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when parameter does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getpar, 1)

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when parameter does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.constrain(self.P.lat(1), 2)
        self.assertRaises(pdffit2.unassignedError, self.P.getpar, 1)


class pselExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()
        self.ip = 1

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, self.ip)

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, self.ip)


class pdeselExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.ip = 1

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, self.ip)

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(pdffit2.unassignedError, self.P.pdesel, self.ip)


class selectAtomTypeExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.iset = 1
        self.i = 1

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.selectAtomType, self.iset, "i", "Ni", True)

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        # selectAtomType should pass with one phase defined
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.selectAtomType(self.iset, "i", "Ni", True)
        self.P.selectAtomType(self.iset, "j", "Ni", False)
        # but fail for phase 2 which is not present
        self.assertRaises(pdffit2.unassignedError, self.P.selectAtomType, 2, "i", "Ca", True)

    def test_ijcharValueError(self):
        """Raise ValueError when ijchar is neither 'i' nor 'j'."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.selectAtomType(self.iset, "i", "Ni", True)
        self.P.selectAtomType(self.iset, "j", "Ni", True)
        self.assertRaises(ValueError, self.P.selectAtomType, self.iset, "k", "Ni", True)


class selectAtomIndexExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.iset = 1
        self.i = 1

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.selectAtomIndex, self.iset, "i", self.i, True)

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        # pass for phase 1
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.selectAtomIndex(self.iset, "i", 1, True)
        self.P.selectAtomIndex(self.iset, "i", 2, False)
        # fail for phase 2
        self.assertRaises(pdffit2.unassignedError, self.P.selectAtomIndex, 2, "i", 1, True)

    def test_ValueError(self):
        """Raise ValueError when selected atom does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(ValueError, self.P.selectAtomIndex, self.iset, "i", 6, True)


class selectAllExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.iset = 1
        self.i = 1

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.selectAll, self.iset, "i")

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        # fail when there is no dataset
        self.assertRaises(pdffit2.unassignedError, self.P.selectAll, self.iset, "i")
        # pass with dataset
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.selectAll(self.iset, "i")
        self.P.selectAll(self.iset, "j")
        # fail for phase 2
        self.assertRaises(pdffit2.unassignedError, self.P.selectAll, 2, "i")
        self.assertRaises(pdffit2.unassignedError, self.P.selectAll, 2, "j")


class selectNoneExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.iset = 1
        self.i = 1

    def tearDown(self):
        del self.P

    def test_unassignedError1(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.selectNone, self.iset, "i")

    def test_unassignedError2(self):
        """Raise pdffit2.unassignedError when set does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        # fail when there is no dataset
        self.assertRaises(pdffit2.unassignedError, self.P.selectNone, self.iset, "i")
        # pass with dataset
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.P.selectNone(self.iset, "i")
        self.P.selectNone(self.iset, "j")
        # fail for phase 2
        self.assertRaises(pdffit2.unassignedError, self.P.selectNone, 2, "i")
        self.assertRaises(pdffit2.unassignedError, self.P.selectNone, 2, "j")


class bangExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.a1 = 1
        self.a2 = 2
        self.a3 = 3

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.bang, self.a1, self.a2, self.a3)

    def test_ValueError1(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(ValueError, self.P.bang, 0, self.a2, self.a3)

    def test_ValueError2(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.bang, self.a1, -1, self.a3)

    def test_ValueError3(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.bang, self.a1, self.a2, 6)


class blenExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()
        self.a1 = 1
        self.a2 = 2

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when no data exists."""
        self.assertRaises(pdffit2.unassignedError, self.P.blen, self.a1, self.a2)

    def test_ValueError1(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.blen, 0, self.a2)

    def test_ValueError2(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.blen, self.a1, 6)

    def test_ValueError3(self):
        """Raise ValueError when selected atom(s) does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.assertRaises(ValueError, self.P.blen, 0, 6)


class show_scatExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.show_scat, "X")


# class set_scatExceptions(unittest.TestCase):
# I'm not sure how to use this function

#    def setUp(self):
#        self.P = PdfFit()
#
#    def test_unassignedError1(self):
#        """raise pdffit2.unassignedError when phase does not exist"""
#        #self.assertRaises(pdffit2.constraintError, self.P.calc)
#
#    def test_unassignedError2(self):
#        """raise pdffit2.unassignedError when phase does not exist"""
#        #self.assertRaises(pdffit2.constraintError, self.P.calc)
#
#    def test_ValueError(self):
#        """raise pdffit2.unassignedError when selected atom does not exist"""
#        #self.assertRaises(pdffit2.constraintError, self.P.calc)


class num_atomsExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when no atoms exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.num_atoms)


class fixparExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when parameter does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(pdffit2.unassignedError, self.P.fixpar, 1)


class freeparExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when parameter does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(pdffit2.unassignedError, self.P.freepar, 1)


class setphaseExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when phase does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(pdffit2.unassignedError, self.P.setphase, 2)


class setdataExceptions(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data set does not exist."""
        self.P.read_struct(self.datafile("Ni.stru"))
        self.P.read_data(self.datafile("Ni.dat"), "X", 25.0, 0.0)
        self.assertRaises(pdffit2.unassignedError, self.P.setdata, 2)


class getcrwExceptions(unittest.TestCase):
    def setUp(self):
        self.P = PdfFit()

    def tearDown(self):
        del self.P

    def test_unassignedError(self):
        """Raise pdffit2.unassignedError when data does not exist."""
        self.assertRaises(pdffit2.unassignedError, self.P.getcrw)


# main
if __name__ == "__main__":
    # suite = unittest.makeSuite(num_atomsExceptions)
    # unittest.TextTestRunner(verbosity=3).run(suite)
    # testcase = calcExceptions('test_unassignedError')
    # unittest.TextTestRunner(verbosity=3).run(testcase)
    unittest.main()

# End of file
