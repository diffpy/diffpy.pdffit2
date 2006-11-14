#!/usr/bin/env python

"""Unit tests for PdfFit.py
"""

# version
__id__ = '$Id$'

import os
import unittest

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
class TestPdfFit(unittest.TestCase):

    places = 6

    def setUp(self):
        self.P = PdfFit()
        return

    def tearDown(self):
        del self.P
        return

#   def test_intro(self):
#       """check PdfFit.intro()
#       """
#       return
#
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
#
#   def test_read_data_string(self):
#       """check PdfFit.read_data_string()
#       """
#       return
#
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
#
#   def test_alloc(self):
#       """check PdfFit.alloc()
#       """
#       return
#
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
#
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
#
#   def test_setpar(self):
#       """check PdfFit.setpar()
#       """
#       return
#
#   def test_setvar(self):
#       """check PdfFit.setvar()
#       """
#       return
#
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

    def test_get_atoms(self):
        """check PdfFit.get_atoms()
        """
        self.P.read_struct(testdata('Ni.stru'))
        self.P.read_struct(testdata('PbScW25TiO3.stru'))
        self.P.setphase(1)
        a1 = self.P.get_atoms()
        a2 = self.P.get_atoms(2)
        self.assertEqual(4*['NI'], a1)
        self.assertEqual(8*['PB']+24*['O']+8*['SC']+8*['W']+8*['TI'], a2)
        return

    def test_get_atom_types(self):
        """check PdfFit.get_atom_types()
        """
        self.P.read_struct(testdata('Ni.stru'))
        self.P.read_struct(testdata('PbScW25TiO3.stru'))
        self.P.setphase(1)
        atp1 = self.P.get_atom_types()
        atp2 = self.P.get_atom_types(2)
        self.assertEqual(['NI'], atp1)
        self.assertEqual(['PB', 'O', 'SC', 'W', 'TI'], atp2)
        return

#   def test_getpar(self):
#       """check PdfFit.getpar()
#       """
#       return
#
#   def test_fixpar(self):
#       """check PdfFit.fixpar()
#       """
#       return
#
#   def test_freepar(self):
#       """check PdfFit.freepar()
#       """
#       return
#
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
#   def test_psel(self):
#       """check PdfFit.psel()
#       """
#       return
#
#   def test_pdesel(self):
#       """check PdfFit.pdesel()
#       """
#       return
#
#   def test_isel(self):
#       """check PdfFit.isel()
#       """
#       return
#
#   def test_idesel(self):
#       """check PdfFit.idesel()
#       """
#       return
#
#   def test_jsel(self):
#       """check PdfFit.jsel()
#       """
#       return
#
#   def test_jdesel(self):
#       """check PdfFit.jdesel()
#       """
#       return
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

    def test_bang(self):
        """check PdfFit.bang()
        """
        self.P.read_struct(testdata('Ni.stru'))
        a = self.P.bang(1, 2, 3)
        self.assertAlmostEqual(60.0, a, self.places)
        return

    def test_blen(self):
        """check PdfFit.blen()
        """
        self.P.read_struct(testdata('Ni.stru'))
        self.P.read_struct(testdata('PbScW25TiO3.stru'))
        dij = self.P.blen(1, 5)
        self.assertAlmostEqual(4.03635, dij, self.places)
        self.P.setphase(1)
        dijs = self.P.blen('Ni', 'ALL', 0.1,  5.0)
        self.assertEqual(216, len(dijs))
        # check if lengths are sorted
        dsteps = [ dijs[i][0] - dijs[i-1][0] for i in range(1, len(dijs)) ]
        minstep = min([0.0] + dsteps)
        self.failUnless(minstep >= 0)
        self.assertAlmostEqual(dijs[0][0], dijs[47][0], self.places)
        self.failIfAlmostEqual(dijs[0][0], dijs[48][0], self.places)
        dijs1 = self.P.blen(1, 'ALL', 0.1,  5.0)
        self.assertEqual(dijs, dijs1)
        self.assertEqual([], self.P.blen(77, 'ALL', 0.1, 5.0))
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
#
#   def test_set_scat(self):
#       """check PdfFit.set_scat()
#       """
#       return
#
#   def test_reset_scat(self):
#       """check PdfFit.reset_scat()
#       """
#       return
#
#   def test_num_atoms(self):
#       """check PdfFit.num_atoms()
#       """
#       return
#
#   def test_lat(self):
#       """check PdfFit.lat()
#       """
#       return
#
#   def test_x(self):
#       """check PdfFit.x()
#       """
#       return
#
#   def test_y(self):
#       """check PdfFit.y()
#       """
#       return
#
#   def test_z(self):
#       """check PdfFit.z()
#       """
#       return
#
#   def test_u11(self):
#       """check PdfFit.u11()
#       """
#       return
#
#   def test_u22(self):
#       """check PdfFit.u22()
#       """
#       return
#
#   def test_u33(self):
#       """check PdfFit.u33()
#       """
#       return
#
#   def test_u12(self):
#       """check PdfFit.u12()
#       """
#       return
#
#   def test_u13(self):
#       """check PdfFit.u13()
#       """
#       return
#
#   def test_u23(self):
#       """check PdfFit.u23()
#       """
#       return
#
#   def test_occ(self):
#       """check PdfFit.occ()
#       """
#       return
#
#   def test_pscale(self):
#       """check PdfFit.pscale()
#       """
#       return
#
#   def test_pfrac(self):
#       """check PdfFit.pfrac()
#       """
#       return
#
#   def test_srat(self):
#       """check PdfFit.srat()
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
#   def test_delta(self):
#       """check PdfFit.delta()
#       """
#       return
#
#   def test_gamma(self):
#       """check PdfFit.gamma()
#       """
#       return
#
#   def test_dscale(self):
#       """check PdfFit.dscale()
#       """
#       return
#
#   def test_qsig(self):
#       """check PdfFit.qsig()
#       """
#       return
#
#   def test_qalp(self):
#       """check PdfFit.qalp()
#       """
#       return
#
#   def test_rcut(self):
#       """check PdfFit.rcut()
#       """
#       return
#
#   def test___init__(self):
#       """check PdfFit.__init__()
#       """
#       return
#
#   def test__PdfFit__getRef(self):
#       """check PdfFit._PdfFit__getRef()
#       """
#       return

# End of class TestPdfFit

if __name__ == '__main__':
    unittest.main()

# End of file
