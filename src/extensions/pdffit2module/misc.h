/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Bindings from python to c++ PdfFit class.
*
* Comments:
*
***********************************************************************/

#ifndef PYPDFFIT2_MISC_H
#define PYPDFFIT2_MISC_H

// copyright
extern char pypdffit2_copyright__name__[];
extern char pypdffit2_copyright__doc__[];
extern "C"
PyObject * pypdffit2_copyright(PyObject *, PyObject *);

// create
extern char pypdffit2_create__name__[];
extern char pypdffit2_create__doc__[];
extern "C"
PyObject * pypdffit2_create(PyObject *, PyObject *);

// read_struct
extern char pypdffit2_read_struct__name__[];
extern char pypdffit2_read_struct__doc__[];
extern "C"
PyObject * pypdffit2_read_struct(PyObject *, PyObject *);

// read_struct_string
extern char pypdffit2_read_struct_string__name__[];
extern char pypdffit2_read_struct_string__doc__[];
extern "C"
PyObject * pypdffit2_read_struct_string(PyObject *, PyObject *);

// read_data
extern char pypdffit2_read_data__name__[];
extern char pypdffit2_read_data__doc__[];
extern "C"
PyObject * pypdffit2_read_data(PyObject *, PyObject *);

// read_data_string
extern char pypdffit2_read_data_string__name__[];
extern char pypdffit2_read_data_string__doc__[];
extern "C"
PyObject * pypdffit2_read_data_string(PyObject *, PyObject *);

// read_data_arrays
extern char pypdffit2_read_data_arrays__name__[];
extern char pypdffit2_read_data_arrays__doc__[];
extern "C"
PyObject * pypdffit2_read_data_arrays(PyObject *, PyObject *);

// pdfrange
extern char pypdffit2_pdfrange__doc__[];
extern char pypdffit2_pdfrange__name__[];
extern "C"
PyObject * pypdffit2_pdfrange(PyObject *, PyObject *);

// reset
extern char pypdffit2_reset__doc__[];
extern char pypdffit2_reset__name__[];
extern "C"
PyObject * pypdffit2_reset(PyObject *, PyObject *);

// alloc
extern char pypdffit2_alloc__doc__[];
extern char pypdffit2_alloc__name__[];
extern "C"
PyObject * pypdffit2_alloc(PyObject *, PyObject *);

// calc
extern char pypdffit2_calc__doc__[];
extern char pypdffit2_calc__name__[];
extern "C"
PyObject * pypdffit2_calc(PyObject *, PyObject *);

// refine
extern char pypdffit2_refine__doc__[];
extern char pypdffit2_refine__name__[];
extern "C"
PyObject * pypdffit2_refine(PyObject *, PyObject *);

// refine_step
extern char pypdffit2_refine_step__doc__[];
extern char pypdffit2_refine_step__name__[];
extern "C"
PyObject * pypdffit2_refine_step(PyObject *, PyObject *);

// save_pdf
extern char pypdffit2_save_pdf__doc__[];
extern char pypdffit2_save_pdf__name__[];
extern "C"
PyObject * pypdffit2_save_pdf(PyObject *, PyObject *);

// save_dif
extern char pypdffit2_save_dif__doc__[];
extern char pypdffit2_save_dif__name__[];
extern "C"
PyObject * pypdffit2_save_dif(PyObject *, PyObject *);

// save_res
extern char pypdffit2_save_res__doc__[];
extern char pypdffit2_save_res__name__[];
extern "C"
PyObject * pypdffit2_save_res(PyObject *, PyObject *);

// save_struct
extern char pypdffit2_save_struct__doc__[];
extern char pypdffit2_save_struct__name__[];
extern "C"
PyObject * pypdffit2_save_struct(PyObject *, PyObject *);

// show_struct
extern char pypdffit2_show_struct__doc__[];
extern char pypdffit2_show_struct__name__[];
extern "C"
PyObject * pypdffit2_show_struct(PyObject *, PyObject *);

// constrain to string
extern char pypdffit2_constrain_str__doc__[];
extern char pypdffit2_constrain_str__name__[];
extern "C"
PyObject * pypdffit2_constrain_str(PyObject *, PyObject *);

// constrain_int to integer
extern char pypdffit2_constrain_int__doc__[];
extern char pypdffit2_constrain_int__name__[];
extern "C"
PyObject * pypdffit2_constrain_int(PyObject *, PyObject *);

// setpar_dbl
extern char pypdffit2_setpar_dbl__doc__[];
extern char pypdffit2_setpar_dbl__name__[];
extern "C"
PyObject * pypdffit2_setpar_dbl(PyObject *, PyObject *);

// setpar with RefVar
extern char pypdffit2_setpar_RV__doc__[];
extern char pypdffit2_setpar_RV__name__[];
extern "C"
PyObject * pypdffit2_setpar_RV(PyObject *, PyObject *);

// setvar
extern char pypdffit2_setvar__doc__[];
extern char pypdffit2_setvar__name__[];
extern "C"
PyObject * pypdffit2_setvar(PyObject *, PyObject *);

// getvar
extern char pypdffit2_getvar__doc__[];
extern char pypdffit2_getvar__name__[];
extern "C"
PyObject * pypdffit2_getvar(PyObject *, PyObject *);

// getR
extern char pypdffit2_getR__doc__[];
extern char pypdffit2_getR__name__[];
extern "C"
PyObject * pypdffit2_getR(PyObject *, PyObject *);

// getpdf_fit
extern char pypdffit2_getpdf_fit__doc__[];
extern char pypdffit2_getpdf_fit__name__[];
extern "C"
PyObject * pypdffit2_getpdf_fit(PyObject *, PyObject *);

// getpdf_obs
extern char pypdffit2_getpdf_obs__doc__[];
extern char pypdffit2_getpdf_obs__name__[];
extern "C"
PyObject * pypdffit2_getpdf_obs(PyObject *, PyObject *);

// getpdf_diff
extern char pypdffit2_getpdf_diff__doc__[];
extern char pypdffit2_getpdf_diff__name__[];
extern "C"
PyObject * pypdffit2_getpdf_diff(PyObject *, PyObject *);

// getcrw
extern char pypdffit2_getcrw__doc__[];
extern char pypdffit2_getcrw__name__[];
extern "C"
PyObject * pypdffit2_getcrw(PyObject *, PyObject *);

// getrw
extern char pypdffit2_getrw__doc__[];
extern char pypdffit2_getrw__name__[];
extern "C"
PyObject * pypdffit2_getrw(PyObject *, PyObject *);

// getpar
extern char pypdffit2_getpar__doc__[];
extern char pypdffit2_getpar__name__[];
extern "C"
PyObject * pypdffit2_getpar(PyObject *, PyObject *);

// fixpar
extern char pypdffit2_fixpar__doc__[];
extern char pypdffit2_fixpar__name__[];
extern "C"
PyObject * pypdffit2_fixpar(PyObject *, PyObject *);

// freepar
extern char pypdffit2_freepar__doc__[];
extern char pypdffit2_freepar__name__[];
extern "C"
PyObject * pypdffit2_freepar(PyObject *, PyObject *);

// setphase
extern char pypdffit2_setphase__doc__[];
extern char pypdffit2_setphase__name__[];
extern "C"
PyObject * pypdffit2_setphase(PyObject *, PyObject *);

// setdata
extern char pypdffit2_setdata__doc__[];
extern char pypdffit2_setdata__name__[];
extern "C"
PyObject * pypdffit2_setdata(PyObject *, PyObject *);

// psel
extern char pypdffit2_psel__doc__[];
extern char pypdffit2_psel__name__[];
extern "C"
PyObject * pypdffit2_psel(PyObject *, PyObject *);

// pdesel
extern char pypdffit2_pdesel__doc__[];
extern char pypdffit2_pdesel__name__[];
extern "C"
PyObject * pypdffit2_pdesel(PyObject *, PyObject *);

// selectAtomType
extern char pypdffit2_selectAtomType__doc__[];
extern char pypdffit2_selectAtomType__name__[];
extern "C"
PyObject * pypdffit2_selectAtomType(PyObject *, PyObject *);

// selectAtomIndex
extern char pypdffit2_selectAtomIndex__doc__[];
extern char pypdffit2_selectAtomIndex__name__[];
extern "C"
PyObject * pypdffit2_selectAtomIndex(PyObject *, PyObject *);

// selectAll
extern char pypdffit2_selectAll__doc__[];
extern char pypdffit2_selectAll__name__[];
extern "C"
PyObject * pypdffit2_selectAll(PyObject *, PyObject *);

// selectNone
extern char pypdffit2_selectNone__doc__[];
extern char pypdffit2_selectNone__name__[];
extern "C"
PyObject * pypdffit2_selectNone(PyObject *, PyObject *);

// bond_angle
extern char pypdffit2_bond_angle__doc__[];
extern char pypdffit2_bond_angle__name__[];
extern "C"
PyObject * pypdffit2_bond_angle(PyObject *, PyObject *);

// bond_length_atoms
extern char pypdffit2_bond_length_atoms__doc__[];
extern char pypdffit2_bond_length_atoms__name__[];
extern "C"
PyObject * pypdffit2_bond_length_atoms(PyObject *, PyObject *);

// bond_length_types
extern char pypdffit2_bond_length_types__doc__[];
extern char pypdffit2_bond_length_types__name__[];
extern "C"
PyObject * pypdffit2_bond_length_types(PyObject *, PyObject *);

// get_scat_string
extern char pypdffit2_get_scat_string__doc__[];
extern char pypdffit2_get_scat_string__name__[];
extern "C"
PyObject * pypdffit2_get_scat_string(PyObject *, PyObject *);

// get_scat
extern char pypdffit2_get_scat__doc__[];
extern char pypdffit2_get_scat__name__[];
extern "C"
PyObject * pypdffit2_get_scat(PyObject *, PyObject *);

// set_scat
extern char pypdffit2_set_scat__doc__[];
extern char pypdffit2_set_scat__name__[];
extern "C"
PyObject * pypdffit2_set_scat(PyObject *, PyObject *);

// reset_scat
extern char pypdffit2_reset_scat__doc__[];
extern char pypdffit2_reset_scat__name__[];
extern "C"
PyObject * pypdffit2_reset_scat(PyObject *, PyObject *);

// lat
extern char pypdffit2_lat__doc__[];
extern char pypdffit2_lat__name__[];
extern "C"
PyObject * pypdffit2_lat(PyObject *, PyObject *);

// x
extern char pypdffit2_x__doc__[];
extern char pypdffit2_x__name__[];
extern "C"
PyObject * pypdffit2_x(PyObject *, PyObject *);

// y
extern char pypdffit2_y__doc__[];
extern char pypdffit2_y__name__[];
extern "C"
PyObject * pypdffit2_y(PyObject *, PyObject *);

// z
extern char pypdffit2_z__doc__[];
extern char pypdffit2_z__name__[];
extern "C"
PyObject * pypdffit2_z(PyObject *, PyObject *);

// u11
extern char pypdffit2_u11__doc__[];
extern char pypdffit2_u11__name__[];
extern "C"
PyObject * pypdffit2_u11(PyObject *, PyObject *);

// u22
extern char pypdffit2_u22__doc__[];
extern char pypdffit2_u22__name__[];
extern "C"
PyObject * pypdffit2_u22(PyObject *, PyObject *);

// u33
extern char pypdffit2_u33__doc__[];
extern char pypdffit2_u33__name__[];
extern "C"
PyObject * pypdffit2_u33(PyObject *, PyObject *);

// u12
extern char pypdffit2_u12__doc__[];
extern char pypdffit2_u12__name__[];
extern "C"
PyObject * pypdffit2_u12(PyObject *, PyObject *);

// u13
extern char pypdffit2_u13__doc__[];
extern char pypdffit2_u13__name__[];
extern "C"
PyObject * pypdffit2_u13(PyObject *, PyObject *);

// u23
extern char pypdffit2_u23__doc__[];
extern char pypdffit2_u23__name__[];
extern "C"
PyObject * pypdffit2_u23(PyObject *, PyObject *);

// occ
extern char pypdffit2_occ__doc__[];
extern char pypdffit2_occ__name__[];
extern "C"
PyObject * pypdffit2_occ(PyObject *, PyObject *);

// pscale
extern char pypdffit2_pscale__doc__[];
extern char pypdffit2_pscale__name__[];
extern "C"
PyObject * pypdffit2_pscale(PyObject *, PyObject *);

// spdiameter
extern char pypdffit2_spdiameter__doc__[];
extern char pypdffit2_spdiameter__name__[];
extern "C"
PyObject * pypdffit2_spdiameter(PyObject *, PyObject *);

// stepcut
extern char pypdffit2_stepcut__doc__[];
extern char pypdffit2_stepcut__name__[];
extern "C"
PyObject * pypdffit2_stepcut(PyObject *, PyObject *);

// sratio
extern char pypdffit2_sratio__doc__[];
extern char pypdffit2_sratio__name__[];
extern "C"
PyObject * pypdffit2_sratio(PyObject *, PyObject *);

// delta2
extern char pypdffit2_delta2__doc__[];
extern char pypdffit2_delta2__name__[];
extern "C"
PyObject * pypdffit2_delta2(PyObject *, PyObject *);

// delta1
extern char pypdffit2_delta1__doc__[];
extern char pypdffit2_delta1__name__[];
extern "C"
PyObject * pypdffit2_delta1(PyObject *, PyObject *);

// dscale
extern char pypdffit2_dscale__doc__[];
extern char pypdffit2_dscale__name__[];
extern "C"
PyObject * pypdffit2_dscale(PyObject *, PyObject *);

// qdamp
extern char pypdffit2_qdamp__doc__[];
extern char pypdffit2_qdamp__name__[];
extern "C"
PyObject * pypdffit2_qdamp(PyObject *, PyObject *);

// qbroad
extern char pypdffit2_qbroad__doc__[];
extern char pypdffit2_qbroad__name__[];
extern "C"
PyObject * pypdffit2_qbroad(PyObject *, PyObject *);

// rcut
extern char pypdffit2_rcut__doc__[];
extern char pypdffit2_rcut__name__[];
extern "C"
PyObject * pypdffit2_rcut(PyObject *, PyObject *);

// get_atoms
extern char pypdffit2_get_atoms__doc__[];
extern char pypdffit2_get_atoms__name__[];
extern "C"
PyObject * pypdffit2_get_atoms(PyObject *, PyObject *);

// num_atoms
extern char pypdffit2_num_atoms__doc__[];
extern char pypdffit2_num_atoms__name__[];
extern "C"
PyObject * pypdffit2_num_atoms(PyObject *, PyObject *);

// get_atom_types
extern char pypdffit2_get_atom_types__doc__[];
extern char pypdffit2_get_atom_types__name__[];
extern "C"
PyObject * pypdffit2_get_atom_types(PyObject *, PyObject *);

// num_phases
extern char pypdffit2_num_phases__doc__[];
extern char pypdffit2_num_phases__name__[];
extern "C"
PyObject * pypdffit2_num_phases(PyObject *, PyObject *);

// num_datasets
extern char pypdffit2_num_datasets__doc__[];
extern char pypdffit2_num_datasets__name__[];
extern "C"
PyObject * pypdffit2_num_datasets(PyObject *, PyObject *);

// phase_fractions
extern char pypdffit2_phase_fractions__doc__[];
extern char pypdffit2_phase_fractions__name__[];
extern "C"
PyObject * pypdffit2_phase_fractions(PyObject *, PyObject *);

// redirect_stdout
extern char pypdffit2_redirect_stdout__doc__[];
extern char pypdffit2_redirect_stdout__name__[];
extern "C"
PyObject * pypdffit2_redirect_stdout(PyObject *, PyObject *);

// restore_stdout
extern char pypdffit2_restore_stdout__doc__[];
extern char pypdffit2_restore_stdout__name__[];
extern "C"
PyObject * pypdffit2_restore_stdout(PyObject *, PyObject *);

// is_element
extern char pypdffit2_is_element__doc__[];
extern char pypdffit2_is_element__name__[];
extern "C"
PyObject * pypdffit2_is_element(PyObject *, PyObject *);

#endif	// PYPDFFIT2_MISC_H
