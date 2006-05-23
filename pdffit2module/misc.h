// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Michael A.G. Aivazis
//                        California Institute of Technology
//                        (C) 1998-2005  All Rights Reserved
// 
//  <LicenseText>
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#if !defined(pypdffit2_misc_h)
#define pypdffit2_misc_h

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

// isel
extern char pypdffit2_isel__doc__[];
extern char pypdffit2_isel__name__[];
extern "C"
PyObject * pypdffit2_isel(PyObject *, PyObject *);

// idesel
extern char pypdffit2_idesel__doc__[];
extern char pypdffit2_idesel__name__[];
extern "C"
PyObject * pypdffit2_idesel(PyObject *, PyObject *);

// jsel
extern char pypdffit2_jsel__doc__[];
extern char pypdffit2_jsel__name__[];
extern "C"
PyObject * pypdffit2_jsel(PyObject *, PyObject *);

// jdesel
extern char pypdffit2_jdesel__doc__[];
extern char pypdffit2_jdesel__name__[];
extern "C"
PyObject * pypdffit2_jdesel(PyObject *, PyObject *);

// bang (bond_angle in c)
extern char pypdffit2_bang__doc__[];
extern char pypdffit2_bang__name__[];
extern "C"
PyObject * pypdffit2_bang(PyObject *, PyObject *);

// blen (bond_length in c)
extern char pypdffit2_blen__doc__[];
extern char pypdffit2_blen__name__[];
extern "C"
PyObject * pypdffit2_blen(PyObject *, PyObject *);

// show_scat
extern char pypdffit2_show_scat__doc__[];
extern char pypdffit2_show_scat__name__[];
extern "C"
PyObject * pypdffit2_show_scat(PyObject *, PyObject *);

// set_scat
extern char pypdffit2_set_scat__doc__[];
extern char pypdffit2_set_scat__name__[];
extern "C"
PyObject * pypdffit2_set_scat(PyObject *, PyObject *);

// set_scat_c
extern char pypdffit2_set_scat_c__doc__[];
extern char pypdffit2_set_scat_c__name__[];
extern "C"
PyObject * pypdffit2_set_scat_c(PyObject *, PyObject *);

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

// srat
extern char pypdffit2_srat__doc__[];
extern char pypdffit2_srat__name__[];
extern "C"
PyObject * pypdffit2_srat(PyObject *, PyObject *);

// delta
extern char pypdffit2_delta__doc__[];
extern char pypdffit2_delta__name__[];
extern "C"
PyObject * pypdffit2_delta(PyObject *, PyObject *);

// gamma
extern char pypdffit2_gamma__doc__[];
extern char pypdffit2_gamma__name__[];
extern "C"
PyObject * pypdffit2_gamma(PyObject *, PyObject *);

// dscale
extern char pypdffit2_dscale__doc__[];
extern char pypdffit2_dscale__name__[];
extern "C"
PyObject * pypdffit2_dscale(PyObject *, PyObject *);

// qsig
extern char pypdffit2_qsig__doc__[];
extern char pypdffit2_qsig__name__[];
extern "C"
PyObject * pypdffit2_qsig(PyObject *, PyObject *);

// qalp
extern char pypdffit2_qalp__doc__[];
extern char pypdffit2_qalp__name__[];
extern "C"
PyObject * pypdffit2_qalp(PyObject *, PyObject *);

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

#endif

// version
// $Id: misc.h,v 1.9 2005/11/04 19:12:06 farrowch Exp $

// End of file
