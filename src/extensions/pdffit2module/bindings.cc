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
* Method table for python module bindings.
*
* Comments:
*
***********************************************************************/
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "bindings.h"
#include "misc.h"

// the method table
struct PyMethodDef pypdffit2_methods[] = {

    //copyright
    {pypdffit2_copyright__name__, pypdffit2_copyright,
     METH_VARARGS, pypdffit2_copyright__doc__},

    //create
    {pypdffit2_create__name__, pypdffit2_create,
     METH_VARARGS, pypdffit2_create__doc__},

    //read_struct
    {pypdffit2_read_struct__name__, pypdffit2_read_struct,
     METH_VARARGS, pypdffit2_read_struct__doc__},

    //read_struct_string
    {pypdffit2_read_struct_string__name__, pypdffit2_read_struct_string,
     METH_VARARGS, pypdffit2_read_struct_string__doc__},

    //read_data
    {pypdffit2_read_data__name__, pypdffit2_read_data,
     METH_VARARGS, pypdffit2_read_data__doc__},

    //read_data_string
    {pypdffit2_read_data_string__name__, pypdffit2_read_data_string,
     METH_VARARGS, pypdffit2_read_data_string__doc__},

    //read_data_arrays
    {pypdffit2_read_data_arrays__name__, pypdffit2_read_data_arrays,
     METH_VARARGS, pypdffit2_read_data_arrays__doc__},

    //pdfrange
    {pypdffit2_pdfrange__name__, pypdffit2_pdfrange,
     METH_VARARGS, pypdffit2_pdfrange__doc__},

    //reset
    {pypdffit2_reset__name__, pypdffit2_reset,
     METH_VARARGS, pypdffit2_reset__doc__},

    //alloc
    {pypdffit2_alloc__name__, pypdffit2_alloc,
     METH_VARARGS, pypdffit2_alloc__doc__},

    //calc
    {pypdffit2_calc__name__, pypdffit2_calc,
     METH_VARARGS, pypdffit2_calc__doc__},

    //refine
    {pypdffit2_refine__name__, pypdffit2_refine,
     METH_VARARGS, pypdffit2_refine__doc__},

    //refine_step
    {pypdffit2_refine_step__name__, pypdffit2_refine_step,
     METH_VARARGS, pypdffit2_refine_step__doc__},

    //save_pdf
    {pypdffit2_save_pdf__name__, pypdffit2_save_pdf,
     METH_VARARGS, pypdffit2_save_pdf__doc__},

    //save_dif
    {pypdffit2_save_dif__name__, pypdffit2_save_dif,
     METH_VARARGS, pypdffit2_save_dif__doc__},

    //save_res
    {pypdffit2_save_res__name__, pypdffit2_save_res,
     METH_VARARGS, pypdffit2_save_res__doc__},

    //save_struct
    {pypdffit2_save_struct__name__, pypdffit2_save_struct,
     METH_VARARGS, pypdffit2_save_struct__doc__},

    //show_struct
    {pypdffit2_show_struct__name__, pypdffit2_show_struct,
     METH_VARARGS, pypdffit2_show_struct__doc__},

    //constrain_int
    {pypdffit2_constrain_int__name__, pypdffit2_constrain_int,
     METH_VARARGS, pypdffit2_constrain_int__doc__},

    //constrain_str
    {pypdffit2_constrain_str__name__, pypdffit2_constrain_str,
     METH_VARARGS, pypdffit2_constrain_str__doc__},

    //setpar_dbl
    {pypdffit2_setpar_dbl__name__, pypdffit2_setpar_dbl,
     METH_VARARGS, pypdffit2_setpar_dbl__doc__},

    //setpar_RV
    {pypdffit2_setpar_RV__name__, pypdffit2_setpar_RV,
     METH_VARARGS, pypdffit2_setpar_RV__doc__},

     //setvar
    {pypdffit2_setvar__name__, pypdffit2_setvar,
     METH_VARARGS, pypdffit2_setvar__doc__},

    //getvar
    {pypdffit2_getvar__name__, pypdffit2_getvar,
     METH_VARARGS, pypdffit2_getvar__doc__},

    //getR
    {pypdffit2_getR__name__, pypdffit2_getR,
     METH_VARARGS, pypdffit2_getR__doc__},

    //getpdf_fit
    {pypdffit2_getpdf_fit__name__, pypdffit2_getpdf_fit,
     METH_VARARGS, pypdffit2_getpdf_fit__doc__},

    //getpdf_obs
    {pypdffit2_getpdf_obs__name__, pypdffit2_getpdf_obs,
     METH_VARARGS, pypdffit2_getpdf_obs__doc__},

    //getpdf_diff
    {pypdffit2_getpdf_diff__name__, pypdffit2_getpdf_diff,
     METH_VARARGS, pypdffit2_getpdf_diff__doc__},

    //getcrw
    {pypdffit2_getcrw__name__, pypdffit2_getcrw,
     METH_VARARGS, pypdffit2_getcrw__doc__},

    //getrw
    {pypdffit2_getrw__name__, pypdffit2_getrw,
     METH_VARARGS, pypdffit2_getrw__doc__},

    //getpar
    {pypdffit2_getpar__name__, pypdffit2_getpar,
     METH_VARARGS, pypdffit2_getpar__doc__},

    //fixpar
    {pypdffit2_fixpar__name__, pypdffit2_fixpar,
     METH_VARARGS, pypdffit2_fixpar__doc__},

    //freepar
    {pypdffit2_freepar__name__, pypdffit2_freepar,
     METH_VARARGS, pypdffit2_freepar__doc__},

    //setphase
    {pypdffit2_setphase__name__, pypdffit2_setphase,
     METH_VARARGS, pypdffit2_setphase__doc__},

    //setdata
    {pypdffit2_setdata__name__, pypdffit2_setdata,
     METH_VARARGS, pypdffit2_setdata__doc__},

    //psel
    {pypdffit2_psel__name__, pypdffit2_psel,
     METH_VARARGS, pypdffit2_psel__doc__},

    //pdesel
    {pypdffit2_pdesel__name__, pypdffit2_pdesel,
     METH_VARARGS, pypdffit2_pdesel__doc__},

    //selectAtomType
    {pypdffit2_selectAtomType__name__, pypdffit2_selectAtomType,
     METH_VARARGS, pypdffit2_selectAtomType__doc__},

    //selectAtomIndex
    {pypdffit2_selectAtomIndex__name__, pypdffit2_selectAtomIndex,
     METH_VARARGS, pypdffit2_selectAtomIndex__doc__},

    //selectAll
    {pypdffit2_selectAll__name__, pypdffit2_selectAll,
     METH_VARARGS, pypdffit2_selectAll__doc__},

    //selectNone
    {pypdffit2_selectNone__name__, pypdffit2_selectNone,
     METH_VARARGS, pypdffit2_selectNone__doc__},

    //bond_angle
    {pypdffit2_bond_angle__name__, pypdffit2_bond_angle,
     METH_VARARGS, pypdffit2_bond_angle__doc__},

    //bond_length_atoms
    {pypdffit2_bond_length_atoms__name__, pypdffit2_bond_length_atoms,
     METH_VARARGS, pypdffit2_bond_length_atoms__doc__},

    //bond_length_types
    {pypdffit2_bond_length_types__name__, pypdffit2_bond_length_types,
     METH_VARARGS, pypdffit2_bond_length_types__doc__},

    //get_scat_string
    {pypdffit2_get_scat_string__name__, pypdffit2_get_scat_string,
     METH_VARARGS, pypdffit2_get_scat_string__doc__},

    //get_scat
    {pypdffit2_get_scat__name__, pypdffit2_get_scat,
     METH_VARARGS, pypdffit2_get_scat__doc__},

    //set_scat
    {pypdffit2_set_scat__name__, pypdffit2_set_scat,
     METH_VARARGS, pypdffit2_set_scat__doc__},

    //reset_scat
    {pypdffit2_reset_scat__name__, pypdffit2_reset_scat,
     METH_VARARGS, pypdffit2_reset_scat__doc__},

    //lat
    {pypdffit2_lat__name__, pypdffit2_lat,
     METH_VARARGS, pypdffit2_lat__doc__},

    //x
    {pypdffit2_x__name__, pypdffit2_x,
     METH_VARARGS, pypdffit2_x__doc__},

    //y
    {pypdffit2_y__name__, pypdffit2_y,
     METH_VARARGS, pypdffit2_y__doc__},

    //z
    {pypdffit2_z__name__, pypdffit2_z,
     METH_VARARGS, pypdffit2_z__doc__},

    //u11
    {pypdffit2_u11__name__, pypdffit2_u11,
     METH_VARARGS, pypdffit2_u11__doc__},

    //u22
    {pypdffit2_u22__name__, pypdffit2_u22,
     METH_VARARGS, pypdffit2_u22__doc__},

    //u33
    {pypdffit2_u33__name__, pypdffit2_u33,
     METH_VARARGS, pypdffit2_u33__doc__},

    //u12
    {pypdffit2_u12__name__, pypdffit2_u12,
     METH_VARARGS, pypdffit2_u12__doc__},

    //u13
    {pypdffit2_u13__name__, pypdffit2_u13,
     METH_VARARGS, pypdffit2_u13__doc__},

    //u23
    {pypdffit2_u23__name__, pypdffit2_u23,
     METH_VARARGS, pypdffit2_u23__doc__},

    //occ
    {pypdffit2_occ__name__, pypdffit2_occ,
     METH_VARARGS, pypdffit2_occ__doc__},

    //pscale
    {pypdffit2_pscale__name__, pypdffit2_pscale,
     METH_VARARGS, pypdffit2_pscale__doc__},

    //spdiameter
    {pypdffit2_spdiameter__name__, pypdffit2_spdiameter,
     METH_VARARGS, pypdffit2_spdiameter__doc__},

    //stepcut
    {pypdffit2_stepcut__name__, pypdffit2_stepcut,
     METH_VARARGS, pypdffit2_stepcut__doc__},

    //sratio
    {pypdffit2_sratio__name__, pypdffit2_sratio,
     METH_VARARGS, pypdffit2_sratio__doc__},

    //delta2
    {pypdffit2_delta2__name__, pypdffit2_delta2,
     METH_VARARGS, pypdffit2_delta2__doc__},

    //delta1
    {pypdffit2_delta1__name__, pypdffit2_delta1,
     METH_VARARGS, pypdffit2_delta1__doc__},

    //dscale
    {pypdffit2_dscale__name__, pypdffit2_dscale,
     METH_VARARGS, pypdffit2_dscale__doc__},

    //qdamp
    {pypdffit2_qdamp__name__, pypdffit2_qdamp,
     METH_VARARGS, pypdffit2_qdamp__doc__},

    //qbroad
    {pypdffit2_qbroad__name__, pypdffit2_qbroad,
     METH_VARARGS, pypdffit2_qbroad__doc__},

    //rcut
    {pypdffit2_rcut__name__, pypdffit2_rcut,
     METH_VARARGS, pypdffit2_rcut__doc__},

    //get_atoms
    {pypdffit2_get_atoms__name__, pypdffit2_get_atoms,
     METH_VARARGS, pypdffit2_get_atoms__doc__},

    //num_atoms
    {pypdffit2_num_atoms__name__, pypdffit2_num_atoms,
     METH_VARARGS, pypdffit2_num_atoms__doc__},

    //get_atom_types
    {pypdffit2_get_atom_types__name__, pypdffit2_get_atom_types,
     METH_VARARGS, pypdffit2_get_atom_types__doc__},

    //num_phases
    {pypdffit2_num_phases__name__, pypdffit2_num_phases,
     METH_VARARGS, pypdffit2_num_phases__doc__},

    //num_datasets
    {pypdffit2_num_datasets__name__, pypdffit2_num_datasets,
     METH_VARARGS, pypdffit2_num_datasets__doc__},

    //phase_fractions
    {pypdffit2_phase_fractions__name__, pypdffit2_phase_fractions,
     METH_VARARGS, pypdffit2_phase_fractions__doc__},

    //redirect_stdout
    {pypdffit2_redirect_stdout__name__, pypdffit2_redirect_stdout,
     METH_VARARGS, pypdffit2_redirect_stdout__doc__},

    //restore_stdout
    {pypdffit2_restore_stdout__name__, pypdffit2_restore_stdout,
     METH_VARARGS, pypdffit2_restore_stdout__doc__},

    //is_element
    {pypdffit2_is_element__name__, pypdffit2_is_element,
     METH_VARARGS, pypdffit2_is_element__doc__},

// Sentinel
    {0, 0}
};

// End of file
