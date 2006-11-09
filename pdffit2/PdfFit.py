########################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farros, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""PdfFit class for fitting pdf data to a model."""

# version
__id__ = "$Id$"

import pdffit2

__intro_message__ = """
        ****************************************************************
        *               P D F F I T   Version   %(version)-14s         *
        *                                       %(date)-11s            *
        * ------------------------------------------------------------ *
        *   (c) 1998-2006 trustees of the Michigan State University.   *
        *   Authors:                                                   *
        *       Thomas Proffen  -   Email: tproffen@lanl.gov           *
        *       Jacques Bloch   -   Email: bloch@pa.msu.edu            *
        *       Chris Farrow    -   Email: farrowch@msu.edu            *
        *       Pavol Juhas     -   Email: juhas@pa.msu.edu            *
        *       Simon Billinge  -   Email: billinge@pa.msu.edu         *
        ****************************************************************
"""

class PdfFit(object):
    """Create PdfFit object."""
    import sys

    # constants and enumerators from pdffit.h:
    # selection of all atoms
    selalias = { 'ALL' : -1 }
    # constraint type identifiers
    FCON = { 'USER' : 0, 'IDENT' : 1, 'FCOMP' : 2, 'FSQR' : 3 }
    # scattering type identifiers
    Sctp = { 'X' : 0, 'N' : 1 }

    def _exportAll(self, namespace):
        """ _exportAll(self, namespace) --> Export all 'public' class methods
            into namespace.

        This function allows for a module-level PdfFit object which doesn't have
        to be referenced when calling a method. This function makes old (python)
        scripts compatible with this class. At the top of the script, create a
        pdffit object, and then call this method. Usually, namespace =
        sys.modules[__name__].__dict__.


        """
        # string aliases (var = "var")
        for a in self.selalias.keys() + self.FCON.keys() + self.Sctp.keys():
            exec("%s = %r" % (a, a), namespace)
        public = [ a for a in dir(self) if "__" not in a and a not in
                ["_handle", "_exportAll", "selalias", "FCON", "Sctp" ] ]
        import sys
        for funcname in public:
            namespace[funcname] = getattr(self, funcname)
        return

    def intro(self):
        """Show introductory message.
        """
        from version import __version__, __date__
        d = { 'version' : __version__,  'date' : __date__ }
        msg = __intro_message__ % d
        print msg
        return

    def read_struct(self, struct):
        """read_struct(struct) --> Read structure from file into memory.

        struct  -- name of file from which to read structure

        Raises:
            pdffit2.calculationError when a lattice cannot be created from the
            given structure
            pdffit2.structureError when a structure file is malformed
            IOError when the file cannot be read from the disk
        """
        pdffit2.read_struct(self._handle, struct)
        self.stru_files.append(struct)
        return


    def read_struct_string(self, struct, name = ""):
        """read_struct_string(struct, name = "") --> Read structure from
        a string into memory.

        struct  -- string containing the contents of the structure file
        name    -- tag with which to label structure

        Raises:
            pdffit2.calculationError when a lattice cannot be created from the
            given structure
            pdffit2.structureError when a structure file is malformed
        """
        pdffit2.read_struct_string(self._handle, struct)
        self.stru_files.append(name)
        return


    def read_data(self, data, stype, qmax, sigmaq):
        """read_data(data, stype, qmax, sigmaq) --> Read pdf data from file into
        memory.

        data    -- name of file from which to read data
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        sigmaq  -- instrumental Q-resolution factor

        Raises: IOError when the file cannot be read from disk
        """
        pdffit2.read_data(self._handle, data, stype, qmax, sigmaq)
        self.data_files.append(data)
        if self.data_server:
            self.data_server.setDataDescriptor(data)
        return


    def read_data_string(self, data, stype, qmax, sigmaq, name = ""):
        """read_data_string(data, stype, qmax, sigmaq, name = "") --> Read
        pdf data from a string into memory.

        data    -- string containing the contents of the data file
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        sigmaq  -- instrumental Q-resolution factor
        name    -- tag with which to label data
        """
        pdffit2.read_data_string(self._handle, data, stype, qmax,
                sigmaq, name)
        name = data
        self.data_files.append(name)
        if self.data_server:
            self.data_server.setDataDescriptor(name)
        return


    def read_data_lists(self, stype, qmax, sigmaq, r_data, Gr_data,
            dGr_data = None, name = "list"):
        """read_data_lists(stype, qmax, sigmaq, r_data, Gr_data, dGr_data =
        None, name = "list") --> Read pdf data into memory from lists.

        All lists must be of the same length.
        stype       -- 'X' (xray) or 'N' (neutron)
        qmax        -- Q-value cutoff used in PDF calculation.
                       Use qmax=0 to neglect termination ripples.
        sigmaq      -- instrumental Q-resolution factor
        r_data      -- list of r-values
        Gr_data     -- list of G(r) values
        dGr_data    -- list of G(r) uncertainty values
        name        -- tag with which to label data

        Raises: ValueError when the data lists are of different length
        """
        pdffit2.read_data_arrays(self._handle, stype, qmax, sigmaq,
                r_data, Gr_data, dGr_data, name)
        self.data_files.append(name)
        return


    def pdfrange(self, iset, rmin, rmax):
        """pdfrange(iset, rmin, rmax) --> Set the range of the fit.

        iset    -- data set to consider
        rmin    -- minimum r-value of fit
        rmax    -- maximum r-value of fit

        Raises: ValueError for bad input values
        """
        pdffit2.pdfrange(self._handle, iset, rmin, rmax)
        return


    def reset(self):
        """reset() --> Clear all stored fit, structure, and parameter data."""
        self.stru_files = []
        self.data_files = []
        pdffit2.reset(self._handle);
        return


    def alloc(self, stype, qmax, sigmaq, rmin, rmax, bin):
        """alloc(stype, qmax, sigmaq, rmin, rmax, bin) --> Allocate space for a
        PDF calculation.

        The structure from which to calculate the PDF must first be imported with
        the read_struct() or read_struct_string() method.
        stype   -- 'X' (xray) or 'N' (neutron)
        qmax    -- Q-value cutoff used in PDF calculation.
                   Use qmax=0 to neglect termination ripples.
        sigmaq  -- instrumental Q-resolution factor
        rmin    -- minimum r-value of calculation
        rmax    -- maximum r-value of calculation
        bin     -- number of data points in calculation

        Raises:
            ValueError for bad input values
            pdffit.unassignedError when no structure has been loaded
        """
        pdffit2.alloc(self._handle, stype, qmax, sigmaq, rmin,
                rmax, bin)
        return


    def calc(self):
        """calc() --> Calculate the PDF of the imported structure.

        Space for the calculation must first be allocated with the alloc()
        method.

        Raises:
            pdffit2.calculationError when allocated space cannot
            accomodate calculation
            pdffit.unassignedError when space for calculation has not been
            allocated
        """
        pdffit2.calc(self._handle)
        return


    def refine(self, toler = 0.00000001):
        """refine(toler = 0.00000001) --> Fit the theory to the imported data.

        toler   --  tolerance of the fit

        Raises:
            pdffit2.calculationError when the model pdf cannot be calculated
            pdffit2.constraintError when refinement fails due to bad
            constraint
            pdffit2.unassigedError when a constraint used but never initialized
            using setpar()
        """
        step = 0
        finished = 0
        while not finished:
            finished = pdffit2.refine_step(self._handle, toler)
            step += 1
            if self.data_server:
                import time
                self.data_server.update(self, step, finished)
                time.sleep(0.1)
        return


    def refine_step(self, toler = 0.00000001):
        """refine_step(toler = 0.00000001) --> Run a single step of the fit.

        toler   --  tolerance of the fit

        Raises:
            pdffit2.calculationError when the model pdf cannot be calculated
            pdffit2.constraintError when refinement fails due to bad
            constraint
            pdffit2.unassigedError when a constraint used but never initialized
            using setpar()

        Returns: 1 (0) if refinement is (is not) finished
        """
        self.finished = pdffit2.refine_step(self._handle, toler)
        return self.finished


    def save_pdf(self, iset, fname):
        """save_pdf(iset, fname) --> Save calculated or fitted PDF to file.

        iset    -- data set to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_pdf(self._handle, iset, fname)
        return


    def save_pdf_string(self, iset):
        """save_pdf_string(iset) --> Save calculated or fitted PDF to string.

        iset    -- data set to save

        Raises:
            pdffit2.unassignedError if the data set is undefined

        Returns: string containing contents of save file
        """
        pdffilestring = pdffit2.save_pdf(self._handle, iset, "")
        return pdffilestring


    def save_dif(self, iset, fname):
        """save_dif(iset, fname) --> Save data and fitted PDF difference to
        file.

        iset    -- data set to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_dif(self._handle, iset, fname)
        return


    def save_dif_string(self, iset):
        """save_dif_string(iset) --> Save data and fitted PDF difference to
        string.

        iset    -- data set to save

        Raises:
            pdffit2.unassignedError if the data set is undefined

        Returns: string containing contents of save file
        """
        diffilestring = pdffit2.save_dif(self._handle, iset, "")
        return diffilestring


    def save_res(self, fname):
        """save_res(fname) --> Save fit-specific data to file.

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if there is no refinement data to save
        """
        pdffit2.save_res(self._handle, fname)
        return


    def save_res_string(self):
        """save_res_string() --> Save fit-specific data to a string.

        Raises:
            pdffit2.unassignedError if there is no refinement data to save

        Returns: string containing contents of save file
        """
        resfilestring = pdffit2.save_res(self._handle, "")
        return resfilestring


    def save_struct(self, ip, fname):
        """save_struct(ip, fname) --> Save structure resulting from fit
        to file.

        ip    -- phase to save

        Raises:
            IOError if file cannot be saved
            pdffit2.unassignedError if the data set is undefined
        """
        pdffit2.save_struct(self._handle, ip, fname)
        return


    def save_struct_string(self, ip):
        """save_struct(ip) --> Save structure resulting from fit to string.

        ip    -- phase to save

        Raises:
            pdffit2.unassignedError if the data set is undefined

        Returns: string containing contents of save file
        """
        structfilestring = pdffit2.save_struct(self._handle, ip, "")
        return structfilestring


    def show_struct(self, ip):
        """show_struct(ip) --> Print structure resulting from fit.

        ip    -- phase to display

        Raises: pdffit2.unassignedError if the phase is undefined
        """
        pdffit2.show_struct(self._handle, ip)
        return


    def constrain(self, var, par, fcon=None):
        """constrain(var, par[, fcon]) --> Constrain a variable to a parameter.

        A variable can be constrained to a number or equation string.
        var     -- variable to constrain, such as x(1)
        par     -- parameter which to constrain the variable. This can be
                   an integer or an equation string containing a reference
                   to another parameter. Equation strings use standard c++
                   syntax. The value of a constrained parameter is accessed
                   as @p in an equation string, where p is the parameter.
                   e.g.
                   >>>  constrain(x(1), 1)
                   >>>  constrain(x(2), "0.5+@1")
        fcon    -- 'USER', 'IDENT', 'FCOMP', or 'FSQR'
                   this is an optional parameter, and I don't know how it is
                   used!

        Raises:
            pdffit2.constraintError if a constraint is bad
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        import types
        var_ref = self.__getRef(var)
        if fcon:
            pdffit2.constrain_int(self._handle, var_ref, par, self.FCON[fcon])
        elif type(par) == types.StringType:
            pdffit2.constrain_str(self._handle, var_ref, par)
        else:
            pdffit2.constrain_int(self._handle, var_ref, par)
        return


    def setpar(self, par, val):
        """setpar(par, val) --> Set value of constrained parameter.

        val     --  Either a numerical value or a reference to a variable

        Raises:
            pdffit2.unassignedError when variable is yet to be assigned
        """
        # people do not use parenthesis, e.g., "setpar(3, qsig)"
        # in such case val is a reference to PdfFit method
        import types
        if type(val) is types.MethodType:
            val = val()
        try:
            val = float(val)
            pdffit2.setpar_dbl(self._handle, par, val)
        except ValueError:
            var_ref = self.__getRef(val)
            pdffit2.setpar_RV(self._handle, par, var_ref)
        return


    def setvar(self, var, val):
        """setvar(var, val) --> Set the value of a variable.

        Raises:
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_ref = self.__getRef(var)
        pdffit2.setvar(self._handle, var_ref, val)
        return


    def getvar(self, var):
        """getvar(var) --> Get stored value of a variable.

        Raises:
            pdffit2.unassignedError if variable does not yet exist
            ValueError if variable index does not exist (e.g. lat(7))
        """
        var_ref = self.__getRef(var)
        retval = pdffit2.getvar(self._handle, var_ref)
        return retval


    def getrw(self):
        """getrw() --> Get goodness of fit value, rw."""
        rw = pdffit2.getrw(self._handle)
        return rw


    def getR(self):
        """getR() --> Get r-points used in the fit.

        This function should only be called after data has been loaded or
        calculated. Before a refinement, the list of r-points will reflect the
        data. Afterwords, they will reflect the fit range.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of equidistance r-points used in fit.
        """
        R = pdffit2.getR(self._handle)
        return R


    def getpdf_fit(self):
        """getpdf_fit() --> Get fitted PDF.

        This function should only be called after a refinement or refinement
        step has been done.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of fitted points, equidistant in r.
        """
        pdfdata = pdffit2.getpdf_fit(self._handle)
        return pdfdata


    def getpdf_obs(self):
        """getpdf_obs() --> Get observed PDF.

        This function should only be called after data has been loaded or
        calculated. Before a refinement, the list of r-points will reflect the
        data. Afterwords, they will reflect the fit range.

        Raises: pdffit2.unassignedError if no data exists

        Returns: List of data points, equidistant in r.
        """
        pdfdata = pdffit2.getpdf_obs(self._handle)
        return pdfdata


    def get_atoms(self, ip=None):
        """get_atoms() --> Get element symbols of all atoms in the structure.

        ip -- index of phase to get the elements from (starting from 1)
              when ip is not given, use current phase

        This function should only be called after a structure has been loaded.

        Raises: pdffit2.unassignedError if no structure exists

        Returns: List of atom names in structure.
        """
        if ip is None:  rv = pdffit2.get_atoms(self._handle)
        else:           rv = pdffit2.get_atoms(self._handle, ip)
        return rv


    def get_atom_types(self, ip=None):
        """get_atom_types() --> Ordered unique element symbols in the structure.

        ip -- index of phase to get the elements from (starting from 1)
              when ip is not given, use current phase

        This function should only be called after a structure has been loaded.

        Raises: pdffit2.unassignedError if no structure exists

        Returns: List of unique atom symbols as they occur in structure.
        """
        if ip is None:  rv = pdffit2.get_atom_types(self._handle)
        else:           rv = pdffit2.get_atom_types(self._handle, ip)
        return rv


    def getpar(self, par):
        """getpar(par) --> Get value of parameter.

        Raises: ValueError if parameter does not exists
        """
        return pdffit2.getpar(self._handle, par)


    def fixpar(self, par):
        """fixpar(par) --> Fix a parameter.

        Fixed parameters are not fitted in a refinement. Passed parameter
        can be 'ALL', in which case all parameters are fixed.

        Raises: pdffit.unassignedError when parameter has not been assigned
        """
        import types
        if type(par) in types.StringTypes and par.upper() in self.selalias:
            par = self.selalias[par.upper()]
        pdffit2.fixpar(self._handle, par)
        return


    def freepar(self, par):
        """freepar(par) --> Free a parameter.

        Freed parameters are fitted in a refinement. Passed parameter
        can be 'ALL', in which case all parameters are freed.

        Raises: pdffit.unassignedError when parameter has not been assigned
        """
        import types
        if type(par) in types.StringTypes and par.upper() in self.selalias:
            par = self.selalias[par.upper()]
        pdffit2.freepar(self._handle, par)
        return


    def setphase(self, ip):
        """setphase(ip) --> Switch to phase ip.

        All parameters assigned after this method is called refer only to the
        current phase.

        Raises: pdffit.unassignedError when phase does not exist
        """
        pdffit2.setphase(self._handle, ip)
        return


    def setdata(self, iset):
        """setdata(iset) --> Set the data in focus.

        Raises: pdffit.unassignedError when data set does not exist
        """
        pdffit2.setdata(self._handle, iset)
        return


    def psel(self, ip):
        """psel(ip) --> Include phase ip in calculation of total PDF

        psel('ALL')     selects all phases for PDF calculation.

        Raises: pdffit2.unassignedError if selected phase does not exist
        """
        import types
        if type(ip) in types.StringTypes and ip.upper() in self.selalias:
            ip = self.selalias[ip.upper()]
        pdffit2.psel(self._handle, ip)
        return


    def pdesel(self, ip):
        """pdesel(ip) --> Exclude phase ip from calculation of total PDF.

        pdesel('ALL')   excludes all phases from PDF calculation.

        Raises: pdffit2.unassignedError if selected phase does not exist
        """
        import types
        if type(ip) in types.StringTypes and ip.upper() in self.selalias:
            ip = self.selalias[ip.upper()]
        pdffit2.pdesel(self._handle, ip)
        return


    def isel(self, ip, element):
        """isel(ip, element) --> Include atoms of given element from phase
        ip as first in pair distance evaluation.  Used for calculation of
        partial PDF.  When element is 'ALL', all elements are included as
        first-in-pair.

        ip      -- phase index starting at 1
        element -- integer index of atom type to be selected starting at 1 or
                   element symbol, such as "Na", "CL", "o", or "ALL"

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if selected atom type does not exist
        """
        import types
        if type(element) is types.IntType:
            atom_types = self.get_atom_types(ip)
            if not 0 <= element-1 < len(atom_types):
                raise ValueError, 'Invalid atom type index %i' % element
            element = atom_types[element - 1]
        # element should be string here
        if element.upper() == 'ALL':
            self.selectAll(ip, 'i')
        else:
            self.selectAtomType(ip, 'i', element, True)
        return


    def idesel(self, ip, element):
        """idesel(ip, element) --> Do not use atoms of given element from
        phase ip as first in pair distance evaluation.  Used for calculation
        of partial PDF.  When element is 'ALL', all atom types are excluded
        from first-in-pair.

        ip      -- phase index starting at 1
        element -- integer index of atom type to be excluded starting at 1 or
                   element symbol, such as "Na", "CL", "o" or "ALL"

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if selected atom type does not exist
        """
        import types
        if type(element) is types.IntType:
            atom_types = self.get_atom_types(ip)
            if not 0 <= element-1 < len(atom_types):
                raise ValueError, 'Invalid atom type index %i' % element
            element = atom_types[element - 1]
        # element should be string here
        if element.upper() == 'ALL':
            self.selectNone(ip, 'i')
        else:
            self.selectAtomType(ip, 'i', element, False)
        return


    def jsel(self, ip, element):
        """jsel(ip, element) --> Include atoms of given element from phase
        ip as second in pair distance evaluation.  Used for calculation of
        partial PDF.  When element is 'ALL', all atom types are included as
        second-in-pair.

        ip      -- phase index starting at 1
        element -- integer index of atom type to be selected starting at 1 or
                   element symbol, such as "Na", "CL", "o", or "ALL"

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if selected atom type does not exist
        """
        import types
        if type(element) is types.IntType:
            atom_types = self.get_atom_types(ip)
            if not 0 <= element-1 < len(atom_types):
                raise ValueError, 'Invalid atom type index %i' % element
            element = atom_types[element - 1]
        # element should be string here
        if element.upper() == 'ALL':
            self.selectAll(ip, 'j')
        else:
            self.selectAtomType(ip, 'j', element, True)
        return


    def jdesel(self, ip, element):
        """jdesel(ip, element) --> Do not use atoms of given element from
        phase ip as second in pair distance evaluation.  Used for calculation
        of partial PDF.  When element is 'ALL', all atom types are excluded
        from second-in-pair.

        ip      -- phase index starting at 1
        element -- integer index of atom type to be excluded starting at 1 or
                   element symbol, such as "Na", "CL", "o" or "ALL"

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if selected atom type does not exist
        """
        import types
        if type(element) is types.IntType:
            atom_types = self.get_atom_types(ip)
            if not 0 <= element-1 < len(atom_types):
                raise ValueError, 'Invalid atom type index %i' % element
            element = atom_types[element - 1]
        # element should be string here
        if element.upper() == 'ALL':
            self.selectNone(ip, 'j')
        else:
            self.selectAtomType(ip, 'j', element, False)
        return

    def selectAtomType(self, ip, ijchar, symbol, flag):
        """Mark given atom type in phase ip as included or excluded in
         first or second in pair for distance evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair
        symbol  -- element symbol
        flag    -- bool flag, True for selection, False for exclusion

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError for invalid value of ijchar
        """
        pdffit2.selectAtomType(ip, ijchar, symbol, flag)
        return

    def selectAtomIndex(self, ip, ijchar, aidx, flag):
        """Mark atom of give index in phase ip as included or excluded in
         first or second in pair for distance evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair
        aidx    -- integer index of atom starting at 1
        flag    -- bool flag, True for selection, False for exclusion

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if atom index or ijchar are invalid
        """
        pdffit2.selectAtomIndex(ip, ijchar, aidx, flag)
        return

    def selectAll(self, ip, ijchar):
        """Include all atoms in phase ip in first or second pair of distance
        evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if ijchar is invalid
        """
        pdffit2.selectAll(ip, ijchar)
        return


    def selectNone(self, ip, ijchar):
        """Exclude all atoms in phase ip from first or second pair of distance
        evaluation.

        ip      -- phase index starting at 1
        ijchar  -- 'i' or 'j' for first or second in pair

        Raises:
            pdffit2.unassignedError if selected phase does not exist
            ValueError if ijchar is invalid
        """
        pdffit2.selectNone(ip, ijchar)
        return

    def bang(self, ia, ja, ka):
        """bang(ia, ja, ka) --> Get the bond angle defined by atoms ia, ja, ka.

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded
        """
        ba = pdffit2.bang(self._handle, ia, ja, ka)
        return ba


    def blen(self, *args):
        """blen(ia, ja) --> Get length of bond defined by atoms ia and ja.

        blen(a1, a2, lb, ub) --> Print length of all a1-a2 bonds in range
        [lb,ub], where a1 and a2 are element names.  Either a1 or a2 can
        be the string "ALL", in which all atom types are used for that end
        of the calculated bonds.

        Raises: ValueError if selected atom(s) does not exist
                pdffit.unassignedError when no structure has been loaded
        """
        if len(args)==2:
            res = pdffit2.blen_atoms(self._handle, args[0], args[1])
            return res
        elif len(args)==4:
            a1 = args[0]
            a2 = args[1]
            lb = args[2]
            ub = args[3]
            pdffit2.blen_types(self._handle, a1, a2, lb, ub)
            return
        else:
            message = "blen() takes 3 or 5 arguments (%i given)" % (len(args)+1)
            raise TypeError, message
        return


    def show_scat(self, stype):
        """show_scat(stype) --> Print scattering length for all atoms in
        the current phase.

        stype -- 'X' (xray) or 'N' (neutron).

        Raises: pdffit2.unassignedError if no phase exists
        """
        print self.get_scat_string(stype)
        return


    def get_scat_string(self, stype):
        """get_scat_string(stype) --> Get string with scattering factors
        of all atoms in the current phase.

        stype -- 'X' (xray) or 'N' (neutron).

        Raises: pdffit2.unassignedError if no phase exists

        Returns: string with all scattering factors.
        """
        return pdffit2.get_scat_string(self._handle, stype)


    def set_scat(self, stype, element, value):
        """set_scat(stype, element, value) --> Set custom scattering factor
        for given element.

        stype   -- 'X' (xray) or 'N' (neutron).
        element -- integer index of atom type or case-insensitive symbol,
                   such as "Na" or "CL"
        value   -- custom value of scattering factor

        Raises:
            pdffit2.unassignedError if no phase exists
            ValueError if atom type does not exist
        """
        pdffit2.set_scat(self._handle, stype, element, value)
        return


    def reset_scat(self, stype, element):
        """reset_scat(stype, element) --> Reset scattering factor for given
        element to standard value.

        stype   -- 'X' (xray) or 'N' (neutron).
        element -- integer index of atom type or case-insensitive symbol,
                   such as "Na" or "CL"
        Raises:
            pdffit2.unassignedError if no phase exists
            ValueError if atom type does not exist
        """
        import types
        if type(element) is types.IntType:
            atom_types = self.get_atom_types()
            if not 0 <= element-1 < len(atom_types):
                raise ValueError, 'Invalid atom type index %i' % 77 #element
            element = atom_types[element - 1]
        pdffit2.reset_scat(self._handle, stype, element)
        return

    def num_atoms(self):
        """num_atoms() --> Get number of atoms in current phase.

        Raises: pdffit2.unassignedError if no atoms exist
        """
        return pdffit2.num_atoms(self._handle)


    # Begin refineable variables.

    def lat(self, n):
        """lat(n) --> Get reference to lattice variable n.

        n can be an integer or a string representing the lattice variable.
        1 <==> 'a'
        2 <==> 'b'
        3 <==> 'c'
        4 <==> 'alpha'
        5 <==> 'beta'
        6 <==> 'gamma'
        """
        LatParams = { 'a':1, 'b':2, 'c':3, 'alpha':4, 'beta':5, 'gamma':6 }
        import types
        if type(n) is types.StringType:
            n = LatParams[n]
        return "lat(%i)" % n


    def x(self, i):
        """x(i) --> Get reference to x-value of atom i."""
        return "x(%i)" % i


    def y(self, i):
        """y(i) --> Get reference to y-value of atom i."""
        return "y(%i)" % i


    def z(self, i):
        """z(i) --> Get reference to z-value of atom i."""
        return "z(%i)" % i


    def u11(self, i):
        """u11(i) --> Get reference to U(1,1) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u11(%i)" % i


    def u22(self, i):
        """u22(i) --> Get reference to U(2,2) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u22(%i)" % i


    def u33(self, i):
        """u33(i) --> Get reference to U(3,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u33(%i)" % i


    def u12(self, i):
        """u12(i) --> Get reference to U(1,2) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u12(%i)" % i


    def u13(self, i):
        """u13(i) --> Get reference to U(1,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u13(%i)" % i


    def u23(self, i):
        """u23(i) --> Get reference to U(2,3) for atom i.

        U is the anisotropic thermal factor tensor.
        """
        return "u23(%i)" % i


    def occ(self, i):
        """occ(i) --> Get reference to occupancy of atom i."""
        return "occ(%i)" % i


    def pscale(self):
        """pscale() --> Get reference to pscale.

        pscale is the fraction of the total structure that the current phase
        represents.
        """
        return "pscale"


    def pfrac(self):
        """pfrac() --> same as pscale.

        pscale is the fraction of the total structure that the current phase
        represents.
        """
        return self.pscale()


    def srat(self):
        """srat() --> Get reference to sigma ratio.

        The sigma ratio determines the reduction in the Debye-Waller factor for
        distances below rcut.
        """
        return "srat"


    def delta1(self):
        """delta1() --> Get reference to 1/R peak sharpening factor.
        """
        return "delta1"


    def delta2(self):
        """delta2() --> Reference to (1/R^2) sharpening factor.
        The phenomenological correlation constant in the Debye-Waller factor.
        The (1/R^2) peak sharpening factor.
        """
        return "delta2"


    def delta(self):
        """delta() --> Reference to (1/R^2) sharpening factor.  Same as delta2.
        The phenomenological correlation constant in the Debye-Waller factor.
        The (1/R^2) peak sharpening factor.
        """
        return self.delta2()


    def gamma(self):
        """gamma() --> Same as delta1().

        gamma varible is deprecated, it has been renamed to delta1.

        1/R peak sharpening factor.
        """
        print >> sys.stderr, "Variable gamma is deprecated, use delta1"
        return self.delta1()


    def dscale(self):
        """dscale() --> Get reference to dscale.

        The data scale factor.
        """
        return "dscale"


    def qsig(self):
        """qsig() --> Get reference to qsig.

        instrument q-resolution factor.
        """
        return "qsig"


    def qalp(self):
        """qalp() --> Get reference to qalp.

        Quadratic peak sharpening factor.
        """
        return "qalp"


    def rcut(self):
        """rcut() --> Get reference to rcut.

        rcut is the value of r below which peak sharpening, defined by the sigma
        ratio (srat), applies.
        """
        return "rcut"


    # End refineable variables.

    def __init__(self, data_server=None):

        self.data_server = data_server
        self.stru_files = []
        self.data_files = []

        self._handle = pdffit2.create()
        self.intro()
        return


    def __getRef(self, var_string):
        """Return the actual reference to the variable in the var_string.

        This function must be called before trying to actually reference an
        internal variable. See the constrain method for an example.

        Raises:
            pdffit2.unassignedError if variable is not yet assigned
            ValueError if variable index does not exist (e.g. lat(7))
        """
        # people do not use parenthesis in their scripts, e.g., "getvar(qsig)"
        # in such case var_string is a reference to PdfFit method
        import types
        if type(var_string) is types.MethodType:
            var_string = var_string()
        arg_int = None
        try:
            method_string, arg_string = var_string.split("(")
            method_string = method_string.strip()
            arg_int = int(arg_string.strip(")").strip())
        except ValueError: #There is no arg_string
            method_string = var_string.strip()

        f = getattr(pdffit2, method_string)
        if arg_int is None:
            retval = f(self._handle)
            return retval
        else:
            retval = f(self._handle, arg_int)
            return retval


    # End of class PdfFit


# End of file
