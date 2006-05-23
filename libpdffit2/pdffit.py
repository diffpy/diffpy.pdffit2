# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _pdffit

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


X = _pdffit.X
N = _pdffit.N
USER = _pdffit.USER
IDENT = _pdffit.IDENT
FCOMP = _pdffit.FCOMP
FSQR = _pdffit.FSQR
ALL = _pdffit.ALL

read_struct = _pdffit.read_struct

read_data = _pdffit.read_data

pdfrange = _pdffit.pdfrange

reset = _pdffit.reset

alloc = _pdffit.alloc

calc = _pdffit.calc

refine = _pdffit.refine

save_pdf = _pdffit.save_pdf

save_dif = _pdffit.save_dif

save_res = _pdffit.save_res

save_struct = _pdffit.save_struct

show_struct = _pdffit.show_struct

getpar = _pdffit.getpar

fixpar = _pdffit.fixpar

freepar = _pdffit.freepar

setphase = _pdffit.setphase

setdata = _pdffit.setdata

psel = _pdffit.psel

pdesel = _pdffit.pdesel

isel = _pdffit.isel

idesel = _pdffit.idesel

jsel = _pdffit.jsel

jdesel = _pdffit.jdesel

bang = _pdffit.bang

show_scat = _pdffit.show_scat

reset_scat = _pdffit.reset_scat

lat = _pdffit.lat

x = _pdffit.x

y = _pdffit.y

z = _pdffit.z

u11 = _pdffit.u11

u22 = _pdffit.u22

u33 = _pdffit.u33

u12 = _pdffit.u12

u13 = _pdffit.u13

u23 = _pdffit.u23

occ = _pdffit.occ
from math import *
import signal, os
import sys

def idle(expr):
    return


sys.ps1 = "pdffit2> "
sys.displayhook = idle

signal.signal(signal.SIGINT,signal.SIG_DFL)

pfrac = _pdffit.cvar.pfrac
srat = _pdffit.cvar.srat
delta = _pdffit.cvar.delta
gamma = _pdffit.cvar.gamma
dscale = _pdffit.cvar.dscale
qsig = _pdffit.cvar.qsig
qalp = _pdffit.cvar.qalp
rcut = _pdffit.cvar.rcut



constrain = _pdffit.constrain

setpar = _pdffit.setpar

setvar = _pdffit.setvar

getvar = _pdffit.getvar

blen = _pdffit.blen

set_scat = _pdffit.set_scat
cvar = _pdffit.cvar

