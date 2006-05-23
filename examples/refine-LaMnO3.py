#
#  This is a template for refining the LaMnO3 strucuture.
#
#
# Last updated: Xiangyun  04/24/2003
# Emil addopted 02/08/2005
#
#
from pdffit import *

import sys

#fstruct = "LaMnO3"
#fdata = "x036t010q35"

fstruct = sys.argv[1]
fdata = sys.argv[2]

reset()

#
# READ the structure file and the Observed PDF file:
read_struct(fstruct+".stru")
read_data(fdata+".gr",N,35.0,0.0)
#
# experiment settings:
#
#
pdfrange(1,1.6,15.0)
#
#
# Set refinement parameters:
#
# Setup lattice parameters
#
constrain(lat(1),1)
constrain(lat(2),2)
constrain(lat(3),3)
#
setpar(1,lat(1))
setpar(2,lat(2))
setpar(3,lat(3))
#
# setpar(3,3.5*getvar(lat(3)))
#
# Set fudge factors
#
constrain(dscale,100)
setpar(100,1.00)
#
constrain(qsig,101)
setpar(101,0.0343)
fixpar(101)
#
constrain(delta,102)
setpar(102,0.00)
#
#
# Refinement step (1)
#
refine()
#
#####################################################
# Positions :(constrained)
#####################################################
#
# LA/CA on (4c)
#
constrain(x(1),"1.0+@21")
constrain(y(1),"0.0+@22")
constrain(x(2),"0.5+@21")
constrain(y(2),"0.5-@22")
constrain(x(3),"0.0-@21")
constrain(y(3),"1.0-@22")
constrain(x(4),"0.5-@21")
constrain(y(4),"0.5+@22")
#
constrain(x(5),"1.0+@21")
constrain(y(5),"0.0+@22")
constrain(x(6),"0.5+@21")
constrain(y(6),"0.5-@22")
constrain(x(7),"0.0-@21")
constrain(y(7),"1.0-@22")
constrain(x(8),"0.5-@21")
constrain(y(8),"0.5+@22")
#
setpar(21,getvar(x(1))-1.00)
setpar(22,y(1))
#
#
# O on (4c)
#
constrain(x(13),"0.0+@23")
constrain(y(13),"0.0+@24")
constrain(x(14),"0.5+@23")
constrain(y(14),"0.5-@24")
constrain(x(15),"1.0-@23")
constrain(y(15),"1.0-@24")
constrain(x(16),"0.5-@23")
constrain(y(16),"0.5+@24")
#
setpar(23,x(13))
setpar(24,y(13))
#
# Refinement step (3)
#
refine()
#
#
# O on (8d)
#
constrain(x(17),"0.0+@25")
constrain(y(17),"0.0+@26")
constrain(z(17),"0.0+@27")

constrain(x(18),"-0.5+@25")
constrain(y(18),"0.5-@26")
constrain(z(18),"1.0-@27")

constrain(x(19),"1.0-@25")
constrain(y(19),"1.0-@26")
constrain(z(19),"0.5+@27")

constrain(x(20),"0.5-@25")
constrain(y(20),"0.5+@26")
constrain(z(20),"0.5-@27")

constrain(x(21),"1.0-@25")
constrain(y(21),"1.0-@26")
constrain(z(21),"1.0-@27")

constrain(x(22),"0.5-@25")
constrain(y(22),"0.5+@26")
constrain(z(22),"0.0+@27")

constrain(x(23),"0.0+@25")
constrain(y(23),"0.0+@26")
constrain(z(23),"0.5-@27")

constrain(x(24),"-0.5+@25")
constrain(y(24),"0.5-@26")
constrain(z(24),"0.5+@27")
#
setpar(25,x(17))
setpar(26,y(17))
setpar(27,z(17))
#
# Refinement step (4)
#
refine()
#
#
# Set up thermals
#
for i in range (1,9):
   constrain(u11(i),7)
   constrain(u22(i),7)
   constrain(u33(i),7)
#
for i in range (9,13):
   constrain(u11(i),8)
   constrain(u22(i),8)
   constrain(u33(i),8)
#
for i in range (13,17):
   constrain(u11(i),9)
   constrain(u22(i),9)
   constrain(u33(i),9)
#
for i in range (17,25):
   constrain(u11(i),10)
   constrain(u22(i),10)
   constrain(u33(i),10)
#
setpar(7,u11(1))
setpar(8,u11(9))
setpar(9,u11(13))
setpar(10,u11(17))
#setpar(8,0.002)
#setpar(9,0.003)
#
#
# Refinement step (2)
#
refine()
#
# save all datas:
#
fname = fstruct+"-"+fdata
save_pdf(1,fname+".calc")
save_dif(1,fname+".dif")
save_struct(1,fname+".rstr")
save_res(fname+".res")
#
