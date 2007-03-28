#
#
from diffpy.pdffit2.PdfFit import *
import sys

P = PdfFit()
fstruct = "cu_bad"
fdata = "cu_q35"

P.reset()

#
# READ the structure file and the Observed PDF file:
P.read_struct(fstruct+".stru")
P.read_data(fdata+".gr",'N',35.0,0.0)
#
# experiment settings:
#
#
P.pdfrange(1,1.6,9.0)
#
#
# Setup lattice parameters
#
P.constrain(P.lat(1),1)
P.constrain(P.lat(2),2)
P.constrain(P.lat(3),3)
#
P.setpar(1,P.getvar(P.lat(1)))
P.setpar(2,P.getvar(P.lat(2)))
P.setpar(3,P.getvar(P.lat(3)))
#
#
# Set fudge factors
#
P.constrain(P.dscale(),100)
P.setpar(100,1.00)
#
P.constrain(P.qsig(),101)
P.setpar(101,0.0343)
P.fixpar(101)
#
P.constrain(P.delta2(),102)
P.setpar(102,0.00)
#
#
P.constrain(P.x(1),21, 'FCOMP')
P.constrain(P.y(1),21, 'FCOMP')
P.constrain(P.z(1),21, 'FCOMP')
P.setpar(21, P.getvar(P.x(1)))
#
#
# Set up thermals
# Same parameter for each
#
for i in range (1,9):
   P.constrain(P.u11(i),7)
   P.constrain(P.u22(i),7)
   P.constrain(P.u33(i),7)
#
P.setpar(7,P.getvar(P.u11(1)))
#
finished = 0
while not finished:
    finished = P.refine_step()
#P.refine()
#
# save all data:
#
fname = fstruct
P.save_pdf(1,fname+".calc")
P.save_dif(1,fname+".dif")
P.save_struct(1,fname+".rstr")
P.save_res(fname+".res")
#
