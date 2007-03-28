#####################################################################
# Ni- (FCC) Refinement - average structure
######################################################################
from diffpy.pdffit2.PdfFit import *

P = PdfFit()
P.read_struct("Ni.stru")
P.alloc('X',30.0,0.03,1.5,20,1001)
P.calc()
P.save_pdf(1,"Ni.calc")
