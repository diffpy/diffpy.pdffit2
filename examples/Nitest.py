from diffpy.pdffit2.PdfFit import *
P = PdfFit()
#####################################################################
# Ni- (FCC) Refinement - average structure
######################################################################

P.read_struct("Ni.stru")

###############################################################
# Experimental and P.lattice parameters
###############################################################

P.constrain(P.lat(1),1)
P.constrain(P.lat(2),1)
P.constrain(P.lat(3),1)

P.setpar(1,P.getvar(P.lat(1)))

P.constrain(P.pscale(),20)
P.constrain(P.delta2(),22)

P.setpar(20,0.0)
P.setpar(22,0.00)
P.fixpar(22)

###############################################################
# Temperature factors
###############################################################

for i in range(1, 5):
    P.constrain(P.u11(i),4)
    P.constrain(P.u22(i),4)
    P.constrain(P.u33(i),4)

P.setpar(4,P.getvar(P.u11(1)))

#####################################################################
# Ni- (FCC) Refinement - average structure -- Phase 2
######################################################################

P.read_struct("Ni.stru")

###############################################################
# Experimental and P.lattice parameters
###############################################################

P.constrain(P.lat(1),1)
P.constrain(P.lat(2),1)
P.constrain(P.lat(3),1)

P.constrain(P.pscale(),20,'FCOMP')
P.constrain(P.delta2(),22)


###############################################################
# Temperature factors
###############################################################

for i in range(1, 5):
    P.constrain(P.u11(i),4)
    P.constrain(P.u22(i),4)
    P.constrain(P.u33(i),4)


###############################################################
# Run fit
###############################################################
P.read_data("Ni.calc",'X',30.0,0.0)

P.constrain(P.qsig(),21)
P.setpar(21,0.03)
P.constrain(P.dscale(),25)
P.setpar(25,1.0)

P.pdfrange(1,1.5,20.0)

P.fixpar('ALL')
P.freepar(25)

P.refine()

###############################################################
# Save data
###############################################################

P.save_pdf(1,"2.pdf")
#save_dif(1,"ahm5.dif")
#save_struct(1,"ahm5.rstr")
P.save_res("2.res")

