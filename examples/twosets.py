# FIXME: this file is obsolete and needs to be cleaned up or removed
# PJ, 2007-03-27

from pdffit import *


#####################################################################
# Ni- (FCC) Refinement - average structure
######################################################################

read_struct("Ni.stru")

###############################################################
# Experimental and lattice parameters
###############################################################

constrain(lat(1), 1)
constrain(lat(2), 1)
constrain(lat(3), 1)

setpar(1, lat(1))

constrain(pscale, 20)
constrain(delta2, 22)

setpar(20, 0.989)
setpar(22, 0.0003)

constrain(sratio, 23)
setpar(23, 0.25)
setvar(rcut, 2.5)

###############################################################
# Temperature factors
###############################################################

for i in range(1, 5):
    constrain(u11(i), 4)
    constrain(u22(i), 4)
    constrain(u33(i), 4)

setpar(4, u11(1))

read_data("Ni_2-8.chi.gr", X, 30.0, 0.0)
constrain(qdamp, 21)
setpar(21, 0.03)
pdfrange(1, 1.5, 15)

# Read second data-set
read_data("Ni_2-8.chi.gr", X, 30.0, 0.0)
constrain(qdamp, 25)
setpar(25, 0.04)
pdfrange(2, 5.0, 19.99)

###############################################################
# Run fit
###############################################################

refine()

###############################################################
# Save data
###############################################################

#save_pdf(1, "ahm4.pdf")
save_struct(1, "twosets.stru")
save_res("twosets.res")

