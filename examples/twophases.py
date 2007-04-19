
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

setpar(22, 3.26)
fixpar(22)

constrain(sratio, 23)
setpar(23, 1.0)
#fixpar(23)
setvar(rcut, 2.5)

###############################################################
# Temperature factors
###############################################################

setpar(4, 0.005544)

for i in range(1, 5):
    constrain(u11(i), 4)
    constrain(u22(i), 4)
    constrain(u33(i), 4)


#######################
# read second phase P2
#######################
read_struct("Ni.stru")

setvar(lat(1), 3.52)
setvar(lat(2), 3.52)
setvar(lat(3), 3.52)
setvar(delta2, 3)
#constrain(delta2, 22)
constrain(pscale, 20, FCOMP)
setpar(20, 0.5)

constrain(sratio, 23)

###############################################################
# Temperature factors P2
###############################################################

for i in range(1, 5):
    constrain(u11(i), 4)
    constrain(u22(i), 4)
    constrain(u33(i), 4)

##########################
#  Dataset
#########################
read_data("Ni_2-8.chi.gr", X, 30.0, 0.0)
constrain(qdamp, 21)
setpar(21, 0.0618)
fixpar(21)

constrain(dscale, 26)
setpar(26, 0.8)

pdfrange(1, 1.5, 19.99)



###############################################################
# Run fit
###############################################################

refine()

###############################################################
# Save data
###############################################################

#save_pdf(1, "twophases.pdf")
#save_struct(1, "twophases-a.stru")
#save_struct(2, "twophases-b.stru")
save_res("twophases.res")

