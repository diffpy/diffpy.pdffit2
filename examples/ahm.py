
reset()

#####################################################################
# Ni- (FCC) Refinement - average structure
######################################################################

read_struct("Ni.stru")

read_data("Ni_2-8.chi.gr", X, 30.0, 0.0)

###############################################################
# Experimental and lattice parameters
###############################################################

constrain(lat(1), "@1")
constrain(lat(2), "@1")
#constrain(lat(2), 1)
constrain(lat(3), "@1")

setpar(1, lat(1))

constrain(pscale, 20)
constrain(qdamp, 21)
constrain(delta2, 22)

setpar(20, 0.989)
setpar(21, 0.03)
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

###############################################################
# Run fit
###############################################################

pdfrange(1, 1.5, 19.99)

ok = refine()

###############################################################
# Save data
###############################################################

if ok:
    save_pdf(1, "ahm5.pdf")
    save_dif(1, "ahm5.dif")
    save_struct(1, "ahm5.rstr")
    save_res("ahm5.res")

