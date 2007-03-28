# FIXME: this file is obsolete and needs to be cleaned up or removed
# PJ, 2007-03-27

from pdffit import *

reset()

read_struct("hj.stru")
#read_struct("hj.init.stru");

read_data("hj.data", N, 29.0, 0.03)

###############################################################
# Experimental and lattice parameters
###############################################################

constrain(lat(1),"@1")
constrain(lat(2),"@2")
constrain(lat(3),"@3")

setpar(1,lat(1))
setpar(2,lat(2))
setpar(3,lat(3))

constrain(lat(5),"tand(atand(@5))")
setpar(5,lat(5))

#fixpar(1)
#fixpar(2)
#fixpar(3)
#fixpar(5)

constrain(pfrac,"cosd(acosd(@20))")
constrain(qsig,"atand(tand(@21))")
constrain(delta2,"exp(log(@22))")
constrain(delta1,"log10(exp10(@23))")

setpar(20,pfrac)
setpar(21,0.017)
#setpar(22,delta2)

#setpar(20,1.2)
#setpar(21,0.03)
setpar(22,delta2)

#setpar(23,delta1)
setpar(23,0.0)
#fixpar(23)

#constrain(srat,24)
#setpar(24,srat)
#setvar(rcut,3.75)

###############################################################
#position
###############################################################
# group 1-13

#for ig in range(1, 14):
ig = 1
constrain(x(2*ig-1), "sqr(sqrt(@101))")
constrain(x(2*ig), 100+2*ig-1)

constrain(x(2*ig+25), 100+2*ig-1, FCOMP)
constrain(x(2*ig+26), 100+2*ig-1, FCOMP)

constrain(z(2*ig-1), 100+2*ig)
constrain(z(2*ig), 100+2*ig)

constrain(z(2*ig+25), 100+2*ig, FCOMP)
constrain(z(2*ig+26), 100+2*ig, FCOMP)

setpar(100+2*ig-1, x(2*ig-1))
setpar(100+2*ig, z(2*ig-1))

#fixpar(100+2*ig-1)
#fixpar(100+2*ig)

# group 14-23

#for ig in range(14, 24):
ig = 14
constrain(y(2*ig+25), 100+2*ig-1)

constrain(occ(2*ig+25), 100+2*ig)

setpar(100+2*ig-1, y(2*ig+25))
setpar(100+2*ig, occ(2*ig+25))

#fixpar(100+2*ig-1)
#fixpar(100+2*ig)


###############################################
#Thermal factors
###############################################

# group 1-13

#for ig in range(1, 14):

ig = 1
ip = 200+3*ig-2
#setpar(ip, sqrt(0.003))
setpar(ip, sqrt(getvar(u11(2*ig-1))))
#fixpar(ip)

for i in range(1, 3):
    constrain(u11(2*ig-2+i), ip, FSQR)
    constrain(u22(2*ig-2+i), ip, FSQR)
    constrain(u33(2*ig-2+i), ip, FSQR)

    constrain(u11(2*ig+24+i), ip, FSQR)
    constrain(u22(2*ig+24+i), ip, FSQR)
    constrain(u33(2*ig+24+i), ip, FSQR)


# group 14-23

#for ig in range(14, 24):

ig = 14
ip = 200+3*ig-2
#setpar(ip, sqrt(0.003))
setpar(ip, sqrt(0.000001))
#fixpar(ip)

for i in range(1, 3):
    constrain(u12(2*ig+24+i), ip)
    constrain(u13(2*ig+24+i), ip)
    constrain(u23(2*ig+24+i), ip)


#fixpar(ALL)

################
# Refinement
################

pdfrange(ALL,1.5,20.0)

refine()
save_pdf(1,"hj15.pdf")
save_dif(1,"hj15.dif")
save_res("hj15.res")
save_struct(1,"hj.stru")

#calc()
#save_pdf(1,"hj.calc")
