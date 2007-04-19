
reset()

read_struct("hj.stru")
#read_struct("hj.init.stru");

read_data("hj.data", N, 29.0, 0.03)

###############################################################
# Experimental and lattice parameters
###############################################################

constrain(lat(1), 1)
constrain(lat(2), 2)
constrain(lat(3), "@3")

setpar(1, lat(1))
setpar(2, lat(2))
setpar(3, lat(3))

constrain(lat(5), 5)
setpar(5, lat(5))

#fixpar(1)
#fixpar(2)
#fixpar(3)
#fixpar(5)

constrain(pscale, 20)
constrain(qdamp, 21)
constrain(delta2, 22)

setpar(20, pscale)
setpar(21, qdamp)
setpar(22, delta2)

setvar(delta1, 0.0)


###############################################################
#position
###############################################################
# group 1-13

for ig in range(1, 14):
    constrain(x(2*ig-1), 100+2*ig-1)
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

for ig in range(14, 24):
    constrain(x(2*ig+25), 100+2*ig-1)
    constrain(x(2*ig+26), 100+2*ig-1)

    constrain(x(2*ig+45), 100+2*ig-1, FCOMP)
    constrain(x(2*ig+46), 100+2*ig-1, FCOMP)

    constrain(z(2*ig+25), 100+2*ig)
    constrain(z(2*ig+26), 100+2*ig)

    constrain(z(2*ig+45), 100+2*ig, FCOMP)
    constrain(z(2*ig+46), 100+2*ig, FCOMP)

    setpar(100+2*ig-1, x(2*ig+25))
    setpar(100+2*ig, z(2*ig+25))

    #fixpar(100+2*ig-1)
    #fixpar(100+2*ig)


###############################################
#Thermal factors
###############################################

# group 1-13

for ig in range(1, 14):

    ip = 200+3*ig-2
    setpar(ip, sqrt(0.003))
    #setpar(ip, sqrt(getvar(u11(2*ig-1))))
    #fixpar(ip)

    for i in range(1, 3):
        constrain(u11(2*ig-2+i), ip, FSQR)
        constrain(u22(2*ig-2+i), ip, FSQR)
        constrain(u33(2*ig-2+i), ip, FSQR)

        constrain(u11(2*ig+24+i), ip, FSQR)
        constrain(u22(2*ig+24+i), ip, FSQR)
        constrain(u33(2*ig+24+i), ip, FSQR)


# group 14-23

for ig in range(14, 24):

    ip = 200+3*ig-2
    setpar(ip, sqrt(0.003))
    #setpar(ip, sqrt(getvar(u11(2*ig+25))))
    #fixpar(ip)

    for i in range(1, 3):
        constrain(u11(2*ig+24+i), ip, FSQR)
        constrain(u22(2*ig+24+i), ip, FSQR)
        constrain(u33(2*ig+24+i), ip, FSQR)

        constrain(u11(2*ig+44+i), ip, FSQR)
        constrain(u22(2*ig+44+i), ip, FSQR)
        constrain(u33(2*ig+44+i), ip, FSQR)

fixpar(ALL)

################
# Refinement
################

pdfrange(ALL, 1.5, 20.0)

freepar(1)
freepar(2)
freepar(3)
freepar(5)
freepar(20)
freepar(21)
freepar(22)

refine()

for ig in range(1, 14):
    freepar(100+2*ig-1)
    freepar(100+2*ig)

refine()

for ig in range(14, 24):
    freepar(100+2*ig-1)
    freepar(100+2*ig)

refine()

for ig in range(1, 14):
    freepar(200+3*ig-2)

refine()

for ig in range(14, 24):
    freepar(200+3*ig-2)

refine()

save_pdf(1, "hj16.pdf")
save_dif(1, "hj16.dif")
save_res("hj16.res")
save_struct(1, "hj.stru")

#calc()
#save_pdf(1, "hj.calc")
