
#####################################################################
# Ni- (FCC) Refinement - average structure
######################################################################

read_struct("ni.stru")

show_scat(X)

set_scat(X,2,1,2,3,4,5,6,7,8,9)

reset_scat(X,2)

#

show_scat(N)

set_scat(N,1,1.2)

reset_scat(N,1)

#blen(1,1,1,3)

#read_struct("Ni.stru")

#read_data("Ni_2-8.chi.gr",X,30.0,0.0)

###############################################################
# Experimental and lattice parameters
###############################################################

#constrain(lat(1),"log(exp(@1))")
#constrain(lat(1),"-log(@23)-1-@22-4.1*sqr(exp(3.5*@23+2))-cos(@22*@23+@22*@23*@23-@22)")
#constrain(lat(1),"sqr(@23)-1-@22-4.1*sqr(exp(3.5*@23+2))-cos(@22+@23+@22+@23+@23-@22)")
#constrain(lat(1),"sqr(@23)-sqr(exp(3.5*@23+2))")
#constrain(lat(1),"-++@23")
#constrain(lat(1),"-@22")
#constrain(lat(1),"    -sqr(@22)")
#constrain(lat(1),"<1>")
#constrain(lat(1),1)
#constrain(lat(2),"@22-4.1*(3.5*@23+2)")

#setpar(1,7.186457)
#setpar(22,0.42)

#setpar(23,0.25)
#fixpar(22)
#calc() 

#sys.exit(0)