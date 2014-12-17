#import necessary modules
import pymstri

#create the mscomplex object. Call help(msc) for info. 
msc = pymstri.mscomplex()

# compute the Morse-Smale complex. Scalar function is given in the 4th field of the off file.
msc.compute_off("3wde.off")

# simplify using persistence upto 100%
msc.simplify_pers(1.0,True,0,0)

# print the homology
print "# components  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 0])
print "# tunnels     : ", len( [cp for cp in msc.cps() if msc.index(cp) == 1])
print "# voids       : ", len( [cp for cp in msc.cps() if msc.index(cp) == 2])
