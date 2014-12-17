# In this script, we will compute the combinatorial structure of
# the Morse-Smale complex after a simplification of 5%. 
#
# We will then save the data and load it back
#
# Finally we will simplify it upto 100% to compute the homology of 
# the underlying domain

#import necessary modules
import pymstri

#create the mscomplex object. Call help(msc) for info. 
msc = pymstri.mscomplex()

# compute the Morse-Smale complex. Scalar function is given in the 4th field of the off file.
msc.compute_off("3wde.off")

# simplify using persistence upto 5% of the function range
msc.simplify_pers(0.05,True,0,0)

# Collect the geometry data
msc.collect_geom()

# print some stats on the cps that survive simplification
print "--------------------"
print "# minima  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 0])
print "# saddles : ", len( [cp for cp in msc.cps() if msc.index(cp) == 1])
print "# maxima  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 2])
print "--------------------"
	
# save the file
msc.save("3wde_mscomplex.bin")

#load back the file
msc.save("3wde_mscomplex.bin")

# print the stats again and check they're the same
print "# minima  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 0])
print "# saddles : ", len( [cp for cp in msc.cps() if msc.index(cp) == 1])
print "# maxima  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 2])
print "--------------------"

# simplify using persistence upto 100%
msc.simplify_pers(1.0,True,0,0)

# print the homology
print "# components  : ", len( [cp for cp in msc.cps() if msc.index(cp) == 0])
print "# tunnels     : ", len( [cp for cp in msc.cps() if msc.index(cp) == 1])
print "# voids       : ", len( [cp for cp in msc.cps() if msc.index(cp) == 2])
print "--------------------"
