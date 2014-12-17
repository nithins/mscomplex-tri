# In this script, we will extract the an extremum graph from the 
# combinatorial Morse-Smale complex.
#
# The extremum graph will consist of only maxima and 1-saddles. 
#
# The extrema that are retained are those that have persistence above 5% 
#
# We will also generate a persistence hierarchy to figure out 
# persistence pairs. 
#
# We will retain only those 1-saddles a persistence pair with a maxima
# and infinite persistent 1-saddles. 

#import necessary modules
import pymstri

#create the mscomplex object. Call help(msc) for info. 
msc = pymstri.mscomplex()

# compute the Morse-Smale complex. Scalar function is given in the 4th field of the off file.
msc.compute_off("3wde.off")

# simplify using persistence upto 5% of the function range
msc.simplify_pers(0.05,True,0,0)
  
# Generate the persistence hierarchy
msc.gen_pers_hierarchy()

# get the set of persistence pairs and form a dictionary for lookups
ppairs     = [msc.canc(i) for i in range(msc.num_canc())]
ppairsDict = dict(ppairs + [(b,a) for a,b in ppairs])

# get the trianulation object for vertex information
tcc        = msc.get_tri_cc()

# compute a list of cps that fulfill the above criteria
cpset = [ cp for cp in msc.cps() 
		if msc.index(cp) == 2 or 
		(msc.index(cp) == 1 and 
			(cp not in ppairsDict or 
			msc.index(ppairsDict[cp]) == 2) )]

# Write the extremum graph data to a file  
fo = open("3wde_extremumGraph.txt","w")

# First, write the number of critical points
fo.write("#numcps\n")
fo.write(str(len(cpset))+"\n")

# Write the per cp information. Note: Each cp is identified by a unique id.
fo.write("#    id       x         y         z         fn        index\n")
for cp in msc.cps():
	x,y,z = tcc.point(msc.vertid(cp))
	fn,index = msc.fn(cp),msc.index(cp)
	fo.write("%5d   %8.4f  %8.4f  %8.4f  %8.6f  %5d\n"%(cp,x,y,z,fn,index))

# Write the combinatorial adjaceny list of each cp
fo.write("#  id  NumAdj adjCp1 adjCp2 ..\n")
for cp in msc.cps():
	# filter connections to retain only relevent connections
	conn = msc.asc(cp) if msc.index(cp) == 1 else msc.des(cp)    
	conn = [i for i in conn if i in cpset]
	
	ln  = ""
	ln += str(cp).ljust(10)
	ln += (str(len(conn))).ljust(10)
	ln += "".join(str(i).ljust(10) for i in conn)
	ln += "\n"

	fo.write(ln)

#finished
fo.write("# The end")
