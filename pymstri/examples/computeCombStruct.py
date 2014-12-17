# In this script, we will compute the combinatorial structure of
# the Morse-Smale complex after a simplification of 5%. 
#
# We will then save the data to a simple ascii file. 

#import necessary modules
import pymstri

#create the mscomplex object. Call help(msc) for info. 
msc = pymstri.mscomplex()

# compute the Morse-Smale complex. Scalar function is given in the 4th field of the off file.
msc.compute_off("3wde.off")

# simplify using persistence upto 5% of the function range
msc.simplify_pers(0.05,True,0,0)

# get the object representing the triangulation. Call help(tcc) for info.
tcc = msc.get_tri_cc()

# Write the combinatorial mscomplex data to a file  
fo = open("3wde_combStruct.txt","w")

# First, write the number of critical points
fo.write("#numcps\n")
fo.write(str(len(msc.cps()))+"\n")

# Write the per cp information. Note: Each cp is identified by a unique id.
fo.write("#    id       x         y         z         fn        index\n")
for cp in msc.cps():
	x,y,z = tcc.point(msc.vertid(cp))
	fn,index = msc.fn(cp),msc.index(cp)
	fo.write("%5d   %8.4f  %8.4f  %8.4f  %8.6f  %5d\n"%(cp,x,y,z,fn,index))

# Write the combinatorial adjaceny list of each cp
fo.write("#  id  NumAdj adjCp1 adjCp2 ..\n")
for cp in msc.cps():
	conn = list(msc.asc(cp)) + msc.des(cp)
	
	ln  = ""
	ln += str(cp).ljust(10)
	ln += (str(len(conn))).ljust(10)
	ln += "".join(str(i).ljust(10) for i in conn)
	ln += "\n"

	fo.write(ln)

#finished
fo.write("# The end")


