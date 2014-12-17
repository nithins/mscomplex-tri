# In this script, we will compute a segmentation of a molecular surface.
#
# Significant protrusions are identified and segmented. 
#
# The segments are then color tagged in a surface file that is 
# visualizable in vtk/paraview. 

#import necessary modules
import pymstri,vtk,random

#create the mscomplex object. 
msc = pymstri.mscomplex()

# compute the Morse-Smale complex. Scalar function is given in the 4th field of the off file.
msc.compute_off("3wde.off")

# simplify using persistence upto 5% of the function range
msc.simplify_pers(0.05,True,0,0)

# collect the msc geometry after simplification
msc.collect_geom()

# get the object representing the triangulation. 
tcc = msc.get_tri_cc()

#get number of cells of each type
nv,ne,nt = tcc.num_dcells(0),tcc.num_dcells(1),tcc.num_dcells(2)

#create and insert the points
pts = vtk.vtkPoints()
for cellid in range(0,nv):
	pts.InsertNextPoint(tcc.point(cellid))
	
#create a vtk_poly_data
pd = vtk.vtkPolyData()
pd.Allocate()

#set its pts
pd.SetPoints(pts)

#iterate over triangle cells  and insert the polys
tcell = vtk.vtkIdList()
tcell.SetNumberOfIds(3)

for cellid in range(nv+ne,nv+ne+nt):
	i,j,k = tcc.cell_vert(cellid)
	tcell.SetId(0,i)
	tcell.SetId(1,j)
	tcell.SetId(2,k)	
	pd.InsertNextCell(vtk.VTK_TRIANGLE,tcell)
	

#create a color array
colors = vtk.vtkDoubleArray()
colors.SetNumberOfComponents(3)
colors.SetNumberOfTuples(nt)
colors.SetName("Color")

# get a list of surviving maxima
survMax = [cp for cp in msc.cps() if msc.index(cp) == 2]

# iterate over the des manifold of each survigin max. 
for cp in survMax:
	#pick a random color
	color = random.random(),random.random(),random.random()

	# The des mfold comprises of a set of triangles. 	
	for cellid in msc.des_geom(cp):
		#triangle cellids are offset by nv+ne in tcc
		tid = cellid - nv - ne 
		colors.SetTuple(tid,color)

#Add the data as cell data
pd.GetCellData().AddArray(colors)


#write the data to a vtk file
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName( "3wde_maxSeg.vtp" )
writer.SetInput( pd)
writer.Write()







