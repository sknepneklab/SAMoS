# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#   Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#
#    Author: Rastko Sknepnek
#   
#    Division of Physics
#    School of Engineering, Physics and Mathematics
#    University of Dundee
#    
#    (c) 2013, 2014, 2015
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

from Configuration import *
from Tesselation import *
import vtk
from glob import glob
from datetime import *

class Writer:
	def __init__(self,nematic=False,connected=False):
		self.nematic=nematic
		self.connected=connected
		
	def writeConfigurationVTK(self,conf,outfile):
		# Data which goes into file: positions, directors, velocities
		# radii
		r = conf.radius
		# positions
		x = conf.rval[:,0]
		y = conf.rval[:,1]
		z = conf.rval[:,2]
		# directors
		nx = conf.nval[:,0]
		ny = conf.nval[:,1]
		nz = conf.nval[:,2]
		# velocities
		vx = conf.vval[:,0]
		vy = conf.vval[:,1]
		vz = conf.vval[:,2]
		
		# Preparking the vtk structure
		Points = vtk.vtkPoints()
		
		Radii = vtk.vtkDoubleArray()
		Radii.SetNumberOfComponents(1)
		Radii.SetName('Radius')

		Velocities = vtk.vtkDoubleArray()
		Velocities.SetNumberOfComponents(3)
		Velocities.SetName("Velocity")

		Directors = vtk.vtkDoubleArray()
		Directors.SetNumberOfComponents(3)
		Directors.SetName("Directors")
		
		if self.nematic:
			NDirectors = vtk.vtkDoubleArray()
			NDirectors.SetNumberOfComponents(3)
			NDirectors.SetName("NDirectors")
		
		# Adding the data to the vtk structures
		for (xx,yy,zz,rr) in zip(x,y,z,r):
			Points.InsertNextPoint(xx,yy,zz)
			Radii.InsertNextValue(rr)
		for (vvx,vvy,vvz) in zip(vx,vy,vz):
			Velocities.InsertNextTuple3(vvx,vvy,vvz)
		for (nnx,nny,nnz) in zip(nx,ny,nz):	
			if self.nematic:
				Directors.InsertNextTuple3(0.5*nnx,0.5*nny,0.5*nnz)
				NDirectors.InsertNextTuple3(-0.5*nnx,-0.5*nny,-0.5*nnz)
			else:
				Directors.InsertNextTuple3(nnx,nny,nnz)
		# Connected, using convex hull (? ask Rastko ...?)
		if self.connected:
			Lines = vtk.vtkCellArray()
			Line = vtk.vtkLine()
			points = np.column_stack((x,y,z)) 
			hull = ConvexHull(points)
			edges = []
			for h in hull.simplices:
				i, j, k = h
				if not sorted([i,j]) in edges: edges.append(sorted([i,j]))
				if not sorted([i,k]) in edges: edges.append(sorted([i,k]))
				if not sorted([j,k]) in edges: edges.append(sorted([j,k]))
				for (i,j) in edges:
					Line.GetPointIds().SetId(0,i)
					Line.GetPointIds().SetId(1,j)
					Lines.InsertNextCell(Line)
			
		# Putting the results into a polydata structure
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(Points)
		if self.connected:
			polydata.SetLines(Lines)
		polydata.GetPointData().AddArray(Radii)
		polydata.GetPointData().AddArray(Velocities)
		polydata.GetPointData().AddArray(Directors)
		if self.nematic:
			polydata.GetPointData().AddArray(NDirectors)
		polydata.Modified()
		
		# Finally, output via binary writer
		writer = vtk.vtkXMLPolyDataWriter()
		#outname = '.'.join(f.split('.')[:-1])
		writer.SetFileName(outfile)
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		writer.SetDataModeToAscii()
		writer.Write()
		
	def writeDefects(self,defects_n, defects_v,numdefect_n,numdefect_v,outfile):
		# Preparing the vtp output
		# Create point structure in vtk
		Points = vtk.vtkPoints()
		print "Created Points"
		# Create (something) associated to the points, with different values for each
		Number = vtk.vtkDoubleArray()
		Number.SetNumberOfComponents(1)
		Number.SetName('Number')
		Size = vtk.vtkDoubleArray()
		Size.SetNumberOfComponents(1)
		Size.SetName('Size')
		print "Created Number"
		# Put one point at the centre, and the ndefect ones around it
		Points.InsertNextPoint(0,0,0)
		Number.InsertNextValue(0)
		Size.InsertNextValue(0)
		for u in range(numdefect_n):
			Points.InsertNextPoint(defects_n[u,1],defects_n[u,2],defects_n[u,3])
			Number.InsertNextValue(1)
			Size.InsertNextValue(1.0)
		for u in range(numdefect_v):
			Points.InsertNextPoint(defects_v[u,1],defects_v[u,2],defects_v[u,3])
			Number.InsertNextValue(2)
			Size.InsertNextValue(1.0)
		print "Added Particles and Numbers"
		
		lines = vtk.vtkCellArray()
		line = vtk.vtkLine()
		for i in range(numdefect_n):
			line = vtk.vtkLine()
			line.GetPointIds().SetId(0,0)
			line.GetPointIds().SetId(1,i+1)
			lines.InsertNextCell(line)
		for i in range(numdefect_v):
			line = vtk.vtkLine()
			line.GetPointIds().SetId(0,0)
			line.GetPointIds().SetId(1,numdefect_n+i+1)
			lines.InsertNextCell(line)
		print "Added lines"
		
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(Points)
		polydata.SetLines(lines)
		polydata.GetPointData().AddArray(Number)
		polydata.GetPointData().AddArray(Size)
		print "Finished Polydata"
		polydata.Modified()
		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outfile)
		# Python 2.7 vs. 3 incompatibility?
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		writer.SetDataModeToAscii()
		writer.Write()
		print "Wrote File"
		
	def writePatches(self,tess,outname):
		print outname
		points = vtk.vtkPoints()
		polygons = vtk.vtkCellArray()
		v=0
		polygon = vtk.vtkPolygon()
		havePoly=[]
		for k in range(len(tess.ParList)):
			nedge=len(tess.ParList[k])
			if nedge<2:
				print nedge
				print k
				print tess.ParList[k]
			else:
				havePoly.append(k)
				#for k in range(300):
				# Create the points of the polygon: the loop centers
				polygon = vtk.vtkPolygon()
				for l in tess.ParList[k]:
					if tess.geom.periodic:
						dl=tess.geom.ApplyPeriodic11(tess.rval[k,:],tess.LoopCen[l])
						dl+=tess.rval[k,:]
						points.InsertNextPoint(dl[0],dl[1],dl[2])
					else:
						points.InsertNextPoint(tess.LoopCen[l][0],tess.LoopCen[l][1],tess.LoopCen[l][2])
				polygon.GetPointIds().SetNumberOfIds(nedge)
				for l in range(nedge):
					#print l
					polygon.GetPointIds().SetId(l,v+l)
				
				polygons.InsertNextCell(polygon)
				v+=nedge
		# Create the matching polydata 
		polygonPolyData = vtk.vtkPolyData()
		polygonPolyData.SetPoints(points)
		polygonPolyData.SetPolys(polygons)
		# Add stresses ...
		eng, press,stress = tess.conf.compute_energy_and_pressure()
		pressure = vtk.vtkDoubleArray()
		pressure.SetNumberOfComponents(1)
		pressure.SetName('Pressure')
		for k in havePoly:
			pressure.InsertNextValue(press[k])
		polygonPolyData.GetCellData().AddArray(pressure)
		
		# Add type
		ptype = vtk.vtkDoubleArray()
		ptype.SetNumberOfComponents(1)
		ptype.SetName('Type')	
		for k in havePoly:
			ptype.InsertNextValue(tess.conf.ptype[k])
		polygonPolyData.GetCellData().AddArray(ptype)
			
		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outname)
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polygonPolyData)
		else:
			writer.SetInputData(polygonPolyData)
		writer.SetDataModeToAscii()
		writer.Write()	