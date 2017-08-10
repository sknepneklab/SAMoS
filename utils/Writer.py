# ***************************************************************************
# *
# *  Copyright (C) 2013-2016 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
# *
# *  SAMoS is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  SAMoS is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *****************************************************************************

from Configuration import *
from Tesselation import *
import vtk
from glob import glob
from datetime import *

class Writer:
	def __init__(self,nematic=False,alpha=0,connected=False):
		self.nematic=nematic
		self.connected=connected
		self.alpha=alpha
		
	def writeConfigurationVTK(self,conf,outfile):
		# Data which goes into file: positions, directors, velocities
		# radii
		r = conf.radius
		# types
		tp = conf.ptype
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
		
		Type = vtk.vtkDoubleArray()
		Type.SetNumberOfComponents(1)
		Type.SetName('Type')

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
		for (xx,yy,zz,rr,tt) in zip(x,y,z,r,tp):
			Points.InsertNextPoint(xx,yy,zz)
			Radii.InsertNextValue(rr)
			Type.InsertNextValue(tt)
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
		polydata.GetPointData().AddArray(Type)
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
		#writer.SetDataModeToAscii()
		writer.SetDataModeToBinary()
                writer.SetCompressorTypeToZLib()
		writer.Write()
		
	def writeDefects(self,defects, numdefect,outfile):
		# Preparing the vtp output
		# Create point structure in vtk
		Points = vtk.vtkPoints()
		print "Created Points"
		Charge = vtk.vtkDoubleArray()
		Charge.SetNumberOfComponents(1)
		Charge.SetName('Charge')
		for u in range(numdefect):
			Points.InsertNextPoint(defects[u][1],defects[u][2],defects[u][3])
			Charge.InsertNextValue(defects[u][0])
		
		#lines = vtk.vtkCellArray()
		#line = vtk.vtkLine()
		#for i in range(numdefect_n):
			#line = vtk.vtkLine()
			#line.GetPointIds().SetId(0,0)
			#line.GetPointIds().SetId(1,i+1)
			#lines.InsertNextCell(line)
		#for i in range(numdefect_v):
			#line = vtk.vtkLine()
			#line.GetPointIds().SetId(0,0)
			#line.GetPointIds().SetId(1,numdefect_n+i+1)
			#lines.InsertNextCell(line)
		#print "Added lines"
		
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(Points)
		#polydata.SetLines(lines)
		polydata.GetPointData().AddArray(Charge)
		
		print "Finished Polydata"
		polydata.Modified()
		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outfile)
		# Python 2.7 vs. 3 incompatibility?
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		#writer.SetDataModeToAscii()
		writer.SetDataModeToBinary()
                writer.SetCompressorTypeToZLib()
		writer.Write()
		print "Wrote File"
		
	def writePatches(self,tess,outname,contractile=False):
		print outname
		points = vtk.vtkPoints()
		polygons = vtk.vtkCellArray()
		v=0
		polygon = vtk.vtkPolygon()
		havePoly=[]
		for k in range(len(tess.ParList)):
			nedge=len(tess.ParList[k])
			if nedge<2:
				huh=0
				#print nedge
				#print k
				#print tess.ParList[k]
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
		try:
                        eng, press,ncon,stress = tess.conf.compute_energy_and_pressure()
                except:
                        pass
		contractile = False
		if contractile:
                        print "Are we actually going here??"
			press_c=tess.computeContractile(self.alpha)
			print press_c
			print np.mean(press_c)
			print np.std(press_c)
			print np.min(press_c)
			print np.max(press_c)
			press+=press_c
		#print press
		#print np.mean(press)
		#print np.std(press)
		#print np.min(press)
		#print np.max(press)
		#pressure = vtk.vtkDoubleArray()
		#pressure.SetNumberOfComponents(1)
		#pressure.SetName('Pressure')
		#for k in havePoly:
			#pressure.InsertNextValue(press[k])
		#polygonPolyData.GetCellData().AddArray(pressure)
		
		# Add type
		ncon = vtk.vtkDoubleArray()
		ncon.SetNumberOfComponents(1)
		ncon.SetName('Z')	
		for k in havePoly:
			ncon.InsertNextValue(len(tess.ParList[k]))
		polygonPolyData.GetCellData().AddArray(ncon)
		
		## Add type
		#ptype = vtk.vtkDoubleArray()
		#ptype.SetNumberOfComponents(1)
		#ptype.SetName('Type')	
		#for k in havePoly:
			#ptype.InsertNextValue(tess.conf.ptype[k])
		#polygonPolyData.GetCellData().AddArray(ptype)
		
                # Add denisity
                tess.ComputePatchArea()
            	density = vtk.vtkDoubleArray()
		density.SetNumberOfComponents(1)
		density.SetName('Density')
		for k in havePoly:
                    if tess.conf.geom.manifold=='sphere':
                        N = tess.conf.N
                        R = tess.conf.geom.R
                        A0 = 4.0*np.pi*R**2/N
                        if tess.area[k] < 0.1*A0:
			   density.InsertNextValue(10.0)
                        else:
			   density.InsertNextValue(A0/tess.area[k])
                    else:
			density.InsertNextValue(1.0/tess.area[k])
		polygonPolyData.GetCellData().AddArray(density)


		writer = vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(outname)
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polygonPolyData)
		else:
			writer.SetInputData(polygonPolyData)
		#writer.SetDataModeToAscii()
		writer.SetDataModeToBinary()
                writer.SetCompressorTypeToZLib()
		writer.Write()	
