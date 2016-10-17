/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file vertex.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Nov-2015
 * \brief Auxiliary functions for Vertex class
 */ 

#include "vertex.hpp"

ostream& operator<<(ostream& out, const Vertex& v)
{
  out << " ---------- VERTEX ------------------ " << endl;
  out << format("id : %d\n") % v.id
  << format("type : %d\n") % v.type
  << format("(x,y,z) : (%15.9e,%15.9e,%15.9e)\n") % v.r.x % v.r.y % v.r.z
  << format("area : %15.9e\n") % v.area
  << format("perimeter : %15.9e\n") % v.perim
  << format("coordination : %d (%d)\n") % v.n_edges % v.z
  << format("(Nx,Ny,Nz) : (%15.9e,%15.9e,%15.9e)\n") % v.N.x % v.N.y % v.N.z;
  out << "neighbours : ";
  for (int i = 0; i < v.n_edges; i++)
    out << v.neigh[i] << " ";
  out << endl << "edges : ";
  for (int i = 0; i < v.n_edges; i++)
    out << v.edges[i] << " ";
  out << endl << "faces : ";
  for (int i = 0; i < v.n_faces; i++)
    out << v.faces[i] << " ";
  out << endl << "dual : ";
  for (unsigned int i = 0; i < v.dual.size(); i++)
    out << v.dual[i] << " ";
  out << endl << "dual_neighbour_map: ";
  for (unsigned int i = 0; i < v.dual_neighbour_map.size(); i++)
    out << v.dual_neighbour_map[i] << " ";
  out << endl;
  if (v.boundary)
    out << "boundary vertex" << endl;
  if (v.ordered)
    out << "vertex is ordered" << endl;
  if (v.attached)
    out << "vertex is attached" << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}
