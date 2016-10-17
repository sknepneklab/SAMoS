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
 * \file face.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Nov-2015
 * \brief Auxiliary functions for Face class
 */ 

#include "face.hpp"

ostream& operator<<(ostream& out, const Face& f)
{
  out << " ---------- FACE ------------------ " << endl;
  out << format("id : %d\n") % f.id
  << format("type : %d\n") % f.type
  << format("number of edges : %d\n") % f.n_sides;
  out << "vertices : ";
  for (int i = 0; i < f.n_sides; i++)
    out << f.vertices[i] << " ";
  out << endl << "edges : ";
  for (int i = 0; i < f.n_sides; i++)
    out << f.edges[i] << " ";
  if (!f.is_hole)
  {
    out << endl << "angles : ";
    for (int i = 0; i < f.n_sides; i++)
      out << f.angles[i] << " ";
  }
  out << endl << "area : " << f.area;
  out << endl << "rc : " << f.rc;
  out << endl;
  if (f.is_hole)
    out << "Face is a hole." << endl;
  if (f.boundary)
    out << "Face is at the boundary." << endl;
  if (f.obtuse)
    out << "Face is obtuse." << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}
