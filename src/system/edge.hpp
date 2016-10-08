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
 * \file edge.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Declaration of Edge class.
 */ 

#ifndef __EDGE_HPP__
#define __EDGE_HPP__

const int NO_FACE = -1;   // value indicting that one of the faces is not set

#include <iostream>
#include <boost/format.hpp>

using std::ostream;
using std::endl;
using boost::format;

/*! Edge class keeps track of the edge information in the mesh.
 *  Edges are directed, so this data structure implement a half-edge.
 */
struct Edge
{
  //! Construct a Edge object
  //! \param id edge id
  //! \param i id of 1st vertex
  //! \param j id of 2nd vertex
  Edge(int id, int i, int j) : id(id), 
                               from(i), 
                               to(j),
                               next(-1), 
                               face(NO_FACE), 
                               visited(false), 
                               pair(-1), 
                               boundary(false), 
                               dual(-1),
                               attempted_removal(false)
                               {   }
  
  //! Check if vertex belongs to the edge
  //! \param v vertex id
  bool vert_in(int v) { return (from == v) || (to == v); }
  
  //! Get the other vertex on the face (-1 if not part of this edge)
  //! \param v vertex index
  int other_vert(int v)
  {
    if (!vert_in(v)) return -1;
    return to;
  }
  
  //! Compare if the id of the edge is equal to an integer value
  //! \val edge index
  bool operator==(int val)
  {
     return (id == val);
  }
    
  int id;                       //!< Edge id
  int from;                     //!< Id of vertex half-edge starts at
  int to;                       //!< Id of vertex half-edge points to
  int next;                     //!< Index of the next edge in the face 

  int face;                     //!< Id of the face (to the left) this half-edge belongs to. NO_FACE if none, i.e. boundary edge.
  bool visited;                 //!< If true, edge visited while building faces.
  int pair;                     //!< Id of the other half of the half-edge pair
  bool boundary;                //!< Edge is boundary if its face is a hole
  int dual;                     //!< Id of the dual vertex (-1 if none)
  bool attempted_removal;       //!< If true, edge removal has been attempted
    
};

ostream& operator<<(ostream&, const Edge&);

#endif
