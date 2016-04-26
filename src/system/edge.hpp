/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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
  Edge(int id, int i, int j) : id(id), from(i), to(j), next(-1), face(NO_FACE), visited(false), pair(-1), boundary(false), dual(-1)  {   }
  
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
    
};

ostream& operator<<(ostream&, const Edge&);

#endif
