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

/*! Edge class keeps track of the edge information in the mesh
 *
 */
struct Edge
{
  //! Construct a Edge object
  //! \param id edge id
  //! \param i id of 1st vertex
  //! \param j id of 2nd vertex
  Edge(int id, int i, int j) : id(id), i(i), j(j), f1(NO_FACE), f2(NO_FACE), boundary(false)  {   }
  
  //! Check if vertex belongs to the edge
  //! \param v vertex id
  bool vert_in(int v) { return (i == v) || (j == v); }
  
  //! Check if a face belongs to the edge
  bool face_of(int f) { return (f == f1) || (f == f2); }
  
  //! Get the other vertex on the face (-1 if not part of this edge)
  //! \param v vertex index
  int other_vert(int v)
  {
    if (!vert_in(v)) return -1;
    if (v == i) return j;
    return i;
  }
  
  //! Get other edge belonging to the face
  //! \param f face id
  int other_face(int f)
  {
    if (!face_of(f)) return -1;
    if (f == f1) return f2;
    return f1;
  }
  
  int id;                      //!< Edge id
  int i;                       //!< 1st vertex
  int j;                       //!< 2nd vertex

  int f1;                      //!< 1st face (NO_FACE if none)
  int f2;                      //!< 2nd face (NO_FACE if none)
  
  bool boundary;               //!< If true, this is a boundary edge
    
};

ostream& operator<<(ostream&, const Edge&);

#endif