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
 * \file face.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Declaration of Face class.
 */ 

#ifndef __FACE_HPP__
#define __FACE_HPP__

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/format.hpp>

#include "vector3d.hpp"

using boost::format;
using std::ostream;
using std::string;
using std::vector;
using std::endl;
using std::find;

/*! Face class keeps track of the face information in the mesh
 *
 */
struct Face
{
  //! Construct a Face object
  //! \param id face id
  Face(int id) : id(id), n_sides(0), edge_face(false), ordered(false), type(1), rc(0,0,0) {   }
  
  //! Add a vertex
  //! \param v adds a vertex
  void add_vertex(int v)
  {
    vertices.push_back(v);
    n_sides++;
  }
  
  //! Add an edge
  //! \param e edge index
  void add_edge(int e)
  {
    edges.push_back(e);
  }
  
  //! Check if vertex belongs to the face
  //! \param v vertex id
  bool has_vertex(int v)
  {
    return (find(vertices.begin(), vertices.end(), v) != vertices.end());
  }
  
  
  int id;                      //!< Face id
  int n_sides;                 //!< Number of sides
  bool edge_face;              //!< Face is an edge face if all its edges are at the boundary
  bool ordered;                //!< if true, vertices in the face are ordered
  int type;                    //!< Face type. This is help determine parameters for interactions that depend in the dual vertex
  
  Vector3d rc;                 //!< Coordiantes of geometric centre of the face
  
  vector<int> vertices;        //!< Contains all vertices
  vector<int> edges;           //!< contains all edges
    
};

ostream& operator<<(ostream&, const Face&);

#endif