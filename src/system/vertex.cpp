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
  << format("coordination : %d\n") % v.n_edges;
  if (v.boundary)
    out << "boundary vertex" << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}