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
  out << endl << "angles : ";
  for (int i = 0; i < f.n_sides; i++)
    out << f.angles[i] << " ";
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
