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
 * \file edge.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Nov-2015
 * \brief Auxiliary functions for Edge class
 */ 

#include "edge.hpp"

ostream& operator<<(ostream& out, const Edge& e)
{
  out << " ---------- EDGE ------------------ " << endl;
  out << format("id : %d\n") % e.id
  << format("from : %d\n") % e.from
  << format("to : %d\n") % e.to
  << format("face : %d\n") % e.face
  << format("pair : %d\n") % e.pair;
  if (e.boundary)
    out << "boundary edge" << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}

