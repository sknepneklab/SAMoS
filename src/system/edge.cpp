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
  << format("i : %d\n") % e.i
  << format("j : %d\n") % e.j
  << format("f1 : %d\n") % e.f1
  << format("f2 : %d\n") % e.f2;
  if (e.boundary)
    out << "boundary edge" << endl;
  out << " ------------------------------------";
  out << endl;
  return out;
}

