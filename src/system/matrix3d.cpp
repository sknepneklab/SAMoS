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
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file vector3d.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Mar-2016
 * \brief Definitions of Matrix3d class auxiliary functions
 */ 

#include "matrix3d.hpp"

//! Output vector components
std::ostream& operator<<(std::ostream& os, const Matrix3d& m)
{
  os << "[";
  os << "[" << m.M[0][0] << "," << m.M[0][1] << "," << m.M[0][2] << "],";
  os << "[" << m.M[1][0] << "," << m.M[1][1] << "," << m.M[1][2] << "],";
  os << "[" << m.M[2][0] << "," << m.M[2][1] << "," << m.M[2][2] << "]";
  os << "]";
  return os;
}