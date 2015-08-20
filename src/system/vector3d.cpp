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
 * \file vector3d.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Aug-2015
 * \brief Definitions of Vector3d class auxiliary functions
 */ 

#include "vector3d.hpp"

//! Output vector components
std::ostream& operator<<(std::ostream& os, const Vector3d& v)
{
  os << "(" << v.x << "," << v.y << "," << v.z << ")";
  return os;
}

//! Compute dot product between two vectors
double dot(const Vector3d& v1, const Vector3d& v2)
{
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

//! Compute cross product between two vectors
Vector3d cross(const Vector3d& v1, const Vector3d& v2)
{
  return Vector3d(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

//! Scale vector by a number
Vector3d operator*(const double c, const Vector3d& v)
{
  return Vector3d(c*v.x, c*v.y, c*v.z);
}