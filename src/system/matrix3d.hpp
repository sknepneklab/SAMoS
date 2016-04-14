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
 * \file matrix3d.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Mar-2016
 * \brief Declaration of Matrix3d class.
 */ 

#ifndef __MATRIX3D_HPP__
#define __MATRIX3D_HPP__

#include <cmath>
#include <iostream>

#include "vector3d.hpp"

/*! Matrix3d class
 *  Handles 3x3 real matrices for vector operations
 */
class Matrix3d
{
public:
  
  //! Deault constructor
  Matrix3d()
  { 
    M[0][0] = 0.0;   M[0][1] = 0.0;   M[0][2] = 0.0;
    M[1][0] = 0.0;   M[1][1] = 0.0;   M[1][2] = 0.0;
    M[2][0] = 0.0;   M[2][1] = 0.0;   M[2][2] = 0.0;  
  }
  
  double M[3][3];              //!< Holds the matrix
};

//! Compute product between matrix and a vector
inline Vector3d operator*(const Matrix3d& m, const Vector3d& v)
{
  double X = m.M[0][0]*v.x + m.M[0][1]*v.y + m.M[0][2]*v.z; 
  double Y = m.M[1][0]*v.x + m.M[1][1]*v.y + m.M[1][2]*v.z;
  double Z = m.M[2][0]*v.x + m.M[2][1]*v.y + m.M[2][2]*v.z;
  return Vector3d(X, Y, Z);
}

//! Compute product between vector and a matrix
inline Vector3d operator*(const Vector3d& v, const Matrix3d& m)
{
  double X = m.M[0][0]*v.x + m.M[1][0]*v.y + m.M[2][0]*v.z; 
  double Y = m.M[0][1]*v.x + m.M[1][1]*v.y + m.M[2][1]*v.z;
  double Z = m.M[0][2]*v.x + m.M[1][2]*v.y + m.M[2][2]*v.z;
  return Vector3d(X, Y, Z);
}

ostream& operator<<(ostream&, const Matrix3d&);

#endif
