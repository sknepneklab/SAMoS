/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
