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
 * \file vector3d.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Aug-2015
 * \brief Declaration of Vector3d class.
 */ 

#ifndef __VECTOR3D_HPP__
#define __VECTOR3D_HPP__

#include <cmath>
#include <iostream>

#include "particle.hpp"

/*! Vector3d class
 *  Handles vectors in 3d Eucledian space
 */
class Vector3d
{
public:
  
  //! Deault constructor
  Vector3d() : x(0.0), y(0.0), z(0.0) { }
  //! Constructor for a vector object
  Vector3d(double x, double y, double z) : x(x), y(y), z(z) { }
  //! Copy constructor  
  Vector3d(const Vector3d& v) { x = v.x; y = v.y; z = v.z; }
  //! Assignment operator
  Vector3d& operator=(const Vector3d& rhs)
  {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return *this;
  }
  
  //! Add two vectors
  Vector3d operator+(const Vector3d& v)
  {
    return Vector3d(x + v.x, y + v.y, z + v.z);
  }
  
  //! Subtract two vectors
  Vector3d operator-(const Vector3d& v)
  {
    return Vector3d(x - v.x, y - v.y, z - v.z);
  }
  
  //! Subtract two vectors (constant version)
  Vector3d operator-(const Vector3d& v) const
  {
    return Vector3d(x - v.x, y - v.y, z - v.z);
  }
  
  //! Negate a vector
  Vector3d operator-()
  {
    return Vector3d(-x, -y, -z);
  }
  
  //! Scale vector by a constant
  Vector3d operator*(const double c)
  {
    return Vector3d(c*x, c*y, c*z);
  }
  
  //! Test equality
  bool operator==(const Vector3d& v)
  {
    return (x == v.x && y == v.y && z == v.z);
  }
  
  //! Add vector to current vector
  Vector3d& operator+=(const Vector3d& v)
  {
    x += v.x;  y += v.y;  z += v.z;
    return *this;
  }
  
  //! Subtract vector from current vector
  Vector3d& operator-=(const Vector3d& v)
  {
    x -= v.x;  y -= v.y;  z -= v.z;
    return *this;
  }
  
  //! Euclidean dot product with another vector
  double dot(const Vector3d& v)
  {
    return x*v.x + y*v.y + z*v.z;
  }
  
  //! Cross prduct with another vector
  Vector3d cross(const Vector3d& v)
  {
    return Vector3d(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
  }
  
  //! Vector length 
  double len() { return std::sqrt(x*x + y*y + z*z); }
  
  //! Vector length squared
  double len2() { return x*x + y*y + z*z; }
  
  //! Rescale vactor
  void scale(double s) { x *= s; y *= s;  z *= s;  }
  
  //! Return rescale vactor
  Vector3d scaled(double s) { return Vector3d(s*x,s*y,s*z);  }
  
  //! Make the vector has unit lenght
  void normalize()
  {
    double len = this->len();
    if (len != double(0))
      this->scale(double(1)/len);
  }
  
  //! Return unit vector in the direction of this vector
  Vector3d unit()
  {
    double len = this->len();
    if (len != double(0))
      return Vector3d(x/len,y/len,z/len);
      return Vector3d(x,y,z);
  }
  
  //! Get the part of the vector parallel to a given vector
  Vector3d parallel_projection(const Vector3d& v)
  {
    Vector3d n = v;
    n.normalize();
    double proj = this->dot(n);
    return Vector3d(proj*n.x, proj*n.y, proj*n.z);
  }
  
  //! Get the part normal to a given vector
  Vector3d perp_projection(const Vector3d& v)
  {
    Vector3d n = this->parallel_projection(v);
    return Vector3d(x-n.x, y-n.y, z-n.z);
  }
  
  //! Rotate vector by \f$ \phi \f$ around another vector
  Vector3d rotate(const double phi, const Vector3d& v)
  {
    double s = std::sin(phi);
    double c = std::cos(phi);
    double k = 1.0 - c;
    
    double nx = x * (c + k * v.x * v.x) + y * (k * v.x * v.y - s * v.z) + z * (k * v.x * v.z + s * v.y);
    double ny = x * (k * v.x * v.y + s * v.z) + y * (c + k * v.y * v.y) + z * (k * v.y * v.z - s * v.x);
    double nz = x * (k * v.x * v.z - s * v.y) + y * (k * v.y * v.z + s * v.x) + z * (c + k * v.z * v.z);
    
    return Vector3d(nx,ny,nz);
  }

  ///@{
  double x, y, z;              //!< Position in the embedding 3d flat space
  //@}

};

//! Compute cross product between two vectors
inline Vector3d cross(const Vector3d& v1, const Vector3d& v2)
{
  return Vector3d(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

//! Scale vector by a number
inline Vector3d operator*(const double c, const Vector3d& v)
{
  return Vector3d(c*v.x, c*v.y, c*v.z);
}

ostream& operator<<(ostream&, const Vector3d&);
double dot(const Vector3d&, const Vector3d&);
//Vector3d cross(const Vector3d&, const Vector3d&);
//Vector3d operator*(const double, const Vector3d&);
double angle(Vector3d&, Vector3d&, const Vector3d&);
Vector3d mirror(Vector3d&, Vector3d&, Vector3d&);

#endif
