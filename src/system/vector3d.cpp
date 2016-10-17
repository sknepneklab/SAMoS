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


//! Compute angle between two vectors. 3rd vector is used to get the right sign. 
//! It is normal to the surface at a given point of orgigin of two points. 
//! \param a first vector
//! \param b second vector
//! \param N normal vector
double angle(Vector3d& a, Vector3d& b, const Vector3d& N)
{
  double a_dot_b = dot(a.unit(),b.unit());
  if (a_dot_b > 1.0) a_dot_b = 1.0;
  else if (a_dot_b < -1.0) a_dot_b = -1.0;
  double phi = std::acos(a_dot_b);
  double sign = dot(cross(a,b),N);
  if (sign >= 0) 
    return phi;
  else
    return -phi;
}

/*! Mirror vector with respect to other vector
 * \param P position of the reference point
 * \param n mirror with respect to this vector
 * \param Q mirror this vector
*/
Vector3d mirror(Vector3d& P, Vector3d& n, Vector3d& Q)
{
  Vector3d N = n.unit();
  double n_xx = N.x*N.x, n_xy = N.x*N.y, n_xz = N.x*N.z;
  double n_yx = N.y*N.x, n_yy = N.y*N.y, n_yz = N.y*N.z;
  double n_zx = N.z*N.x, n_zy = N.z*N.y, n_zz = N.z*N.z;
  
  Vector3d P_m_Q = P - Q;
  double qx = Q.x + 2.0*((1.0 - n_xx)*P_m_Q.x         - n_xy*P_m_Q.y -         n_xz*P_m_Q.z);
  double qy = Q.y + 2.0*(      - n_yx*P_m_Q.x + (1.0 - n_yy)*P_m_Q.y -         n_yz*P_m_Q.z); 
  double qz = Q.z + 2.0*(      - n_zx*P_m_Q.x         - n_zy*P_m_Q.y + (1.0 - n_zz)*P_m_Q.z); 
  
  return Vector3d(qx,qy,qz);
}
