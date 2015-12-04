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

//! Compute angle between two vectors. 3rd vector is used to get the right sign. 
//! It is normal to the surface at a given point of orgigin of two points. 
//! \param a first vector
//! \param b second vector
//! \param N normal vector
double angle(Vector3d& a, Vector3d& b, const Vector3d& N)
{
  double phi = std::acos(dot(a.unit(),b.unit()));
  double sign = dot(cross(a,b),N);
  if (sign >= 0) 
    return phi;
  else 
    return -phi + 2.0*M_PI;
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