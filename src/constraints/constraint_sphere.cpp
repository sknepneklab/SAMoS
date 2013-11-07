/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_sphere.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Implementation of the spherical constraint
 */ 

#include "constraint_sphere.hpp"

/*! Force particle to be confined to the surface of the sphere and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down using simple radial projection.
 *  \param p Particle which is to be projected onto the sphere
 */
void ConstraintSphere::enforce(Particle& p)
{
  double x = p.x, y = p.y, z = p.z;
  double R = sqrt(x*x + y*y + z*z);
  double s = m_r/R;
  // Scale back to the surface
  p.x *= s; p.y *= s; p.z *= s;
  // Compute unit normal
  double nx = p.x/m_r, ny = p.y/m_r, nz = p.z/m_r;
  // Compute tangent component of the velocity
  p.vx -= p.vx*nx; p.vy -= p.vy*ny; p.vz -= p.vz*nz;
}

/*! Rotate velocity vector of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto sphere and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintSphere::rotate_velocity(Particle& p, double phi)
{
  double x = p.x, y = p.y, z = p.z;
  double R = sqrt(x*x + y*y + z*z);
  // Compute unit normal
  double nx = x/R, ny = y/R, nz = z/R;
  // Compute angle sins and cosines
  double c = cos(phi), s = sin(phi);
  // Assemble rotation matrix  
  double Rxx = c + nx*nx*(1.0-c),    Rxy = nx*ny*(1.0-c) - nz*s,  Rxz = nx*nz*(1.0-c) + ny*s;
  double Ryx = ny*nx*(1.0-c) + nz*s, Ryy = c + ny*ny*(1.0-c),     Ryz = ny*nz*(1.0-c) - nx*s;
  double Rzx = nz*nx*(1.0-c) - ny*s, Rzy = nz*ny*(1.0-c) - nx*s,  Rzz = c + ny*ny*(1.0-c);
  // Apply rotation matrix to velocity 
  double vx = Rxx*p.vx + Rxy*p.vy + Rxz*p.vz;
  double vy = Ryx*p.vx + Ryy*p.vy + Ryz*p.vz;
  double vz = Rzx*p.vx + Rzy*p.vy + Rzz*p.vz;
  // Update particle velocity
  p.vx = vx;
  p.vy = vy;
  p.vz = vz;
  
}