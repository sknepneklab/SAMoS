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
 * \file constraint_plane_walls.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Aug-2014
 * \brief Implementation of the planar walls constraint
 */ 

#include "constraint_plane_walls.hpp"

/*! Force all particles to be confined to the surface of the xy plane 
 *  by setting z = 0, vz = 0, fz = 0.
 *  At the walls, particle is stopped at the wall by setting \f$ v_x = 0 \f$
 *  and x-coordinate is set back to l or -l.
 *  \param p particle to project onto the sphere
 */
void ConstraintPlaneWalls::enforce(Particle& p)
{
  bool periodic = m_system->get_periodic();
  double ylo = m_system->get_box()->ylo, yhi = m_system->get_box()->yhi;
  double ly = yhi - ylo;
  p.z = 0.0;
  p.vz = 0.0;
  p.fz = 0.0;
  // Check periodic boundary conditions 
  if (periodic)
  {
    if (p.y <= ylo) p.y += ly;
    else if (p.y >= yhi) p.y -= ly;
  }
  if (p.x >= m_l)
  {
    p.x = m_l;
    p.vx = -p.vx;
  }
  else if (p.x <= -m_l)
  {
    p.x = -m_l;
    p.vx = -p.vx;
  }
}

/*! Rotate director of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintPlaneWalls::rotate_director(Particle& p, double phi)
{
  // Sine and cosine of the rotation angle
  double c = cos(phi), s = sin(phi);
  // Rotation matrix around z axis
  double Rxx = c, Rxy = -s;
  double Ryx = s, Ryy = c;
  // Apply rotation matrix
  double nx = Rxx*p.nx + Rxy*p.ny;
  double ny = Ryx*p.nx + Ryy*p.ny;
  double len = sqrt(nx*nx + ny*ny);
  // Update particle director
  p.nx = nx/len;
  p.ny = ny/len;
}

/*! Rotate velocity of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintPlaneWalls::rotate_velocity(Particle& p, double phi)
{
  // Sine and cosine of the rotation angle
  double c = cos(phi), s = sin(phi);
  // Rotation matrix around z axis
  double Rxx = c, Rxy = -s;
  double Ryx = s, Ryy = c;
  // Apply rotation matrix
  double vx = Rxx*p.vx + Rxy*p.vy;
  double vy = Ryx*p.vx + Ryy*p.vy;
  // Update particle director
  p.vx = vx;
  p.vy = vy;
}

/*! Project particle torque onto the normal vector (z axis). The assumption here is that 
 *  the particle's director and velocity are all already in the xy plane
 *  and that it is constrained to the plane.
 *  \param p Particle whose torque to project
*/ 
double ConstraintPlaneWalls::project_torque(Particle& p)
{
  return p.tau_z;  
}