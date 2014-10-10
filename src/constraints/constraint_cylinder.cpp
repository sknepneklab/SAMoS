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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_cylinder.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Oct-2014
 * \brief Implementation of the cylindrical constraint
 */ 

#include "constraint_cylinder.hpp"

/*! Force particle to be confined to the surface of a cylinder and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface of the cylinder
 *  \param p Particle which is to be projected onto the cylinder
 */
void ConstraintCylinder::enforce(Particle& p)
{
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double x = p.x, y = p.y;
  double R = sqrt(x*x + y*y);
  double s = m_r/R;
  // Scale back to the surface
  p.x *= s; p.y *= s; 
  // Compute unit normal
  double Nx = p.x/m_r, Ny = p.y/m_r;
  // compute v.N
  double v_dot_N = p.vx*Nx + p.vy*Ny;
  // compute n.N
  double n_dot_N = p.nx*Nx + p.ny*Ny;
  // Project velocity onto tangent plane
  p.vx -= v_dot_N*Nx; p.vy -= v_dot_N*Ny;
  // Project director onto tangent plane
  p.nx -= n_dot_N*Nx; p.ny -= n_dot_N*Ny;
  // normalize director
  double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
  p.nx *= inv_len;  p.ny *= inv_len;  p.nz *= inv_len;
  if (periodic)
  {
    if (p.z > box->zhi) p.z -= box->Lz;
    else if (p.z < box->zlo) p.z += box->Lz;
  }
}

/*! Rotate director of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto cylinder and that its director is laying in 
 *  the tangent plane.
 *  \param p Particle whose director to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintCylinder::rotate_director(Particle& p, double phi)
{
  double x = p.x, y = p.y;
  double R = sqrt(x*x + y*y);
  // Compute unit normal
  double U = x/R, V = y/R, W = 0.0;
  // Compute angle sins and cosines
  double c = cos(phi), s = sin(phi);
  // Compute new velocity coordinates
  double nx = U*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nx*c + (-W*p.ny+V*p.nz)*s;
  double ny = V*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.ny*c + ( W*p.nx-U*p.nz)*s;
  double nz = W*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nz*c + (-V*p.nx+U*p.ny)*s;
  double len = sqrt(nx*nx + ny*ny + nz*nz);
  // Update particle director (normalize it along the way to collect for any numerical drift that may have occurred)
  p.nx = nx/len;
  p.ny = ny/len;
  p.nz = nz/len;  
}

/*! Rotate velocity of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto cylinder and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintCylinder::rotate_velocity(Particle& p, double phi)
{
  double x = p.x, y = p.y, z = p.z;
  double R = sqrt(x*x + y*y + z*z);
  // Compute unit normal
  double U = x/R, V = y/R, W = z/R;
  // Compute angle sins and cosines
  double c = cos(phi), s = sin(phi);
  // Compute new velocity coordinates
  double vx = U*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vx*c + (-W*p.vy+V*p.vz)*s;
  double vy = V*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vy*c + ( W*p.vx-U*p.vz)*s;
  double vz = W*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vz*c + (-V*p.vx+U*p.vy)*s;
  // Update particle velocity
  p.vx = vx;
  p.vy = vy;
  p.vz = vz;  
}


/*! Project particle torque onto the normal vector. The assumption here is that 
 *  the particle's director and velocity are all already in the tangent plane 
 *  and that it is constrained to the cylinder.
 *  \param p Particle whose torque to project
*/ 
double ConstraintCylinder::project_torque(Particle& p)
{
  double x = p.x, y = p.y;
  // Compute unit normal
  double Nx = x/m_r, Ny = y/m_r;
  return (p.tau_x*Nx + p.tau_y*Ny);  
}