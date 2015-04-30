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
 * \file constraint_plane.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of the planar constraint
 */ 

#include "constraint_plane.hpp"

/*! Force all particles to be confined to the surface of the xy plane 
 *  by setting z = 0, vz = 0, fz = 0
 *  \param p particle to project onto the sphere
 */
void ConstraintPlane::enforce(Particle& p)
{
  bool periodic = m_system->get_periodic();
  double Lx = m_system->get_box()->Lx;
  double Ly = m_system->get_box()->Ly;
  double xlo = -0.5*Lx, xhi = 0.5*Lx;
  double ylo = -0.5*Ly, yhi = 0.5*Ly;
  p.z = 0.0;
  p.vz = 0.0;
  p.fz = 0.0;
  // Check periodic boundary conditions 
  if (periodic)
    m_system->enforce_periodic(p);
  else // reflective boundary conditions
  {
    if (p.x < xlo) 
    {
      p.x = xlo;
      p.vx = -p.vx;
    }
    else if (p.x > xhi)
    {
      p.x = xhi;
      p.vx = -p.vx;
    }
    if (p.y < ylo) 
    {
      p.y = ylo;
      p.vy = -p.vy;
    }
    else if (p.y > yhi)
    {
      p.y = yhi;
      p.vy = -p.vy;
    }
  }
  // normalize director
  p.nz = 0.0;
  double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny);
  p.nx *= inv_len;  p.ny *= inv_len;  
}

/*! Rotate director of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintPlane::rotate_director(Particle& p, double phi)
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
void ConstraintPlane::rotate_velocity(Particle& p, double phi)
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
double ConstraintPlane::project_torque(Particle& p)
{
  return p.tau_z;  
}

/*! Rescale box size and make sure that all particles fit in it.
 *  Rescaling is done only at certain steps and only if rescale 
 *  factor is not equal to 1.
*/
bool ConstraintPlane::rescale()
{
  if (m_rescale != 1.0)
  {
    int step = m_system->get_step();
    if ((step % m_rescale_freq == 0) && (step < m_rescale_steps))
    {
      m_system->get_box()->rescale(m_scale);
      for  (int i = 0; i < m_system->size(); i++)
      {
        Particle& p = m_system->get_particle(i);
        p.x *= m_scale; 
        p.y *= m_scale; 
        p.z *= m_scale;
        this->enforce(p);
      }
      return true;
    }
  }
  return false;
}

