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
  if (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end())
  {
    double x = p.x, y = p.y, z = p.z;
    double R = sqrt(x*x + y*y + z*z);
    double s = m_r/R;
    // Scale back to the surface
    p.x *= s; p.y *= s; p.z *= s;
    // Compute unit normal
    double Nx = p.x/m_r, Ny = p.y/m_r, Nz = p.z/m_r;
    // compute v.N
    double v_dot_N = p.vx*Nx + p.vy*Ny + p.vz*Nz;
    // compute n.N
    double n_dot_N = p.nx*Nx + p.ny*Ny + p.nz*Nz;
    // Project velocity onto tangent plane
    p.vx -= v_dot_N*Nx; p.vy -= v_dot_N*Ny; p.vz -= v_dot_N*Nz;
    // Project director onto tangent plane
    p.nx -= n_dot_N*Nx; p.ny -= n_dot_N*Ny; p.nz -= n_dot_N*Nz;
    // normalize director
    double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
    p.nx *= inv_len;  p.ny *= inv_len;  p.nz *= inv_len;
  }
}

/*! Rescale sphere size and make sure that all particles are still on it.
 *  Rescaling is done only at certain steps and only if rescale 
 *  factor is not equal to 1.
 */
bool ConstraintSphere::rescale()
{
  if (m_rescale != 1.0)
  {
    int step = m_system->get_step();
    if ((step % m_rescale_freq == 0) && (step < m_rescale_steps))
    {
      m_r *= m_scale;
      for  (int i = 0; i < m_system->size(); i++)
      {
        Particle& p = m_system->get_particle(i);
        this->enforce(p);
      }
      return true;
    }
  }
  return false;
}

