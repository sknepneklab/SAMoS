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
 * \file integrator_nve.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Feb-2014
 * \brief Implementation of the NVE integrator.
 */ 

#include "integrator_nve.hpp"

/*! This is standard velocity Verlet NVE integrator.
 *  It has ability to limit maximum particle displacement 
 *  which is useful for relaxation simulations intended to
 *  remove unphysical overlaps between particles. 
**/
void IntegratorNVE::integrate()
{
  int N = m_system->size();
  double dt_2 = 0.5*m_dt;
  // Perform first half step for velocity
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    p.vx += dt_2*p.fx;
    p.vy += dt_2*p.fy;
    p.vz += dt_2*p.fz;
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constraint->project_torque(p);
  }
  // update position
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    double dx = m_dt*p.vx;
    double dy = m_dt*p.vy;
    double dz = m_dt*p.vz;
    if (m_has_limit)
    {
      double dr = sqrt(dx*dx + dy*dy +dz*dz);
      if (dr > m_limit)
      {
        dx = m_limit*dx/dr;
        dy = m_limit*dy/dr;
        dz = m_limit*dz/dr;
      }
    }
    p.x += dx;
    p.y += dy; 
    p.z += dz;
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Change orientation of the director (in the tangent plane) according to eq. (1b)
    double dtheta = m_dt*p.omega;
    m_constraint->rotate_director(p,dtheta);
  }
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute();
  // compute torques in the current configuration
  if (m_align)
    m_align->compute();
  // Enforce constraints and update alignment
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    // Change orientation of the velocity (in the tangent plane) according to eq. (1b)
    double dtheta = m_dt*m_constraint->project_torque(p);
    if (m_has_theta_limit)
      if (fabs(dtheta) > m_theta_limit)
        dtheta = SIGN(dtheta)*m_theta_limit;
    m_constraint->rotate_director(p,dtheta);
    p.omega = dtheta*m_dt;
  }
  // Perform second half step for velocity only if there is no limit on particle move
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    p.vx += dt_2*p.fx;
    p.vy += dt_2*p.fy;
    p.vz += dt_2*p.fz;
    // we also need to limit velocity to it does not go crazy
    if (m_has_limit)
    {
      double v = sqrt(p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
      if (v*m_dt > m_limit)
      {
        p.vx = p.vx/v*m_limit/m_dt;
        p.vy = p.vx/v*m_limit/m_dt;
        p.vz = p.vx/v*m_limit/m_dt;
      }
    }
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constraint->project_torque(p);
  }
}