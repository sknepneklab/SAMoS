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
 * \file integrator_vicsek.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Dec-2013
 * \brief Implementation of the Vicsek dynamics integrator.
 */ 

#include "integrator_vicsek.hpp"

/*! This function integrates equations of motion as introduced in Vicsek model
 *  see, e.g., Eqs. (13) and (14) in T. Vicsek, A. Zafeiris, Phys. Reports 517 (2012) p. 71-140.
 *  Actual equations are a bit modified in order to introduce spacial interaction term \f$ \vec F_i = \sum_{i\neq j}\vec F_{ij} \f$.
 *  \f$ \vec v_i(t+1) = v_0\frac{\langle\vec v_j(t)\rangle}{\left|\langle\vec v_j(t)\rangle\right|} + \mu \sum_{i\neq j}\vec F_{ij} + \mathrm{noise} \f$ and
 *  \f$ \vec r_i(t+1) = \vec r_i(t) + \vec v_i(t+1) \f$. \f$ \langle\vec v_j(t)\rangle \f$ is the average 
 *   velocity of the neighbouring particles in some radius and \f$ v_0 \f$ is magnitude of the velocity.
 *   \f$ \mu \f$ is the mobility.
 *   We note that Viscek model has no explicit time step, i.e. it is set to one.
 *   Noise term is essentially a random rotation of the particle's velocity around its normal vector 
 *   by angle \f$ \delta\vartheta \f$ where \f$ \delta\vartheta \f$ is a random variable drawn uniformly 
 *   from the interval \f$ [-\eta\pi\dots\eta\pi] \f$, where \f$ \eta \f$ is a constant.
**/
void IntegratorVicsek::integrate()
{
  double noise = m_eta*sqrt(m_dt);
  int N = m_system->get_group(m_group_name)->get_size();
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // reset forces and torques
  m_system->reset_forces();
  m_system->reset_torques();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  // compute torques in the current configuration
  if (m_align)
    m_align->compute();
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    // Update particle velocity 
    p.vx = m_v0*p.tau_x + m_mu*p.fx;
    p.vy = m_v0*p.tau_y + m_mu*p.fy;
    p.vz = m_v0*p.tau_z + m_mu*p.fz;
    // Project everything back to the manifold
    m_constraint->enforce(p);
    // Change orientation of the velocity (in the tangent plane) 
    double theta = 2.0*noise*M_PI*(m_rng->drnd() - 0.5);
    m_constraint->rotate_velocity(p,theta);
    // Update particle position 
    p.x += m_dt*p.vx;
    p.y += m_dt*p.vy;
    p.z += m_dt*p.vz;
    // Project everything back to the manifold
    m_constraint->enforce(p);
    p.age += m_dt;
  }
}
