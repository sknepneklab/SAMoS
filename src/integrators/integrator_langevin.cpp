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
 * \file integrator_langevin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Implementation of the Langevin dynamics integrator for particle position.
 */ 

#include "integrator_langevin.hpp"

/*! Integrates stochastic equations of motion using Langevin dynamics.
 *  The implementation is based on the BAOAB algorithm in the book
 *  "Molecular Dynamics: With Deterministic and Stochastic Numerical Methods" 
 *  by Ben Leimkuhler and Charles Matthews
 *  Interdisciplinary Applied Mathematics, Vol. 39, Springer 
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/
void IntegratorLangevin::integrate()
{
  int N = m_system->get_group(m_group_name)->get_size();
  double T = m_temp->get_val(m_system->get_run_step());
  double B = sqrt(T*(1.0-exp(-2.0*m_gamma*m_dt)));
  double exp_dt = exp(-m_gamma*m_dt);
  double dt2 = 0.5*m_dt;
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // BAO steps
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double fact = dt2/p.mass;
    // B step
    p.vx += fact*p.fx;
    p.vy += fact*p.fy;
    p.vz += fact*p.fz;
    // A step
    p.x += dt2*p.vx;
    p.y += dt2*p.vy;
    p.z += dt2*p.vz;
    // O step
    p.vx *= exp_dt;
    p.vy *= exp_dt;
    p.vz *= exp_dt;
    if (B != 0.0)
    {
      double stoch_fact = B/sqrt(p.mass);
      p.vx += stoch_fact*m_rng->gauss_rng(1.0);
      p.vy += stoch_fact*m_rng->gauss_rng(1.0);
      p.vz += stoch_fact*m_rng->gauss_rng(1.0);
    }
    // Project everything back to the manifold
    m_constrainer->enforce(p);
  }
  
  // reset forces 
  m_system->reset_forces();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  
  // AB steps
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    // A step
    p.x += dt2*p.vx;
    p.y += dt2*p.vy;
    p.z += dt2*p.vz;
    // B step
    double fact = dt2/p.mass;
    p.vx += fact*p.fx;
    p.vy += fact*p.fy;
    p.vz += fact*p.fz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    p.age += m_dt;
  }
  // Update vertex mesh
  m_system->update_mesh();
}
