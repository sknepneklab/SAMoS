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
 * \file integrator_langevin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Implementation of the Langevin dynamics integrator for particle position.
 */ 

#include "integrator_langevin.hpp"

void IntegratorLangevin::integrate()
{
  if (m_method == "baoab")
    this->integrate_baoab();
  else if (m_method == "spv")
    this->integrate_spv();
  else if (m_method == "bbk")
    this->integrate_bbk();
  else
    throw runtime_error("Unknown integration method");
}

// Private methods

/*! Integrates stochastic equations of motion using Langevin dynamics.
 *  The implementation is based on the BAOAB algorithm in the book
 *  "Molecular Dynamics: With Deterministic and Stochastic Numerical Methods" 
 *  by Ben Leimkuhler and Charles Matthews
 *  Interdisciplinary Applied Mathematics, Vol. 39, Springer (page 271)
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/
void IntegratorLangevin::integrate_baoab()
{
  int N = m_system->get_group(m_group_name)->get_size();
  double T = m_temp->get_val(m_system->get_run_step());
  double B = sqrt(T*(1.0-exp(-2.0*m_gamma*m_dt)));
  double exp_dt = exp(-m_gamma*m_dt);
  double dt2 = 0.5*m_dt;
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();

  // BAOA steps
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
    // Project everything back to the manifold
    m_constrainer->enforce(p);
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
    // A step
    p.x += dt2*p.vx;
    p.y += dt2*p.vy;
    p.z += dt2*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
  }
  
  // reset forces 
  m_system->reset_forces();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  
  // B step
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double fact = dt2/p.mass;
    p.vx += fact*p.fx;
    p.vy += fact*p.fy;
    p.vz += fact*p.fz;
    p.age += m_dt;
  }
  // Update vertex mesh
  m_system->update_mesh();
}

/*! Integrates stochastic equations of motion using Langevin dynamics.
 *  The implementation is based on the SPV algorithm in the book
 *  "Molecular Dynamics: With Deterministic and Stochastic Numerical Methods" 
 *  by Ben Leimkuhler and Charles Matthews
 *  Interdisciplinary Applied Mathematics, Vol. 39, Springer (page 272)
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/
void IntegratorLangevin::integrate_spv()
{
  int N = m_system->get_group(m_group_name)->get_size();
  double T = m_temp->get_val(m_system->get_run_step());
  double zeta = sqrt(T*(1.0-exp(-2.0*m_gamma*m_dt)));
  double eta;
  if (m_gamma == 0.0) eta = 0.0;
  else eta = (1.0-exp(-m_dt*m_gamma))/m_gamma;
  double exp_dt = exp(-m_dt*m_gamma);
  double dt2 = 0.5*m_dt;
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // Step 1
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    p.x += dt2*p.vx;
    p.y += dt2*p.vy;
    p.z += dt2*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
  }

  // reset forces 
  m_system->reset_forces();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);

  // Steps 2 and 3
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double fact = eta/p.mass;
    p.vx = exp_dt*p.vx + fact*p.fx;
    p.vy = exp_dt*p.vy + fact*p.fy;
    p.vz = exp_dt*p.vz + fact*p.fz;
    if (zeta != 0.0)
    {
      double stoch_fact = zeta/sqrt(p.mass);
      p.vx += stoch_fact*m_rng->gauss_rng(1.0);
      p.vy += stoch_fact*m_rng->gauss_rng(1.0);
      p.vz += stoch_fact*m_rng->gauss_rng(1.0);
    }
    // Step 3
    p.x += dt2*p.vx;
    p.y += dt2*p.vy;
    p.z += dt2*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    p.age += m_dt;
  }
}

/*! Integrates stochastic equations of motion using Langevin dynamics.
 *  The implementation is based on the BBK algorithm in the book
 *  "Molecular Dynamics: With Deterministic and Stochastic Numerical Methods" 
 *  by Ben Leimkuhler and Charles Matthews
 *  Interdisciplinary Applied Mathematics, Vol. 39, Springer (page 273)
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/
void IntegratorLangevin::integrate_bbk()
{
  int N = m_system->get_group(m_group_name)->get_size();
  double T = m_temp->get_val(m_system->get_run_step());
  double B = sqrt(2.0*T*m_dt*m_gamma);
  double dt2 = 0.5*m_dt;
  double one_m_dt2 = 1.0 - m_gamma*dt2;
  double one_div_one_p_dt2 = 1.0/(1.0 + m_gamma*dt2);
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // Steps 1 and 2
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double dt2_over_mass = dt2/p.mass;
    p.vx = one_m_dt2*p.vx + dt2_over_mass*p.fx;
    p.vy = one_m_dt2*p.vy + dt2_over_mass*p.fy;
    p.vz = one_m_dt2*p.vz + dt2_over_mass*p.fz;
    if (B != 0.0)
    {
      double stoch_fact = 0.5*B/sqrt(p.mass);
      p.vx += stoch_fact*m_Rx[i];
      p.vy += stoch_fact*m_Ry[i];
      p.vz += stoch_fact*m_Rz[i];
    }
    p.x += m_dt*p.vx;
    p.y += m_dt*p.vy;
    p.z += m_dt*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
  }

  // reset forces 
  m_system->reset_forces();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);

  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double dt2_over_mass = dt2/p.mass;
    p.vx = one_div_one_p_dt2*(p.vx + dt2_over_mass*p.fx);
    p.vy = one_div_one_p_dt2*(p.vy + dt2_over_mass*p.fy);
    p.vz = one_div_one_p_dt2*(p.vz + dt2_over_mass*p.fz);
    if (B != 0.0)
    {
      double stoch_fact = 0.5*B*one_div_one_p_dt2/sqrt(p.mass);
      m_Rx[i] = m_rng->gauss_rng(1.0);
      m_Ry[i] = m_rng->gauss_rng(1.0);
      m_Rz[i] = m_rng->gauss_rng(1.0);
      p.vx += stoch_fact*m_Rx[i];
      p.vy += stoch_fact*m_Ry[i];
      p.vz += stoch_fact*m_Rz[i];
      p.age += m_dt;
    }
  }
}
