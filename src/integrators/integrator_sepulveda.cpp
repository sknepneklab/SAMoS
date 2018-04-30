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
 * \file integrator_sepulveda.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Apr-2018
 * \brief Implementation of the Sepulveda integrator.
 */ 

#include "integrator_sepulveda.hpp"

/*! This is standard velocity Verlet Sepulveda integrator.
 *  It has ability to limit maximum particle displacement 
 *  which is useful for relaxation simulations intended to
 *  remove unphysical overlaps between particles. 
**/
void IntegratorSepulveda::integrate()
{
  int N = m_system->get_group(m_group_name)->get_size();
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  double dt_2 = 0.5*m_dt;
  double B = sqrt(m_tau*m_dt);
  double theta = 1.0 - m_dt;
  
  // Generalized Langevin Dynamics, pg. 377 step 1
  // Perform first half step for velocity
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    //Particle& p = m_system->get_particle(i);
    p.vx += dt_2*(p.fx + m_eta_x[i])/p.mass;
    p.vy += dt_2*(p.fy + m_eta_y[i])/p.mass;
    p.vz += dt_2*(p.fz + m_eta_z[i])/p.mass;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constrainer->project_torque(p);
  }
  // Generalized Langevin Dynamics, pg. 377 step 2 
  // update position
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    //Particle& p = m_system->get_particle(i);
    p.x += m_dt*p.vx;
    p.y += m_dt*p.vy;
    p.z += m_dt*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
  }
  // Generalized Langevin Dynamics, pg. 377 step 3
  // Orstein-Uhlenbeck for etas
  // Note: in the book what we call m_alpha is called c_k
  // and what we call B is called alpha sqrt(2KbT c_k)
  for (int i = 0; i < N; i++)
  {
    // This computes realtive velocity with all neighbours 
    // ------------------------------
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    double dvx = 0.0, dvy = 0.0, dvz = 0.0;
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      dvx += pj.vx - p.vx;
      dvy += pj.vy - p.vy;
      dvz += pj.vz - p.vz;
    }
    // ------------------------------
    m_eta_x[i] = theta*m_eta_x[i] + (1-theta)*(-m_alpha*p.vx + dvx) + B*m_rng->gauss_rng(1.0);
    m_eta_y[i] = theta*m_eta_y[i] + (1-theta)*(-m_alpha*p.vy + dvy) + B*m_rng->gauss_rng(1.0);
    m_eta_z[i] = theta*m_eta_z[i] + (1-theta)*(-m_alpha*p.vz + dvz) + B*m_rng->gauss_rng(1.0);
  }

  // Enforce constraints and update alignment
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    //Particle& p = m_system->get_particle(i);
    // Change orientation of the velocity (in the tangent plane) according to eq. (1b)
    double dtheta = m_dt*p.omega; //m_constraint->project_torque(p);
    m_constrainer->rotate_director(p,dtheta);
    //p.omega = dtheta*m_dt;
  }

  // reset forces and torques
  m_system->reset_forces();
  m_system->reset_torques();


  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  // compute torques in the current configuration
  if (m_align)
    m_align->compute();
  
  // Generalized Langevin Dynamics, pg. 377 step 4
  // Perform second half step for velocity only if there is no limit on particle move
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    //Particle& p = m_system->get_particle(i);
    p.vx += dt_2*(p.fx + m_eta_x[i])/p.mass;
    p.vy += dt_2*(p.fy + m_eta_y[i])/p.mass;
    p.vz += dt_2*(p.fz + m_eta_z[i])/p.mass;
    
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constrainer->project_torque(p);
    p.age += m_dt;
  }
  // Update vertex mesh
  m_system->update_mesh();
}
