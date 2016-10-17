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
 * \file integrator_brownian_pos.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Implementation of the Brownian dynamics integrator for particle position.
 */ 

#include "integrator_brownian_pos.hpp"

/*! Integrates equation of motion in the over-damped limit using a first order 
 *  scheme. 
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/
void IntegratorBrownianPos::integrate()
{
  int N = m_system->get_group(m_group_name)->get_size();
  double T = m_temp->get_val(m_system->get_run_step());
  double B = sqrt(2.0*m_mu*T);
  double sqrt_dt = sqrt(m_dt);
  double fr_x = 0.0, fr_y = 0.0, fr_z = 0.0;  // Random part of the force
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  
  // reset forces 
  m_system->reset_forces();
  
  // compute forces in the current configuration
  if (m_potential)
    m_potential->compute(m_dt);
  // iterate over all particles 
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    // Update velocity
    p.vx = m_mu*p.fx;
    p.vy = m_mu*p.fy;
    p.vz = m_mu*p.fz;
    // Update particle position 
    p.x += m_dt*p.vx;
    p.y += m_dt*p.vy;
    p.z += m_dt*p.vz;
    // Check is non-zero T and if non-zero add stochastic part
    if (T > 0.0)
    {
      fr_x = B*m_rng->gauss_rng(1.0);
      fr_y = B*m_rng->gauss_rng(1.0);
      fr_z = B*m_rng->gauss_rng(1.0);
      p.vx += fr_x; 
      p.vy += fr_y;
      p.vz += fr_z;  
      p.x += sqrt_dt*fr_x;
      p.y += sqrt_dt*fr_y;
      p.z += sqrt_dt*fr_z;
    }
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    p.age += m_dt;
  }
  // Update vertex mesh
  m_system->update_mesh();
}
