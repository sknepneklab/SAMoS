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
 * \file integrator_fire.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Feb-2014
 * \brief Implementation of the FIRE integrator.
 */ 

#include "integrator_fire.hpp"

/*! This is the FIRE minimiser.  
**/
void IntegratorFIRE::integrate()
{
  double P = 0.0;
  double E; 
  double Fnorm = 0.0;
  double Vnorm = 0.0;
  double vx, vy, vz;
  double ax, ay, az;

  int N = m_system->get_group(m_group_name)->get_size();
  double sqrt_ndof = sqrt(3*N);
  vector<int> particles = m_system->get_group(m_group_name)->get_particles();
  double dt_2 = 0.5*m_dt;
  
  // Perform first half step for velocity
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    p.vx += dt_2*p.fx;
    p.vy += dt_2*p.fy;
    p.vz += dt_2*p.fz;
    // Update position
    p.x += m_dt*p.vx;
    p.y += m_dt*p.vy;
    p.z += m_dt*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constrainer->project_torque(p);
    double dtheta = m_dt*p.omega; //m_constraint->project_torque(p);
    m_constrainer->rotate_director(p,dtheta);
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

  // Perform second half step for velocity 
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    p.vx += dt_2*p.fx;
    p.vy += dt_2*p.fy;
    p.vz += dt_2*p.fz;
    // Compute FIRE related quantities
    P += p.vx*p.fx + p.vy*p.fy + p.vz*p.fz;
    Fnorm += p.fx*p.fx + p.fy*p.fy + p.fz*p.fz;
    Vnorm += p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
    // Project everything back to the manifold
    m_constrainer->enforce(p);
    // Update angular velocity
    p.omega += dt_2*m_constrainer->project_torque(p);
  }
  // Update vertex mesh
  m_system->update_mesh();

  Fnorm = sqrt(Fnorm);
  Vnorm = sqrt(Vnorm);  

  E = m_potential->compute_potential_energy();
  cout << "E = " << E << "  Fnorm  = " << Fnorm/sqrt_ndof << "  Vnorm = " << Vnorm << endl;

  if (Fnorm/sqrt_ndof < m_F_tol && fabs(E - m_old_energy) < m_E_tol)
  {
    m_converged = true;
    m_msg->msg(Messenger::INFO,"FIRE minimisation converged.");
    return;
  }
  
  double inv_Fnorm = 1.0/Fnorm;
  double fact_1 = 1.0 - m_alpha;
  double fact_2 = m_alpha * inv_Fnorm * Vnorm;
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_system->get_particle(pi);
    p.vx = fact_1*p.vx + fact_2 * p.fx;
    p.vy = fact_1*p.vy + fact_2 * p.fy;
    p.vz = fact_1*p.vz + fact_2 * p.fz;
  }
  
  if (P > 0.0)
  {
    m_last_neg++;
    if (m_last_neg > m_N_min)
    {
      m_dt = min(m_dt*m_f_inc, m_dt_max);
      m_alpha *= m_f_alpha;
    }
  }
  else
  {
    m_dt *= m_f_dec;
    m_alpha = m_alpha_init;
    m_last_neg = 0;
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
    }
  }
  
  m_converged = false;
  m_old_energy = E;

}
