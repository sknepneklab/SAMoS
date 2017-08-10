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
 * \file pair_coulomb_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Implementation of PairCoulombPotential class
 */ 

#include "pair_coulomb_potential.hpp"

void PairCoulombPotential::compute(double dt)
{
  int N = m_system->size();
  double sigma = m_sigma;
  double alpha = m_alpha;
  double sigma_sq = sigma*sigma;
  double phase_fact_i = 1.0;  // phase in factor for particle i
  double phase_fact_j = 1.0;  // phase in factor for particle j
  double phase_fact = 1.0;    // phase in factor for pair interaction (see below)
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("coulomb",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      phase_fact_i = 0.5*(1.0 + m_val->get_val(static_cast<int>(pi.age/dt)));
    for (int j = i+1; j < N; j++)
    {
      Particle& pj = m_system->get_particle(j);
      if (m_phase_in)
      {
        phase_fact_j = 0.5*(1.0 + m_val->get_val(static_cast<int>(pj.age/dt)));
        // Determine global phase in factor: particles start at 0.5 strength (both daughters of a division replace the mother)
        // Except for the interaction between daughters which starts at 0
        if ( phase_fact_i < 1.0 && phase_fact_j < 1.0)
	        phase_fact=phase_fact_i + phase_fact_j - 1.0;
		    else 
	        phase_fact = phase_fact_i*phase_fact_j;
      }
      if (m_has_pair_params)
      {
        int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
        alpha = m_pair_params[pi_t][pj_t].alpha;
        sigma = m_pair_params[pi_t][pj_t].sigma;
        sigma_sq = sigma*sigma;
      }
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r_sq);
      double inv_r_sq = sigma_sq/r_sq;
      double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
      // Handle potential 
      double potential_energy = phase_fact*(alpha/r + 4.0*fabs(alpha)*inv_r_6*inv_r_6);
      m_potential_energy += potential_energy;
      // Handle force
      double r_3 = r*r_sq;
      double force_factor = phase_fact*(alpha/r_3 + 48.0*fabs(alpha)*inv_r_6*inv_r_6*inv_r_sq);
      pi.fx -= force_factor*dx;
      pi.fy -= force_factor*dy;
      pi.fz -= force_factor*dz;
      // Use 3d Newton's law
      pj.fx += force_factor*dx;
      pj.fy += force_factor*dy;
      pj.fz += force_factor*dz;
      if (m_system->compute_per_particle_energy())
      {
        pi.add_pot_energy("coulomb",potential_energy);
        pj.add_pot_energy("coulomb",potential_energy);
      }
    }
  }
}
