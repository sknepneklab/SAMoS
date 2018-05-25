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
 * \file pair_yukawa_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-May-2018
 * \brief Implementation of PairYukawaPotential class
 */ 

#include "pair_yukawa_potential.hpp"

void PairYukawaPotential::compute(double dt)
{
  int N = m_system->size();
  double g = m_g;
  double kappa = m_kappa;
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
  
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("yukawa",0.0);
    }
  }

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= rcut_sq)
      {
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          g = m_pair_params[pi_t][pj_t].g;
          kappa = m_pair_params[pi_t][pj_t].kappa;
        }
        double r = sqrt(r_sq);
        double exp_fact = g*exp(-kappa*r);
        // Handle potential 
        double potential_energy = exp_fact/r;
        m_potential_energy += potential_energy;
        // Handle force
        double force_factor = exp_fact*(1 + kappa*r)/(r_sq*r);
        pi.fx += force_factor*dx;
        pi.fy += force_factor*dy;
        pi.fz += force_factor*dz;
        // Use 3d Newton's law
        pj.fx -= force_factor*dx;
        pj.fy -= force_factor*dy;
        pj.fz -= force_factor*dz;
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("yukawa",potential_energy);
          pj.add_pot_energy("yukawa",potential_energy);
        }
      }
    }
  }
}
