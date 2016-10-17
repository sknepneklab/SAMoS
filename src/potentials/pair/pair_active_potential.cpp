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
 * \file pair_active_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Mar-2015
 * \brief Implementation of PairActivePotential class
 */ 

#include "pair_active_potential.hpp"

void PairActivePotential::compute(double dt)
{
  int N = m_system->size();
  double alpha = m_alpha;
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
  double potential_energy = 0.0;
   
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("active",0.0);
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
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= rcut_sq)
      {
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          alpha = m_pair_params[pi_t][pj_t].alpha;
        }
        double dot = pi.nx*pj.nx + pi.ny*pj.ny + pi.nz*pj.nz;
        if (dot < 0.0)
        {
          double dnx = pi.nx - pj.nx;
          double dny = pi.ny - pj.ny;
          double dnz = pi.nz - pj.nz;
          double len_dn = sqrt(dnx*dnx + dny*dny + dnz*dnz);
          // Handle force
          double force_factor = 0.5*alpha*dot/len_dn;
          pi.fx += force_factor*dnx;
          pi.fy += force_factor*dny;
          pi.fz += force_factor*dnz;
          // Use 3d Newton's law
          pj.fx -= force_factor*dnx;
          pj.fy -= force_factor*dny;
          pj.fz -= force_factor*dnz;
          if (m_system->compute_per_particle_energy())
          {
            pi.add_pot_energy("active",potential_energy);
            pj.add_pot_energy("active",potential_energy);
          } 
        }
      }
    }
  }
}
