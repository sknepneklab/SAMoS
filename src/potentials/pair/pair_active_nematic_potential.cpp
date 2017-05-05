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
 * \file pair_active_nematic_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-May-2017
 * \brief Implementation of PairActivePotential class
 */ 

#include "pair_active_nematic_potential.hpp"

void PairActiveNematicPotential::compute(double dt)
{
  int N = m_system->size();
  double alpha = -m_alpha; // minus comes from the defintion 
  double rcut = m_rcut;
  double rcut_sq = rcut*rcut;
  double potential_energy = 0.0;
   
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("active_nematic",0.0);
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
          alpha = -m_pair_params[pi_t][pj_t].alpha; // minus comes from the defintion 
        }
        double Qxx = pj.nx*pj.nx - 1.0/3.0, Qxy = pj.nx*pj.ny,           Qxz = pj.nx*pj.nz;
        double Qyx = Qxy,                   Qyy = pj.ny*pj.ny - 1.0/3.0, Qyz = pj.ny*pj.nz;
        double Qzx = Qxz,                   Qzy = Qyz,                   Qzz = pj.nz*pj.nz - 1.0/3.0; 
        pi.fx += alpha*(Qxx*dx + Qxy*dy + Qxz*dz)/r_sq;
        pi.fy += alpha*(Qyx*dx + Qyy*dy + Qyz*dz)/r_sq;
        pi.fz += alpha*(Qzx*dx + Qzy*dy + Qzz*dz)/r_sq;
      }
    }
  }
}
