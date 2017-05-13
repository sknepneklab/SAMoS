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
  double alpha = m_alpha;   
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

  if (m_system->record_force_type())
    this->reset_force_types("active_nematic");

  // Reset total potential energy to zero
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    double Qi_xx = pi.nx*pi.nx - 1.0/3.0, Qi_xy = pi.nx*pi.ny,           Qi_xz = pi.nx*pi.nz;
    double Qi_yx = Qi_xy,                 Qi_yy = pi.ny*pi.ny - 1.0/3.0, Qi_yz = pi.ny*pi.nz;
    double Qi_zx = Qi_xz,                 Qi_zy = Qi_yz,                 Qi_zz = pi.nz*pi.nz - 1.0/3.0;
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
        double Qj_xx = pj.nx*pj.nx - 1.0/3.0, Qj_xy = pj.nx*pj.ny,           Qj_xz = pj.nx*pj.nz;
        double Qj_yx = Qj_xy,                 Qj_yy = pj.ny*pj.ny - 1.0/3.0, Qj_yz = pj.ny*pj.nz;
        double Qj_zx = Qj_xz,                 Qj_zy = Qj_yz,                 Qj_zz = pj.nz*pj.nz - 1.0/3.0; 
        double factor = alpha/r_sq;
        double fx = factor*((Qj_xx - Qi_xx)*dx + (Qj_xy - Qi_xy)*dy + (Qj_xz - Qi_xz)*dz);
        double fy = factor*((Qj_yx - Qi_yx)*dx + (Qj_yy - Qi_yy)*dy + (Qj_yz - Qi_yz)*dz);
        double fz = factor*((Qj_zx - Qi_zx)*dx + (Qj_zy - Qi_zy)*dy + (Qj_zz - Qi_zz)*dz);
        pi.fx += fx;
        pi.fy += fy;
        pi.fz += fz;
        if (m_system->record_force_type())
          pi.add_force_type("active_nematic",fx,fy,fz);
      }
    }
  }
}
