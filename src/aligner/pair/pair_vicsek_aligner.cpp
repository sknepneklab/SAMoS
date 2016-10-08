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
 * \file pair_vicsek_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Jan-2014
 * \brief Implementation of PairVicsekPotential class
 */ 

#include "pair_vicsek_aligner.hpp"

void PairVicsekAlign::compute()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double rcut = m_rcut;
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      if (periodic)
      {
        if (dx > box->xhi) dx -= box->Lx;
        else if (dx < box->xlo) dx += box->Lx;
        if (dy > box->yhi) dy -= box->Ly;
        else if (dy < box->ylo) dy += box->Ly;
        if (dz > box->zhi) dz -= box->Lz;
        else if (dz < box->zlo) dz += box->Lz;
      }
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= rcut*rcut)
      {
        pi.tau_x += pj.vx;
        pi.tau_y += pj.vy;
        pi.tau_z += pj.vz;
        // Since we use 3d Newton's and only have half the the neighbour list we need handle the other part as well
        pj.tau_x += pi.vx;
        pj.tau_y += pi.vy;
        pj.tau_z += pi.vz;
      }
    }
    pi.tau_x += pi.vx;
    pi.tau_y += pi.vy;
    pi.tau_z += pi.vz;
  }
  // Now we need to normalize all taus
  for  (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    double v = sqrt(p.tau_x*p.tau_x + p.tau_y*p.tau_y + p.tau_z*p.tau_z);
    if (v == 0.0) v = 1.0;
    double inv_v = 1.0/v;
    p.tau_x *= inv_v;
    p.tau_y *= inv_v;
    p.tau_z *= inv_v;
  }
}
