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
 * \file bond_harmonic_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Implementation of BondHarmonicPotential class
 */ 

#include "bond_harmonic_potential.hpp"

void BondHarmonicPotential::compute()
{
  int Nbonds = m_system->num_bonds();
  BoxPtr box = m_system->get_box();
  double k = m_k;
  double l0 = m_l0;
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
    m_system->apply_periodic(dx,dy,dz);
    if (m_has_bond_params)
    {
      k = m_bond_params[b.type-1].k;
      l0 = m_bond_params[b.type-1].l0;
     }
      
    double r_sq = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r_sq);
    double dl = r - l0;
    // Handle potential 
    double potential_energy = 0.5*k*dl*dl;
    m_potential_energy += potential_energy;
    // Handle force
    double force_factor = k*dl/r;
    pi.fx += force_factor*dx;
    pi.fy += force_factor*dy;
    pi.fz += force_factor*dz;
    // Use 3d Newton's law
    pj.fx -= force_factor*dx;
    pj.fy -= force_factor*dy;
    pj.fz -= force_factor*dz;
  }
}
