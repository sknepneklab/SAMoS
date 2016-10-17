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
 * \file bond_active_force.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Aug-2015
 * \brief Implementation of BondActiveForce class
 */ 

#include "bond_active_force.hpp"

void BondActiveForce::compute()
{
  int Nbonds = m_system->num_bonds();
  bool periodic = m_system->get_periodic();
  double f = m_f;
    
  for  (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    Particle& pi = m_system->get_particle(b.i);
    Particle& pj = m_system->get_particle(b.j);
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
    if (periodic)
      m_system->apply_periodic(dx,dy,dz);
    if (m_has_bond_params)
    {
      f = m_bond_params[b.type-1].f;
    }
      
    double r_sq = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r_sq);
    // Handle force
    double force_factor = f/r;
    pi.fx += force_factor*dx;
    pi.fy += force_factor*dy;
    pi.fz += force_factor*dz;
    // Note that both beads get the same force. This is not a pair force and thus 
    // 3d Newton's law does not hold for the pair of beads.
    pj.fx += force_factor*dx;
    pj.fy += force_factor*dy;
    pj.fz += force_factor*dz;
  }
}
