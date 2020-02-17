/* ***************************************************************************
 *
 *  Copyright (C) 2013-2020 University of Dundee
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
 * \file population_region.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 14-Feb-2020
 * \brief Implementation of PopulationRegion class.
 */ 

#include "population_region.hpp"


/*! Remove particle. 
 * 
 *  Particles that are outside the "allowed" region are marked for removal and removed from the system.
 * 
 *  \param t current time step 
 * 
*/
void PopulationRegion::remove(int t)
{
  bool has_region = m_has_left_bound || m_has_right_bound || m_has_top_bound || m_has_bottom_bound;
  if (m_freq > 0 && t % m_freq == 0 && has_region) // Attempt removal only at certain time steps
  { 
    double fact = m_freq*m_system->get_integrator_step();
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before cells remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    vector<int> to_remove;
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (p.in_tissue && !p.boundary)
      {
        bool remove_particle = false;
        if (m_has_left_bound && p.x < m_xmin)
          remove_particle = true;
        if (m_has_right_bound && p.x > m_xmax)
          remove_particle = true;
        if (m_has_bottom_bound && p.y < m_ymin)
          remove_particle = true;
        if (m_has_top_bound && p.y > m_ymax)
          remove_particle = true;
        if (remove_particle)
          to_remove.push_back(p.get_id());
      }
    }
    int offset = 0;
    for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end(); it++)
    {
      m_system->remove_particle((*it)-offset);
      offset++;      
    }
    if (m_system->size() == 0)
    {
      m_msg->msg(Messenger::ERROR,"Region population control. No cells left in the system. Please make sure that the allowed region is large enough.");
      throw runtime_error("No cells left in the system.");
    }
    if (!m_system->group_ok(m_group_name))
    {
      cout << "After cells remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    m_system->set_force_nlist_rebuild(true);
  }
}
