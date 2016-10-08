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
 * \file population_actomyosin_poisson.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Jul-2016
 * \brief Implementation of PopulationActomyosinPoisson class.
 */ 

#include "population_actomyosin_poisson.hpp"


/*! This function controls attachement and detachement of myosin beads. This is a
 *  achieved via chaning types of myosin head groups from attached (type "A") to detached (type "D")a
 *  and vice versa. In this implementation we use a Poisson process. 
 *
 *  \note In order to avoid head groups, e.g. chainging  their type from "A" to "D" and back to "A" in the 
 *  single time step, we implement both processes in the same function using lists of indices.
 *  \param t current time step
 *  
*/
void PopulationActomyosinPoisson::divide(int t)
{
  if ((m_freq > 0) && (t % m_freq == 0) && (m_attach_rate > 0.0 || m_detach_rate > 0.0))  // Attempt D to A transition only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before actomyosin division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->size();
    vector<int> to_attach;
    vector<int> to_detach;
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(i);                
      if (pi.get_type() == m_type_actin) 
      {
        vector<int>& neigh = m_nlist->get_neighbours(i);
        double min_dist = 1e10;
        int min_id = 0; 
        for (unsigned int j = 0; j < neigh.size(); j++)  // Loop over all neigbours of actin bead i
        {
          Particle& pj = m_system->get_particle(neigh[j]);         
          if (pj.get_type() == m_type_d)   // if the neighbour if of type "detached"
          {
            double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
            m_system->apply_periodic(dx,dy,dz);
            double r = sqrt(dx*dx + dy*dy + dz*dz);  // compute distance to the neighbour
            if (r < min_dist)   // if this distance is below current minimum distance record it and the id of the bead itself
            {
              min_dist = r;
              min_id = pj.get_id();
            }
          }
        }
        if (min_dist < m_re)  // If the closest "detached" bead is within cutoff distance re
          to_attach.push_back(min_id);
          //if (find(to_attach.begin(),to_attach.end(),min_id) == to_attach.end())
          //  to_attach.push_back(min_id);
      }
      if (pi.get_type() == m_type_a)
        to_detach.push_back(pi.get_id());
    }
    for (vector<int>::iterator it_a = to_attach.begin(); it_a != to_attach.end(); it_a++)
    {
      Particle& pi = m_system->get_particle(*it_a);
      double prob = 1.0 - exp(-m_attach_rate*pi.age);
      if (m_rng->drnd() < prob)  // flip its type to "attached" with probability attach_prob.
      {
        pi.set_type(m_type_a);
        pi.age = 0.0;
      }
    }
    for (vector<int>::iterator it_d = to_detach.begin(); it_d != to_detach.end(); it_d++)
    {
      Particle& pi = m_system->get_particle(*it_d);
      double f = sqrt(pi.fx*pi.fx + pi.fy*pi.fy + pi.fz*pi.fz);
      double prob = 1.0 - exp(-(m_detach_rate+m_lambda*f)*pi.age);
      if (m_rng->drnd() < prob)  // flip its type to "attached" with probability attach_prob.
      {
        pi.set_type(m_type_d);  
        pi.age = 0.0;
      }
    }
  }
}

