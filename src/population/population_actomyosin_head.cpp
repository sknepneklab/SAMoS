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
 * \file population_actomyosin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Oct-2016
 * \brief Implementation of PopulationActomyosinHead class.
 */ 

#include "population_actomyosin_head.hpp"


/*! This function controls attachment and detachment of myosin beads. This is a
 *  achieved via changing types of myosin head groups from attached (type "A") to detached (type "D")a
 *  and vice versa. 
 *  This is base on the exact time of attachment/detachment 
 *  
*/
void PopulationActomyosinHead::divide(int t)
{
  if ((m_freq > 0) && (t % m_freq == 0) && (m_attach_rate > 0.0 || m_detach_rate > 0.0))  // Attempt D to A transition only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before actomyosin head division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->size();
    vector<int> to_attach;
    vector<int> to_detach;
    // Probability of attachment/detachment for a given particle is rate per particle multiplied with time,
    // where time is equal to m_freq*integrator_time_step.   
    //double attach_prob = m_attach_rate*m_freq*m_system->get_integrator_step();
    //double detach_prob = m_detach_rate*m_freq*m_system->get_integrator_step();
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
          if ((pj.get_type() == m_type_d) && (pj.bind <= t))   // if the neighbour if of type "detached" and is about to be attached; 
          // I think <= is important to make sure that there are not "left-over"" beads. 
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
      if ((pi.get_type() == m_type_a) && (pi.unbind <= t)) // Same logic here for the <= instead of ==
        to_detach.push_back(pi.get_id());
    }
    for (vector<int>::iterator it_a = to_attach.begin(); it_a != to_attach.end(); it_a++)
      {
        Particle& pi = m_system->get_particle(*it_a);
        pi.set_type(m_type_a);
        pi.unbind = t+static_cast<int>(-m_mean_detach*log(m_rng->drnd()));    // next unbind is at the current time + log distribution 
      }
    for (vector<int>::iterator it_d = to_detach.begin(); it_d != to_detach.end(); it_d++)
    {
      Particle& pi = m_system->get_particle(*it_d);
      // TODO: Figure out how to incorporate forces
      //double f = sqrt(pi.fx*pi.fx + pi.fy*pi.fy + pi.fz*pi.fz);
      //double prob = m_detach_prob*exp(m_lambda*f);
      //double f = sqrt(pi.fx*pi.fx + pi.fy*pi.fy + pi.fz*pi.fz);
      //double detach_prob_fdep = detach_prob*exp(m_lambda*f);
      pi.set_type(m_type_d);
      pi.bind = t+static_cast<int>(-m_mean_attach*log(m_rng->drnd()));  // next bind is at the current time + log distribution 
    }
  }
}

