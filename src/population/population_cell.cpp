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

/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file population_cell.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Dec-2015
 * \brief Implementation of PopulationCell class.
 */ 

#include "population_cell.hpp"


/*! Divide cells according to their current area.
 *  We draw a uniform random number \f$ \zeta \f$ between 0 and 1
 *  if \f$ \zeta < \exp\left(\gamma(A_i-A_0)\right) \f$, split particle.
 *  Here, \f$ A_i \f$ is current cell area, \f$ A_0 \f$ if the native cell area 
 *  (which changes during growth) and \f$ \gamma \f$ is a constant.  
 * 
 *  After the division, each daughter cell has native area equal to the half of the 
 *  native area of the mother cell. 
 *
 *  Position of daughter cells is determined by the orientation vector n.
 *  \param t current time step
 *  
*/
void PopulationCell::divide(int t)
{
  if (m_freq > 0 && t % m_freq == 0)  // Attempt division only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    double fact = m_freq*m_div_rate*m_system->get_integrator_step();
    Mesh& mesh = m_system->get_mesh();
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi); 
      Vertex& V = mesh.get_vertices()[p.get_id()];
	    // actual probability of dividing now: (attempt_freq * dt) * m_div_rate * (V.area - max_area) 
	    // here m_div_rate is a scaling factor which controls the actual division rate
      if (p.in_tissue && !p.boundary && V.area > m_max_A0)
      { 
        double prob_div = fact*(V.area - m_max_A0); // Bell model of division
        if (m_rng->drnd() < prob_div)  // Only internal verices can divide
        {
          //cout << t << " " << V.area << " " << p.A0 << " " << exp((V.area-p.A0)/m_div_rate) << endl;
          Particle p_new(m_system->size(), p.get_type(), p.get_radius());
          p_new.x = p.x + m_alpha*m_split_distance*p.get_radius()*p.nx;
          p_new.y = p.y + m_alpha*m_split_distance*p.get_radius()*p.ny;
          p_new.z = p.z + m_alpha*m_split_distance*p.get_radius()*p.nz;
          p_new.set_parent(p.get_flag());
          m_system->apply_periodic(p_new.x,p_new.y,p_new.z);
          
          p.x -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nx;
          p.y -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.ny;
          p.z -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nz;
          p.age = 0.0;
          p.A0 = p.get_A0();
          m_system->apply_periodic(p.x,p.y,p.z);
          
          p_new.nx = p.nx; p_new.ny = p.ny; p_new.nz = p.nz;
          p_new.vx = p.vx; p_new.vy = p.vy; p_new.vz = p.vz;
          p_new.Nx = p.Nx; p_new.Ny = p.Ny; p_new.Nz = p.Nz;
          p_new.age = 0.0;
          p_new.A0 = p.get_A0();
          p_new.set_radius(p.get_radius());
          p_new.set_type(p.get_type());
          p_new.set_default_area(p.get_A0());
          p_new.in_tissue = true;
          for(list<string>::iterator it_g = p.groups.begin(); it_g != p.groups.end(); it_g++)
            p_new.add_group(*it_g);
          m_system->add_particle(p_new);
        }
      }
    }
    if (!m_system->group_ok(m_group_name))
    {
      cout << "After cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    m_system->set_force_nlist_rebuild(true);
  }
}

/*! Remove particle. 
 * 
 *  In biological systems all cells die. Here we assume that the death is random, but 
 *  proportional to the particles age.
 * 
 *  \param t current time step
 * 
*/
void PopulationCell::remove(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_death_rate > 0.0)  // Attempt removal only at certain time steps
  { 
    double fact = m_freq*m_system->get_integrator_step();
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before Cell Remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    vector<int> to_remove;
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      // actual probability of dying now: (attempt_freq * dt) exp[(age-max_age)*m_death_rate] 
      // In the limits where we can linearise this it gives:
      // P_div = (attempt_freq * dt) [1+m_death_rate*(age-age_max)]. This death rate is an inverse time scale
      //double prob_death =m_freq*m_system->get_integrator_step()*exp((p.age-m_max_age)*m_death_rate); 
      // Trying a very simple, linearly increasing death chance
      //double prob_death = fact*p.age/m_max_age;
      double prob_death = m_death_rate*m_freq*m_system->get_integrator_step(); // actual probability of dividing now: rate * (attempt_freq * dt)
      if (p.in_tissue && !p.boundary && m_rng->drnd() < prob_death)
          to_remove.push_back(p.get_id());
    }
    int offset = 0;
    for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end(); it++)
    {
      m_system->remove_particle((*it)-offset);
      offset++;      
    }
    if (m_system->size() == 0)
    {
      m_msg->msg(Messenger::ERROR,"Cell population control. No cells left in the system. Please reduce that death rate.");
      throw runtime_error("No cells left in the system.");
    }
    if (!m_system->group_ok(m_group_name))
    {
      cout << "After Cell Remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    m_system->set_force_nlist_rebuild(true);
  }
}

/*! Grow (rescale) cell native area by a given amount.
 * 
 *  \param t current time step
 *  
*/
void PopulationCell::grow(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_growth_rate > 0.0) 
  { 
    double fact = m_freq*m_system->get_integrator_step()*m_growth_rate;
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    for (int i = 0; i < N; i++)
    {
	  // Growth probability stays dimensionless, between 0 and 1. Instead, the actual growth rate is no an inverse time
      if (m_rng->drnd() < m_growth_prob)
      {
        int pi = particles[i];
        Particle& p = m_system->get_particle(pi); 
        if (p.in_tissue)
          p.A0 *= (1.0+fact);
      }
    }
    m_system->set_force_nlist_rebuild(true);
    if (m_rescale_contacts)
      m_system->set_nlist_rescale(sqrt(1.0+fact));
  }
}
