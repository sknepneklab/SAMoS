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
 * \file population_density.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Implementation of PopulationDensity class.
 */ 

#include "population_density.hpp"


/*! Divide particles according to their age.
 *  We draw a uniform random number between 0 and 1
 *  if this number is less than particle age divided by
 *  the division rate, split particle. 
 * 
 *  Particle is split along the direction vector n
 *  New particles are placed one radius apart. Original 
 *  particle is pushed back by 1/2r and the new one is
 *  placed at 1/2r.
 * 
 *  \param t current time step
 *  
*/
void PopulationDensity::divide(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_div_rate > 0.0)  // Attempt division only at certain time steps
  { 
    cout << "Handling divisions for group " << m_group_name << endl;
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before divide P: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int new_type;   // type of new particle
    double new_r;   // radius of newly formed particle
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    BoxPtr box = m_system->get_box();
    double prob_div = m_div_rate*m_freq*m_system->get_integrator_step(); // actual probability of dividing now: rate * (attempt_freq * dt)
    if (prob_div > 1.0)
    {
      m_msg->msg(Messenger::ERROR,"Division rate "+lexical_cast<string>(prob_div)+" is too large for current time step and attempt rate.");
      m_msg->msg(Messenger::INFO,"We have: division "+lexical_cast<string>(m_div_rate)+" frequency "+lexical_cast<string>(m_freq)+" time step "+lexical_cast<string>(m_system->get_step())+".");
      throw runtime_error("Too high division.");
    }
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi); 
      if (m_rng->drnd() < prob_div*(1.0-p.coordination/m_rho_max))
      {
        cout << " dividing particle of type " << p.get_type() << endl;
        Particle p_new(m_system->size(), p.get_type(), p.get_radius());
        p_new.x = p.x + m_alpha*m_split_distance*p.get_radius()*p.nx;
        p_new.y = p.y + m_alpha*m_split_distance*p.get_radius()*p.ny;
        p_new.z = p.z + m_alpha*m_split_distance*p.get_radius()*p.nz;
        p_new.set_parent(p.get_flag());
        m_system->apply_periodic(p_new.x,p_new.y,p_new.z);
        
        p.x -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nx;
        p.y -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.ny;
        p.z -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nz;
        m_system->apply_periodic(p.x,p.y,p.z);
        
        p_new.nx = p.nx; p_new.ny = p.ny; p_new.nz = p.nz;
        p_new.vx = p.vx; p_new.vy = p.vy; p_new.vz = p.vz;
        p_new.Nx = p.Nx; p_new.Ny = p.Ny; p_new.Nz = p.Nz;
        p.age = 0.0;     // age of parent is 0
        p_new.age = 0.0; // age of child is 0
        for(list<string>::iterator it_g = p.groups.begin(); it_g != p.groups.end(); it_g++)
          p_new.add_group(*it_g);
        if (p.in_tissue) p_new.in_tissue = true;
        p_new.set_radius(p.get_radius());
        p_new.set_length(p.get_length());
        p_new.set_default_area(p.get_A0());
        p_new.A0 = p.A0;
        if (m_rng->drnd() < m_type_change_prob_1)  // Attempt to change type and group for first child
        {
          if (m_new_type == 0)
            new_type = p.get_type();
          else
            new_type = m_new_type;
          p.set_type(new_type);
          m_system->change_group(p.get_id(),m_old_group,m_new_group);
        }
        // For the polydispersity function: Change radius of second child.
        if (m_new_radius == 0.0)
          new_r = p_new.get_radius(); 
        else 
        {
          if (m_poly == 0.0)
            new_r = m_new_radius;
          else 
            new_r = m_new_radius*(1.0 + m_poly*(m_rng->drnd() - 0.5));
        }
        p_new.set_radius(new_r);
        if (m_rng->drnd() < m_type_change_prob_2)  // Attempt to change type and group for second child
        {
          if (m_new_type == 0)
            new_type = p_new.get_type();
          else
            new_type = m_new_type;
          p_new.set_type(new_type);
          m_system->add_particle(p_new);
          m_system->change_group(p_new.get_id(),m_old_group,m_new_group);
        }
        else
          m_system->add_particle(p_new);
        cout << "old particle: " << p << endl;
        cout << "new particle: " << p << endl;
        cout << "groups of new particle: " << endl;
        for (list<string>::const_iterator it = p_new.groups.begin(); it != p_new.groups.end(); it++)
            cout << format(" %s ") % (*it);
        cout << endl << "groups of old particle: " << endl;
        for (list<string>::const_iterator it = p.groups.begin(); it != p.groups.end(); it++)
            cout << format(" %s ") % (*it);
        //p.print_groups();
        //p_new.print_groups();
      }
    }
    if (!m_system->group_ok(m_group_name))
    {
      cout << "After Divide P: Group info mismatch for group : " << m_group_name << endl;
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
void PopulationDensity::remove(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_death_rate > 0.0)  // Attempt removal only at certain time steps
  { 
    cout << "Handling deaths for group " << m_group_name << endl;
    if (!m_system->group_ok(m_group_name))
    {
      cerr << "Before Remove P: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    vector<int> to_remove;
    double prob_death = m_death_rate*m_freq*m_system->get_integrator_step(); // actual probability of dividing now: rate * (attempt_freq * dt)
    if (prob_death>1.0)
    {
      cerr << "Error: death rate " << prob_death << " is too large for current time step and attempt rate!" << endl;
      throw runtime_error("Too high death.");
    }
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (m_rng->drnd() < prob_death) {
        cout << "particle of type " << p.get_type() << " died " << endl;
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
      m_msg->msg(Messenger::ERROR,"Density population control. No particles left in the system. Please reduce that death rate.");
      throw runtime_error("No particles left in the system.");
    }
    if (!m_system->group_ok(m_group_name))
    {
      cout << "After Remove P: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    m_system->set_force_nlist_rebuild(true);
  }
}
