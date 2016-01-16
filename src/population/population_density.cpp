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
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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
    double prob_div = m_div_rate*m_freq*m_system->get_step(); // actual probability of dividing now: rate * (attempt_freq * dt)
    if (prob_div>1.0)
      {
	cout << "Error: division rate " << prob_div << " is too large for current time step and attempt rate!" << endl;
	cout << "We have: division " << m_div_rate << " frequency " << m_freq << " time step " << m_system->get_step() << endl;
	throw runtime_error("Too high division.");
      }
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi); 
      if (m_rng->drnd() < m_div_rate*(1.0-p.coordination/m_rho_max) )
      {
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
        p.age = 0.0; // age of parent is 0
        p_new.age = 0.0; // age of child is 0: make sure age-dependent potentials are consistent with this
        for(list<string>::iterator it_g = p.groups.begin(); it_g != p.groups.end(); it_g++)
          p_new.groups.push_back(*it_g);
        if (m_rng->drnd() < m_type_change_prob_1)  // Attempt to change type, radius and group for first child
        {
          if (m_new_type == 0)
            new_type = p.get_type();
          else
            new_type = m_new_type;
          p.set_type(new_type);
          if (m_new_radius == 0.0)
            new_r = p.get_radius();
          else
            new_r = m_new_radius;
          p.set_radius(new_r);
          m_system->change_group(p,m_old_group,m_new_group);
        }
        m_system->add_particle(p_new);
        Particle& pr = m_system->get_particle(p_new.get_id());
        if (m_rng->drnd() < m_type_change_prob_2)  // Attempt to change type, radius and group for second child
        {
          if (m_new_type == 0)
            new_type = pr.get_type();
          else
            new_type = m_new_type;
          pr.set_type(new_type);
          if (m_new_radius == 0.0)
            new_r = pr.get_radius();
          else
            new_r = m_new_radius;
          pr.set_radius(new_r);
          m_system->change_group(pr,m_old_group,m_new_group);
        }
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
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before Remove P: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    vector<int> to_remove;
    double prob_death = m_death_rate*m_freq*m_system->get_step(); // actual probability of dividing now: rate * (attempt_freq * dt)
    if (prob_death>1.0)
      {
	cout << "Error: death rate " << prob_death << " is too large for current time step and attempt rate!" << endl;
	throw runtime_error("Too high death.");
      }
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (m_rng->drnd() < prob_death)
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