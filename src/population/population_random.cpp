/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file population_random.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Implementation of PopulationRandom class.
 */ 

#include "population_random.hpp"


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
void PopulationRandom::divide(int t)
{
  if (t % m_freq == 0)  // Attempt division only at certain time steps
  { 
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    bool periodic = m_system->get_periodic();
    BoxPtr box = m_system->get_box();
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (m_rng->drnd() < p.age/m_div_rate)
      {
        Particle p_new(m_system->size(), p.get_type(), p.get_radius());
        p_new.x = p.x + 0.5*p.get_radius()*p.nx;
        p_new.y = p.y + 0.5*p.get_radius()*p.ny;
        p_new.z = p.z + 0.5*p.get_radius()*p.nz;
        if (periodic)
        {
          if (p_new.x > box->xhi) p_new.x -= box->Lx;
          else if (p_new.x < box->xlo) p_new.x += box->Lx;
          if (p_new.y > box->yhi) p_new.y -= box->Ly;
          else if (p_new.y < box->ylo) p_new.y += box->Ly;
          if (p_new.z > box->zhi) p_new.z -= box->Lz;
          else if (p_new.z < box->zlo) p_new.z += box->Lz;
        }
        p.x -= 0.5*p.get_radius()*p.nx;
        p.y -= 0.5*p.get_radius()*p.ny;
        p.z -= 0.5*p.get_radius()*p.nz;
        if (periodic)
        {
          if (p.x > box->xhi) p.x -= box->Lx;
          else if (p.x < box->xlo) p.x += box->Lx;
          if (p.y > box->yhi) p.y -= box->Ly;
          else if (p.y < box->ylo) p.y += box->Ly;
          if (p.z > box->zhi) p.z -= box->Lz;
          else if (p.z < box->zlo) p.z += box->Lz;
        }
        p_new.nx = p.nx; p_new.ny = p.ny; p_new.nz = p.nz;
        p_new.vx = p.vx; p_new.vy = p.vy; p_new.vz = p.vz;
        p.age = 0.0;
        p_new.age = 0.0;
        for(list<string>::iterator it_g = p.groups.begin(); it_g != p.groups.end(); it_g++)
          p_new.groups.push_back(*it_g);
        m_system->add_particle(p_new);
      }
    }
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
void PopulationRandom::remove(int t)
{
  if (t % m_freq == 0)  // Attempt removal only at certain time steps
  { 
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi);
      if (m_rng->drnd() < p.age/m_death_rate)
        m_system->remove_particle(p);
      if (m_system->size() == 0)
      {
        m_msg->msg(Messenger::ERROR,"Random population control. No particles left in the system. Please reduce that death rate.");
        throw runtime_error("No particles left in the system.");
      }
    }
  }
}