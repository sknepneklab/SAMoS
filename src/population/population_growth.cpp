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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file population_growth.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Jun-2015
 * \brief Implementation of PopulationGrowth class.
 */ 

#include "population_growth.hpp"


/*! Grow (rescale) particle radius by a given amount.
 * 
 *  \param t current time step
 *  
*/
void PopulationGrowth::grow(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && t < m_rescale_steps && m_rescale != 1.0) 
  { 
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi); 
      p.scale_radius(m_scale);
    }
    m_system->set_force_nlist_rebuild(true);
    m_system->set_nlist_rescale(m_scale);
  }
}