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
 * \file population_elongation.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Jun-2015
 * \brief Implementation of PopulationElongation class.
 */ 

#include "population_elongation.hpp"


/*! Grow (rescale) rod length by a given amount.
 * 
 *  \param t current time step
 *  
*/
void PopulationElongation::elongate(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && t < m_rescale_steps && m_rescale != 1.0) 
  { 
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    for (int i = 0; i < N; i++)
    {
      int pi = particles[i];
      Particle& p = m_system->get_particle(pi); 
      p.scale_length(m_scale);
    }
    m_system->set_force_nlist_rebuild(true);
    m_system->set_nlist_rescale(m_scale);
  }
}