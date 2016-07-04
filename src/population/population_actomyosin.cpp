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
 * \file population_actomyosin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 08-Jun-2016
 * \brief Implementation of PopulationActomyosin class.
 */ 

#include "population_actomyosin.hpp"


/*! This function controls attachement, that is transition from the
 *  detached state (D) to the attached state (A). For simiplicity, we
 *  assume that this happens with proability set by parameter attachment_prob
 *
 *  \param t current time step
 *  
*/
void PopulationActomyosin::divide(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_attach_prob > 0.0)  // Attempt D to A transition only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before actomyosin division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int N = m_system->size();
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(i); 
      if (pi.get_type() == m_type_d)
        if (m_rng->drnd() < m_attach_prob)
          pi.set_type(m_type_a);
    }
  }
}

/*! This function handles detachement. If myosin particle is more than \f$ r_e \f$ away from an
 *  actin bead, it will transition A to D with proability set by paramter detachement_prob.
 *  If it is within \f$ r_e \f$ from an actin bead the detachement proability will be
 *  \f$ a\exp(\lambda f) \f$ where \f$ a \f$ is the bare proability set by "detachement_prob",
 *  \f$ \lambda \f$ is the paramter set in the input file and \f$ f \f$ is the magnitude of the 
 *  force acting on the myosin bead.
 * 
 *  \param t current time step
 * 
*/
void PopulationActomyosin::remove(int t)
{
  if (!m_has_nlist)
    throw runtime_error("Actomyosin population control requires neighbour list to be defined.");
  if (m_freq > 0 && t % m_freq == 0 && m_detach_prob > 0.0)  // Attempt of the A to D transition at certain time steps
  { 
    int N = m_system->size();
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(i);
      if (pi.get_type() == m_type_a)
      {
        vector<int>& neigh = m_nlist->get_neighbours(i);
        double min_dist = 1e10;
        for (unsigned int j = 0; j < neigh.size(); j++)
        {
          Particle& pj = m_system->get_particle(neigh[j]);
          if (pj.get_type() == m_type_actin)
          {
            double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
            m_system->apply_periodic(dx,dy,dz);
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (r < min_dist)
              min_dist = r;
          }
        }
        double prob = m_detach_prob;
        if (min_dist < m_re)
        {
          double f = sqrt(pi.fx*pi.fx + pi.fy*pi.fy + pi.fz*pi.fz);
          prob *= exp(m_lambda*f);
        }
        if (m_rng->drnd() < prob)
          pi.set_type(m_type_d);
      }
    }
  }
}
